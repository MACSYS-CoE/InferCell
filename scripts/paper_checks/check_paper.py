"""Mechanical checks on the manuscript (the /write-paper skill).

    python3 scripts/paper_checks/check_paper.py              # build, then check
    python3 scripts/paper_checks/check_paper.py --no-build   # check the last build's log
    python3 scripts/paper_checks/check_paper.py --final      # also fail on stubs, DRAFT markers, unused numbers, uncited bib entries
    python3 scripts/paper_checks/check_paper.py --paper docs/paper   # paper directory, if not ./paper
    python3 scripts/paper_checks/check_paper.py --selftest   # prove the static checks catch breakage

The manuscript is paper/main.tex with every \\input/\\include expanded (sections/*.tex,
tables/*.tex). Checks:
  1. `latexmk -pdf main.tex` exits 0 (unless --no-build);
  2. main.log has no undefined reference or citation, no multiply defined label, and no
     overfull hbox over 5 pt;
  3. every \\includegraphics target exists, searched in the \\graphicspath directories;
  4. every figure, table and appendix label (fig: tab: app:) is \\cref'd from the text;
  5. every \\pn{key} used is defined in paper/numbers.tex;
  6. paper/numbers.tex is current: the project's numbers generator (scripts/paper_numbers.py,
     or scripts/make_paper_numbers.py), run with `--out <tmp>`, regenerates it byte for byte.
     A missing generator is a failure: every number the paper quotes is a \\pn{} macro.
     A generator without --out is a note. --python (or $PAPER_PYTHON) picks its interpreter.
     If scripts/paper_tables.py exists, paper/tables/ is checked the same way (`<tmpdir>`).
With --final, also: no \\stub{} or \\DRAFT{} remains, every number in numbers.tex is used,
and every references.bib entry is cited. Without --final those are printed as notes.
If the project README.md has a figures table, a figure it does not name is a note.
Exits 1 on any failure. Standard library only; imports helpers from check_prose.py.
"""
import argparse
import os
import re
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from check_prose import body, find_paper, flatten, read, strip_comments  # noqa: E402

OVERFULL_PT = 5.0
REF_PREFIXES = ("fig:", "tab:", "app:")
GRAPHIC_EXT = ("", ".pdf", ".png", ".jpg", ".jpeg", ".eps")
GENERATORS = ("paper_numbers.py", "make_paper_numbers.py")


def build(paper):
    r = subprocess.run(["latexmk", "-pdf", "-interaction=nonstopmode", "-halt-on-error", "main.tex"],
                       cwd=paper, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    return r.returncode


def log_problems(log):
    out = []
    for m in re.finditer(r"Warning: (Reference|Citation) `([^']+)' on page", log):
        out.append("undefined %s %s" % (m.group(1).lower(), m.group(2)))
    for m in re.finditer(r"Label `([^']+)' multiply defined", log):
        out.append("label %s is multiply defined" % m.group(1))
    if re.search(r"There were undefined (references|citations)", log) and not out:
        out.append("undefined references or citations (see main.log)")
    for m in re.finditer(r"Overfull \\hbox \(([\d.]+)pt too wide\)[^\n]*lines? (\d+)", log):
        if float(m.group(1)) > OVERFULL_PT:
            out.append("overfull hbox %spt at line %s" % (m.group(1), m.group(2)))
    return sorted(set(out), key=out.index)


def graphics_dirs(tex, paper):
    m = re.search(r"\\graphicspath\{((?:\{[^}]*\})+)\}", tex)
    dirs = re.findall(r"\{([^}]*)\}", m.group(1)) if m else []
    return [os.path.normpath(os.path.join(paper, d)) for d in dirs] + [paper]


def graphics_problems(doc, dirs, readme):
    fails, notes = [], []
    for name in re.findall(r"\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}", doc):
        if os.path.basename(name).startswith("orcid"):
            continue
        hit = next((os.path.join(d, name + e) for d in dirs for e in GRAPHIC_EXT
                    if os.path.isfile(os.path.join(d, name + e))), None)
        if not hit:
            fails.append("missing figure %s (searched %s)" % (name, ", ".join(os.path.relpath(d) for d in dirs)))
        elif readme and "figures/" in readme and os.path.basename(hit) not in readme:
            notes.append("figure %s is not named in README.md" % os.path.basename(hit))
    return fails, notes


def ref_problems(doc):
    labels = [l for l in re.findall(r"\\label\{([^}]+)\}", doc) if l.startswith(REF_PREFIXES)]
    refs = set()
    for group in re.findall(r"\\(?:[cC]ref|[cC]refrange|ref|autoref|eqref)\{([^}]+)\}", doc):
        refs.update(x.strip() for x in group.split(","))
    return ["label %s is never \\cref'd" % l for l in labels if l not in refs]


def defined_numbers(numbers):
    return re.findall(r"^\\pndef\{([^}]+)\}", numbers, flags=re.M)


def cited_keys(doc):
    cited = set()
    for group in re.findall(r"\\(?:no)?cite[a-zA-Z]*\*?(?:\[[^\]]*\]){0,2}\{([^}]*)\}", doc):
        cited.update(x.strip() for x in group.split(","))
    return cited


def uncited(doc, bib):
    cited = cited_keys(doc)
    if "*" in cited:
        return []
    keys = re.findall(r"^\s*@\w+\s*\{\s*([^,\s]+)\s*,", bib, flags=re.M)
    return ["references.bib entry %s is never cited" % k for k in keys if k not in cited]


def static(tex, numbers, bib, dirs, readme, final):
    """(failures, notes) of every check that reads only the sources. tex is flattened."""
    doc = body(tex)
    fails, notes = graphics_problems(doc, dirs, readme)
    fails += ref_problems(doc)
    defined = set(defined_numbers(numbers))
    used = set(re.findall(r"\\pn\{([^}]+)\}", doc))
    fails += ["undefined number \\pn{%s}" % k for k in sorted(used - defined)]
    unused = sorted(defined - used)
    stubs, drafts = len(re.findall(r"\\stub\{", tex)), len(re.findall(r"\\DRAFT\{", tex))
    markers = len(re.findall(r"\[CITE[ -](NEEDED|PROPOSED)", tex))
    late = ["%d \\stub{}, %d \\DRAFT{} and %d [CITE ...] markers remain" % (stubs, drafts, markers)] \
        if stubs or drafts or markers else []
    late += ["unused number %s" % k for k in unused]
    late += uncited(doc, bib)
    if final:
        fails += late
    else:
        notes += late
    return fails, notes


def generator(root):
    return next((os.path.join(root, "scripts", g) for g in GENERATORS
                 if os.path.isfile(os.path.join(root, "scripts", g))), None)


def stale_numbers(root, paper, python):
    """(failures, notes) for numbers.tex against its generator."""
    gen = generator(root)
    if not gen:
        return ["no numbers generator: write scripts/paper_numbers.py (see the write-paper skill's "
                "assets/paper_numbers_template.py)"], []
    if "--out" not in read(gen):
        return [], ["%s has no --out option, so numbers.tex staleness is not checked "
                    "(add --out as in the write-paper template)" % os.path.basename(gen)]
    fd, tmp = tempfile.mkstemp(suffix=".tex")
    os.close(fd)
    try:
        r = subprocess.run([python, gen, "--out", tmp], cwd=root,
                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        if r.returncode:
            hint = " (wrong interpreter? pass --python or set PAPER_PYTHON)" if b"ModuleNotFoundError" in r.stdout else ""
            return ["%s --out failed%s: %s" % (os.path.basename(gen), hint, r.stdout.decode()[-300:])], []
        here = os.path.join(paper, "numbers.tex")
        if not os.path.isfile(here):
            return ["no paper/numbers.tex: run scripts/%s" % os.path.basename(gen)], []
        return ([] if read(tmp) == read(here) else
                ["paper/numbers.tex is stale: rerun scripts/%s" % os.path.basename(gen)]), []
    finally:
        os.remove(tmp)


def stale_tables(root, paper, python):
    gen = os.path.join(root, "scripts", "paper_tables.py")
    if not os.path.isfile(gen):
        return []
    tmp = tempfile.mkdtemp()
    try:
        r = subprocess.run([python, gen, tmp], cwd=root, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        if r.returncode:
            return ["paper_tables.py failed: %s" % r.stdout.decode()[-300:]]
        here = os.path.join(paper, "tables")
        names = set(os.listdir(tmp)) | set(os.listdir(here) if os.path.isdir(here) else [])
        return ["paper/tables/%s is stale or missing: rerun python3 scripts/paper_tables.py" % n
                for n in sorted(names)
                if not os.path.isfile(os.path.join(here, n)) or not os.path.isfile(os.path.join(tmp, n))
                or read(os.path.join(here, n)) != read(os.path.join(tmp, n))]
    finally:
        for n in os.listdir(tmp):
            os.remove(os.path.join(tmp, n))
        os.rmdir(tmp)


def selftest():
    """Break a minimal document in each way the static checks guard, and expect a failure."""
    figdir = tempfile.mkdtemp()
    open(os.path.join(figdir, "a.pdf"), "w").close()
    numbers = "\\pndef{x}{1}%\n"
    bib = "@article{k1,\n title={t}}\n"
    good = (r"\begin{document}\begin{figure}\includegraphics{a}\label{fig:a}\end{figure}"
            r"See \cref{fig:a}, \pn{x} \citep{k1}.\end{document}")
    cases = {
        "clean": (good, False, 0),
        "broken cref": (good.replace(r"\cref{fig:a}", r"\cref{fig:b}"), False, 1),
        "missing figure": (good.replace("{a}", "{z}"), False, 1),
        "undefined number": (good.replace(r"\pn{x}", r"\pn{y}"), False, 1),
        "stub, not final": (good.replace("See", r"\stub{later} See"), False, 0),
        "stub under --final": (good.replace("See", r"\stub{later} See"), True, 1),
        "cite marker --final": (good.replace("See", "[CITE NEEDED: x] See"), True, 1),
        "uncited bib --final": (good.replace(r"\citep{k1}", ""), True, 1),
    }
    ok = True
    for name, (tex, final, want) in cases.items():
        fails, _ = static(tex, numbers, bib, [figdir], "", final)
        got = 1 if fails else 0
        ok &= got == want
        print("%s  %-24s %s" % ("ok  " if got == want else "FAIL", name, fails or "no failures"))
    log = ("Overfull \\hbox (12.0pt too wide) in paragraph at lines 10--12\n"
           "LaTeX Warning: Reference `fig:q' on page 1 undefined"
           "\nPackage natbib Warning: Citation `x2000' on page 2 undefined on input l"
           "\nLaTeX Warning: Label `fig:a' multiply defined.")
    lp = log_problems(log)
    ok &= len(lp) == 4
    print("%s  %-24s %s" % ("ok  " if len(lp) == 4 else "FAIL", "log problems", lp))
    d = tempfile.mkdtemp()
    os.makedirs(os.path.join(d, "sections"))
    with open(os.path.join(d, "main.tex"), "w") as f:
        f.write("\\input{numbers}\\begin{document}\\input{sections/results}\\end{document}")
    with open(os.path.join(d, "sections", "results.tex"), "w") as f:
        f.write("\\begin{figure}\\includegraphics{a}\\label{fig:b}\\end{figure} \\pn{y} % \\cref{fig:b}\n")
    fails, _ = static(flatten(os.path.join(d, "main.tex")), numbers, "", [figdir], "", False)
    want = ["label fig:b is never \\cref'd", "undefined number \\pn{y}"]
    ok &= fails == want
    print("%s  %-24s %s" % ("ok  " if fails == want else "FAIL", "\\input expanded", fails))
    return 0 if ok else 1


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--paper", help="paper directory (default: ../paper relative to this script, or ./paper)")
    ap.add_argument("--no-build", action="store_true")
    ap.add_argument("--python", default=os.environ.get("PAPER_PYTHON", sys.executable),
                    help="interpreter for the numbers/tables generators (default: $PAPER_PYTHON or this one)")
    ap.add_argument("--final", action="store_true")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        sys.exit(selftest())
    paper = find_paper(a.paper)
    root = os.path.dirname(paper) if os.path.basename(os.path.dirname(paper)) != "docs" \
        else os.path.dirname(os.path.dirname(paper))
    fails, notes0 = stale_numbers(root, paper, a.python)
    fails += stale_tables(root, paper, a.python)
    if not a.no_build:
        rc = build(paper)
        if rc:
            fails.append("latexmk exited %d (see %s/main.log)" % (rc, os.path.relpath(paper)))
    log = os.path.join(paper, "main.log")
    if os.path.exists(log):
        fails += log_problems(read(log))
    else:
        fails.append("no main.log: build first")
    tex = flatten(os.path.join(paper, "main.tex"))
    numbers = read(os.path.join(paper, "numbers.tex")) if os.path.isfile(os.path.join(paper, "numbers.tex")) else ""
    bibs = re.findall(r"\\bibliography\{([^}]+)\}", tex)
    bib = "".join(read(os.path.join(paper, b.strip() + ("" if b.strip().endswith(".bib") else ".bib")))
                  for g in bibs for b in g.split(",")
                  if os.path.isfile(os.path.join(paper, b.strip() + ("" if b.strip().endswith(".bib") else ".bib"))))
    readme = read(os.path.join(root, "README.md")) if os.path.isfile(os.path.join(root, "README.md")) else ""
    f2, notes = static(tex, numbers, bib, graphics_dirs(tex, paper), readme, a.final)
    fails += f2
    notes = notes0 + notes
    for n in notes:
        print("note: " + n)
    for f in fails:
        print("FAIL: " + f)
    print("check_paper: %s" % ("FAIL (%d)" % len(fails) if fails else "ok"))
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()

"""Prose metrics for the manuscript, per section (the /write-paper skill).

    python3 scripts/paper_checks/check_prose.py                    # report per section
    python3 scripts/paper_checks/check_prose.py --final            # also fail when a budget is missed
    python3 scripts/paper_checks/check_prose.py --dump-text Results  # plain text of one section (for Pangram)
    python3 scripts/paper_checks/check_prose.py --dump-text list   # list the section titles
    python3 scripts/paper_checks/check_prose.py --paper docs/paper # paper directory, if not ./paper
    python3 scripts/paper_checks/check_prose.py --selftest

The manuscript is paper/main.tex with every \\input/\\include expanded (numbers.tex is not
expanded). A rough LaTeX stripper turns it into plain sentences: it drops comments,
floats, display maths, citations, \\stub{} text and the preamble; inline maths becomes
one word ("X"); \\pn{} becomes "N"; \\cref{} becomes "Ref"; text commands (\\emph,
\\textbf, ...) keep their argument. Sentences split at . ! ? before a capital letter or
at the end of a paragraph.

Sections: the abstract; every numbered \\section before \\appendix is MAIN TEXT;
\\section* (author summary, acknowledgements, data availability) is unnumbered and not
budgeted; every \\section after \\appendix (or inside \\begin{appendices}) is an appendix.

Per section it reports words, sentences, the mean sentence length and the sentences over
LONG words, and lists bare numbers (digits typed into prose instead of \\pn{} macros).

Configuration comes from paper/story.md (written at Stage 1 of /write-paper):
  * the "Not" column of the term table (a markdown table whose header row has the cells
    "Use" and "Not") gives the avoided terms. Parenthetical notes are ignored, and cells
    are split on commas and " or ";
  * a "## Budgets" section of `- key: value` lines overrides the defaults below.
    Keys: main_words, mean_max, long, long_max, abstract_max (integers);
    term_sections (comma-separated section titles, or "main"); allowed_numbers
    (comma-separated tokens that may appear as bare numbers in Results and Abstract).

--final fails if: the main text is over main_words; any main-text section has a mean
sentence over mean_max; more than long_max main-text sentences are over `long` words;
the abstract is over abstract_max; an avoided term appears in a term_sections section;
or Results or the Abstract has a bare number with a decimal point or at least three
digits that is not in allowed_numbers (a measured quantity typed by hand). Without
--final, the same findings are printed as notes and the exit code is 0.

check_paper.py and check_invariants.py import the helpers here. Standard library only.
"""
import argparse
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def project_root(start=HERE):
    """The nearest directory above `start` that has paper/main.tex or docs/paper/main.tex.
    The scripts live in <project>/scripts/paper_checks/, so this is usually two levels up."""
    d = start
    while True:
        if any(os.path.isfile(os.path.join(d, sub, "main.tex")) for sub in ("paper", os.path.join("docs", "paper"))):
            return d
        parent = os.path.dirname(d)
        if parent == d:
            return os.path.normpath(os.path.join(start, ".."))
        d = parent


ROOT = project_root()
DEFAULTS = {"main_words": 6000, "mean_max": 20, "long": 35, "long_max": 5, "abstract_max": 300,
            "term_sections": "Methods, Results, Discussion", "allowed_numbers": ""}


# ---------------------------------------------------------------- shared helpers

def read(path):
    with open(path, encoding="utf-8", errors="replace") as f:  # TeX logs are not always UTF-8
        return f.read()


def strip_comments(tex):
    return re.sub(r"(?<!\\)%.*", "", tex)


def flatten(path, skip=("numbers",)):
    """main.tex with every \\input{x} / \\include{x} replaced by the file's text, recursively.
    Paths resolve against the directory of the top file, as LaTeX does. Files whose stem is in
    `skip` are left as they are; a missing file is left as the \\input command."""
    base = os.path.dirname(os.path.abspath(path))

    def expand(tex, level):
        def sub(m):
            name = m.group(2).strip()
            stem = os.path.splitext(os.path.basename(name))[0]
            if stem in skip or level > 10:
                return m.group(0)
            for cand in (name, name + ".tex"):
                p = os.path.join(base, cand)
                if os.path.isfile(p):
                    return "\n" + expand(strip_comments(read(p)), level + 1) + "\n"
            return m.group(0)
        return re.sub(r"\\(input|include)\{([^}]+)\}", sub, tex)

    return expand(strip_comments(read(path)), 0)


def balanced(s, i):
    """Index just past the {...} group that opens at s[i] == '{'."""
    depth = 0
    for j in range(i, len(s)):
        if s[j] == "{":
            depth += 1
        elif s[j] == "}":
            depth -= 1
            if depth == 0:
                return j + 1
    return len(s)


def drop_command(s, name, repl=""):
    """Replace every \\name{...} (balanced braces) with repl."""
    out, i, pat = [], 0, "\\" + name + "{"
    while True:
        k = s.find(pat, i)
        if k < 0:
            return "".join(out) + s[i:]
        out.append(s[i:k] + repl)
        i = balanced(s, k + len(pat) - 1)


def body(tex):
    """The text between \\begin{document} and \\end{document} (or all of it)."""
    tex = tex.split(r"\begin{document}", 1)[-1]
    return tex.split(r"\end{document}", 1)[0]


SECTION_RE = re.compile(r"\\(section|appendix)(\*?)(?:\s*\{)?")


def split_sections(tex):
    """[(title, kind, raw tex)] of a document body, kind in abstract / main / unnumbered /
    appendix / front. Raw LaTeX is kept, so callers can strip it or scan it for keys."""
    out = []
    ab = re.search(r"\\begin\{abstract\}(.*?)\\end\{abstract\}", tex, flags=re.S)
    if ab:
        out.append(("Abstract", "abstract", ab.group(1)))
        tex = tex[:ab.start()] + tex[ab.end():]
    tex = re.sub(r"\\begin\{thebibliography\}.*?\\end\{thebibliography\}", " ", tex, flags=re.S)
    tex = re.sub(r"\\(bibliography|bibliographystyle)\{[^}]*\}|\\printbibliography", " ", tex)
    tex = tex.replace("\\begin{appendices}", "\\appendix").replace("\\end{appendices}", " ")
    title, kind, appendix, start = "Front matter", "front", False, 0
    for m in SECTION_RE.finditer(tex):
        if m.group(1) == "section" and not m.group(0).endswith("{"):
            continue  # \sectionmark and the like
        chunk = tex[start:m.start()]
        if title is not None and chunk.strip():
            out.append((title, kind, chunk))
        if m.group(1) == "appendix":
            appendix, title, start = True, None, m.end()
            continue
        close = balanced(tex, m.end() - 1)
        title = " ".join(tex[m.end():close - 1].split())
        kind = "appendix" if appendix else ("unnumbered" if m.group(2) else "main")
        start = close
    chunk = tex[start:]
    if title is not None and chunk.strip():
        out.append((title, kind, chunk))
    return [s for s in out if s[1] != "front" or sentences(strip(s[2]))]


def strip(tex):
    """Plain text from LaTeX prose: floats, display maths and citations go; inline maths
    becomes X, \\pn{} N, \\cref{} Ref."""
    tex = strip_comments(tex)
    for env in ("figure", "figure\\*", "table", "table\\*", "equation", "equation\\*", "align",
                "align\\*", "gather", "gather\\*", "multline", "multline\\*", "eqnarray",
                "tabular", "tabularx", "tikzpicture", "algorithm", "verbatim", "lstlisting"):
        tex = re.sub(r"\\begin\{%s\}.*?\\end\{%s\}" % (env, env), " ", tex, flags=re.S)
    tex = re.sub(r"\\\[.*?\\\]", " ", tex, flags=re.S)
    tex = re.sub(r"\$\$.*?\$\$", " ", tex, flags=re.S)
    for name in ("stub", "DRAFT", "label", "subsection", "subsection*", "subsubsection",
                 "subsubsection*", "paragraph", "paragraph*", "keywords", "footnote", "footnotetext",
                 "vspace", "hspace", "FloatBarrier"):
        tex = drop_command(tex, name)
    tex = drop_command(tex, "pn", "N")
    for name in ("cref", "Cref", "ref", "eqref", "autoref", "crefrange", "Crefrange"):
        tex = drop_command(tex, name, "Ref")
    tex = re.sub(r"\\cite[a-zA-Z]*\*?(\[[^\]]*\]){0,2}\{[^}]*\}", "", tex)
    tex = re.sub(r"\\\(.*?\\\)", "X", tex, flags=re.S)
    tex = re.sub(r"(?<!\\)\$\s*\\(to|rightarrow|longrightarrow)\s*\$", "to", tex)
    tex = re.sub(r"(?<!\\)\$\s*\\times\s*\$", "x", tex)
    tex = re.sub(r"(?<!\\)\$[^$]*\$", "X", tex)
    tex = re.sub(r"\\item\b(\[[^\]]*\])?", "\n\n", tex)  # each list item is its own paragraph
    tex = re.sub(r"\\(begin|end)\{(itemize|enumerate|description|center|quote)\}", "\n\n", tex)
    for _ in range(3):  # nested text commands
        tex = re.sub(r"\\(emph|textbf|textit|textrm|textsf|textsc|code|texttt|mbox|text|"
                     r"textsuperscript|textsubscript|url|href\{[^}]*\})\{([^{}]*)\}", r"\2", tex)
    return tex


def paragraphs(text):
    """Plain paragraphs: leftover commands and braces removed, whitespace collapsed."""
    text = re.sub(r"\\[a-zA-Z]+\*?(\[[^\]]*\])?", " ", text)
    text = text.replace("~", " ").replace("\\,", " ").replace("--", "–")
    text = re.sub(r"\\[%&#_$]", lambda m: m.group(0)[1], text)
    text = re.sub(r"[{}]", "", text).replace("``", "\"").replace("''", "\"")
    paras = [re.sub(r"\s+([.,;:)])", r"\1", re.sub(r"\(\s*\)", "", " ".join(p.split())))
             for p in re.split(r"\n\s*\n", text)]
    return [p for p in paras if re.search(r"[A-Za-z]", p)]


ABBREV = r"\b(e\.g|i\.e|et al|cf|Fig|Figs|Eq|Eqs|Sec|Ref|vs|approx|resp)\."


def sentences(text):
    out = []
    for para in paragraphs(text):
        para = re.sub(ABBREV, r"\1", para)
        for s in re.split(r"(?<=[.!?])\s+(?=[A-Z(])", para):
            words = re.findall(r"[A-Za-z0-9][\w'’\-\.]*", s)
            if words:
                out.append((s, len(words)))
    return out


def find_paper(arg):
    """The paper directory: --paper, else ROOT/paper, else ./paper."""
    for cand in ([arg] if arg else []) + [os.path.join(ROOT, "paper"), os.path.join(ROOT, "docs", "paper"), "paper"]:
        if cand and os.path.isfile(os.path.join(cand, "main.tex")):
            return os.path.abspath(cand)
    sys.exit("no main.tex found: pass --paper <dir>")


# ---------------------------------------------------------------- story.md configuration

def story_config(paper_dir):
    """(avoided-term patterns, budgets) from paper/story.md; defaults when absent."""
    cfg = dict(DEFAULTS)
    path = os.path.join(paper_dir, "story.md")
    if not os.path.isfile(path):
        return [], cfg
    return parse_story(read(path), cfg)


def parse_story(md, cfg):
    avoid, in_terms, not_col = [], False, None
    for line in md.splitlines():
        cells = [c.strip() for c in line.strip().strip("|").split("|")] if line.strip().startswith("|") else None
        if cells and "Use" in cells and "Not" in cells:
            in_terms, not_col = True, cells.index("Not")
            continue
        if in_terms and cells:
            if set("".join(cells)) <= set("-: "):
                continue
            if not_col < len(cells):
                cell = re.sub(r"\([^)]*\)", "", cells[not_col]).replace("`", "").replace("<", "").replace(">", "")
                for term in re.split(r",|\bor\b|;", cell):
                    term = term.strip().strip("*").strip()
                    if term and term.lower() not in ("—", "-", "n/a"):
                        avoid.append(term_pattern(term))
            continue
        in_terms = False
    budgets = re.search(r"^#+\s*Budgets?\s*$(.*?)(?=^#|\Z)", md, flags=re.M | re.S | re.I)
    if budgets:
        for m in re.finditer(r"^\s*[-*]\s*`?(\w+)`?\s*:\s*(.+?)\s*$", budgets.group(1), flags=re.M):
            key, val = m.group(1), m.group(2).strip().strip("`")
            if key in cfg:
                cfg[key] = int(val) if isinstance(DEFAULTS[key], int) and re.fullmatch(r"\d+", val) else val
    return avoid, cfg


def term_pattern(term):
    """A case-insensitive pattern for a term and its plural; hyphen/space variants match."""
    words = [re.escape(w) for w in re.split(r"[\s\-]+", term) if w]
    core = r"[\s\-]?".join(words)
    return r"\b%s(s|es)?\b" % core


# ---------------------------------------------------------------- report

def bare_numbers(text, allowed):
    """Pure-number tokens in plain text, excluding years and allowed tokens."""
    out = []
    for m in re.finditer(r"(?<![\w.\-])\d[\d,{}]*(\.\d+)?(?![\w\-])", text):
        tok = m.group(0).rstrip(",")
        if re.fullmatch(r"(19|20)\d\d", tok) or tok in allowed:
            continue
        out.append(tok)
    return out


def analyse(tex, avoid, cfg):
    """(rows, findings) for a flattened document. findings are messages that fail --final."""
    secs = split_sections(body(tex))
    allowed = {t.strip() for t in re.sub(r"\([^)]*\)", "", str(cfg["allowed_numbers"])).split(",") if t.strip()}
    scope = str(cfg["term_sections"])
    term_secs = None if scope.strip().lower() == "main" else {t.strip().lower() for t in scope.split(",")}
    rows, findings, main_words, main_long = [], [], 0, []
    for title, kind, raw in secs:
        plain = strip(raw)
        ss = sentences(plain)
        words = sum(n for _, n in ss)
        long_ = [s for s, n in ss if n > cfg["long"]]
        mean = words / len(ss) if ss else 0.0
        nums = bare_numbers(" ".join(paragraphs(plain)), allowed)
        rows.append((title, kind, words, len(ss), mean, long_, nums))
        if kind == "main":
            main_words += words
            main_long += long_
            if mean > cfg["mean_max"]:
                findings.append("%s: mean sentence %.1f words > %d" % (title, mean, cfg["mean_max"]))
        if kind == "abstract" and words > cfg["abstract_max"]:
            findings.append("abstract %d words > %d" % (words, cfg["abstract_max"]))
        in_scope = (kind == "main") if term_secs is None else (title.lower() in term_secs)
        if in_scope:
            for pat in avoid:
                for m in re.finditer(pat, plain, flags=re.I):
                    findings.append("%s uses avoided term '%s'" % (title, m.group(0)))
        if title.lower().startswith("results") or kind == "abstract":
            for tok in nums:
                if "." in tok or len(re.sub(r"\D", "", tok)) >= 3:
                    findings.append("%s has a typed number '%s' (use \\pn{} or add it to allowed_numbers)"
                                    % (title, tok))
    if main_words > cfg["main_words"]:
        findings.append("main text %d words > %d" % (main_words, cfg["main_words"]))
    if len(main_long) > cfg["long_max"]:
        findings.append("%d main-text sentences over %d words > %d" % (len(main_long), cfg["long"], cfg["long_max"]))
    return rows, findings, main_words, main_long


def report(tex, avoid, cfg):
    rows, findings, main_words, main_long = analyse(tex, avoid, cfg)
    print("%-36s %-10s %6s %5s %6s %5s %5s" % ("section", "kind", "words", "sent", "mean", ">%d" % cfg["long"], "nums"))
    for title, kind, words, n, mean, long_, nums in rows:
        print("%-36s %-10s %6d %5d %6.1f %5d %5d" % (title[:36], kind, words, n, mean, len(long_), len(nums)))
        if kind == "main":
            for s in long_:
                print("    long: " + s[:150])
        if nums:
            print("    bare numbers: " + ", ".join(nums[:20]) + (" ..." if len(nums) > 20 else ""))
    print("main text: %d words (budget %d), %d sentences over %d words (budget %d)"
          % (main_words, cfg["main_words"], len(main_long), cfg["long"], cfg["long_max"]))
    return findings


def dump_text(tex, which):
    secs = split_sections(body(tex))
    if which.lower() == "list":
        for title, kind, raw in secs:
            words = sum(n for _, n in sentences(strip(raw)))
            print("%-10s %6d  %s" % (kind, words, title))
        return 0
    hits = [s for s in secs if s[0].lower() == which.lower()] or \
           [s for s in secs if s[0].lower().startswith(which.lower())]
    if not hits:
        print("no section titled %r; try --dump-text list" % which, file=sys.stderr)
        return 1
    for title, kind, raw in hits:
        print("\n\n".join(paragraphs(strip(raw))))
    return 0


# ---------------------------------------------------------------- selftest

def selftest():
    import tempfile
    ok = []
    d = tempfile.mkdtemp()
    os.makedirs(os.path.join(d, "sections"))
    with open(os.path.join(d, "main.tex"), "w") as f:
        f.write("\\documentclass{article}\\input{numbers}\\begin{document}\\begin{abstract}One two three.\\end{abstract}\n"
                "\\section*{Author summary}Plain words here.\n\\input{sections/results}\n"
                "\\section{Discussion}Fine.\n\\begin{appendices}\n\\section{Extra}Appendix text 12.5 here.\n\\end{appendices}\n"
                "\\bibliography{references}\\end{document}\n")
    with open(os.path.join(d, "sections", "results.tex"), "w") as f:
        f.write("\\section{Results}\\label{sec:results}This speedup is \\pn{x} (see \\cref{fig:a}). "
                "We measured 3.27 per cent.\\stub{ignored text}\n\n"
                "Second paragraph $a_0$ here, e.g. this one. % a comment 999\n")
    with open(os.path.join(d, "story.md"), "w") as f:
        f.write("# Story\n\n## Terms\n\n| Use | Not |\n|---|---|\n| ratio R | speedup, gain (in Results) |\n"
                "| throughput | performance or speed |\n\n## Budgets\n\n- abstract_max: 2\n- mean_max: 30\n")
    tex = flatten(os.path.join(d, "main.tex"))
    secs = {t: (k, strip(r)) for t, k, r in split_sections(body(tex))}
    ok.append(("sections found", set(secs) == {"Abstract", "Author summary", "Results", "Discussion", "Extra"}))
    ok.append(("kinds", secs["Author summary"][0] == "unnumbered" and secs["Extra"][0] == "appendix"
               and secs["Results"][0] == "main"))
    ok.append(("stub dropped", "ignored" not in secs["Results"][1]))
    ok.append(("comment dropped", "999" not in secs["Results"][1]))
    ok.append(("e.g. not a break", len(sentences(secs["Results"][1])) == 3))
    avoid, cfg = story_config(d)
    ok.append(("terms parsed", len(avoid) == 4 and cfg["abstract_max"] == 2 and cfg["mean_max"] == 30))
    rows, findings, _, _ = analyse(tex, avoid, cfg)
    ok.append(("avoided term found", any("speedup" in x for x in findings)))
    ok.append(("typed number found", any("'3.27'" in x for x in findings)))
    ok.append(("appendix number not failed", not any("12.5" in x for x in findings)))
    ok.append(("abstract budget", any(x.startswith("abstract 3 words") for x in findings)))
    for name, good in ok:
        print("%s  %s" % ("ok  " if good else "FAIL", name))
    return 0 if all(g for _, g in ok) else 1


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--paper", help="paper directory (default: ../paper relative to this script, or ./paper)")
    ap.add_argument("--final", action="store_true")
    ap.add_argument("--dump-text", metavar="SECTION")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        sys.exit(selftest())
    paper = find_paper(a.paper)
    tex = flatten(os.path.join(paper, "main.tex"))
    if a.dump_text:
        sys.exit(dump_text(tex, a.dump_text))
    avoid, cfg = story_config(paper)
    findings = report(tex, avoid, cfg)
    label = "FAIL" if a.final else "note"
    for x in findings:
        print("%s: %s" % (label, x))
    print("check_prose: %s" % ("FAIL (%d)" % len(findings) if a.final and findings else "ok"))
    sys.exit(1 if a.final and findings else 0)


if __name__ == "__main__":
    main()

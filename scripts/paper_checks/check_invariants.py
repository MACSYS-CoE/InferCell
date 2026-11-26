"""Check that an editing pass changed the prose and nothing else (the /write-paper skill).

    python3 scripts/paper_checks/check_invariants.py snapshot            # record, before the pass
    python3 scripts/paper_checks/check_invariants.py compare             # after the pass; exit 1 on a change
    python3 scripts/paper_checks/check_invariants.py compare --section Results   # only that section's keys
    python3 scripts/paper_checks/check_invariants.py --paper docs/paper snapshot
    python3 scripts/paper_checks/check_invariants.py --selftest

An editing pass (/deslop, /prose-style, the joins pass) may reword sentences. It may not
add or drop a number, a citation, a cross-reference, a label or a display equation. The
snapshot records, for each section of the flattened manuscript (see check_prose.py):
  pn     the \\pn{} keys used            cite   the citation keys
  ref    the \\cref/\\ref/\\eqref targets   label  the \\label keys
  eq     each display equation, whitespace removed
as multisets, so a number quoted twice and then once is a change.

compare fails when the WHOLE PAPER's multiset of any kind differs from the snapshot, and
prints what was added and what was dropped. Keys that moved between sections (a sentence
moved by the joins pass) are printed as notes, not failures. Inline maths is not
checked: an edit pass may legitimately turn a phrase into a symbol.

The snapshot lives in paper/.write-paper/invariants.json. Standard library only.
"""
import argparse
import collections
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from check_prose import body, find_paper, flatten, split_sections  # noqa: E402

KINDS = ("pn", "cite", "ref", "label", "eq")
EQ_ENVS = r"equation\*?|align\*?|gather\*?|multline\*?|eqnarray\*?"


def keys(group):
    return [x.strip() for x in group.split(",") if x.strip()]


def extract(raw):
    out = {k: [] for k in KINDS}
    for g in re.findall(r"\\pn\{([^}]+)\}", raw):
        out["pn"].append(g.strip())
    for g in re.findall(r"\\(?:no)?cite[a-zA-Z]*\*?(?:\[[^\]]*\]){0,2}\{([^}]*)\}", raw):
        out["cite"] += keys(g)
    for g in re.findall(r"\\(?:[cC]ref|[cC]refrange|ref|eqref|autoref|pageref)\{([^}]+)\}", raw):
        out["ref"] += keys(g)
    out["label"] = [g.strip() for g in re.findall(r"\\label\{([^}]+)\}", raw)]
    eqs = re.findall(r"\\begin\{(%s)\}(.*?)\\end\{\1\}" % EQ_ENVS, raw, flags=re.S)
    eqs = [e for _, e in eqs] + re.findall(r"\\\[(.*?)\\\]", raw, flags=re.S)
    out["eq"] = ["".join(re.sub(r"\\label\{[^}]*\}", "", e).split()) for e in eqs]
    return out


def snapshot_of(tex):
    """{section title: {kind: sorted list}} for a flattened manuscript."""
    return {title: {k: sorted(v) for k, v in extract(raw).items()}
            for title, _, raw in split_sections(body(tex))}


def total(snap, section=None):
    out = {k: collections.Counter() for k in KINDS}
    for title, kinds in snap.items():
        if section and title.lower() != section.lower():
            continue
        for k in KINDS:
            out[k].update(kinds.get(k, []))
    return out


def short(item, kind):
    return (item[:70] + "...") if kind == "eq" and len(item) > 70 else item


def compare(before, after, section=None):
    """(failures, notes)."""
    fails, notes = [], []
    if section and not any(t.lower() == section.lower() for t in after):
        return ["section %r is missing after the pass" % section], []
    b, a = total(before, section), total(after, section)
    for k in KINDS:
        dropped, added = b[k] - a[k], a[k] - b[k]
        for item, n in sorted(dropped.items()):
            fails.append("%s dropped: %s%s" % (k, short(item, k), " (x%d)" % n if n > 1 else ""))
        for item, n in sorted(added.items()):
            fails.append("%s added: %s%s" % (k, short(item, k), " (x%d)" % n if n > 1 else ""))
    if not section:
        for title in sorted(set(before) | set(after)):
            for k in KINDS:
                x = collections.Counter(before.get(title, {}).get(k, []))
                y = collections.Counter(after.get(title, {}).get(k, []))
                if x != y and b[k] == a[k]:
                    notes.append("%s: %s keys moved in or out (%s)" % (title, k, ", ".join(
                        short(i, k) for i in sorted(((x - y) + (y - x)).keys()))[:160]))
        for title in sorted(set(before) - set(after)):
            notes.append("section %r no longer exists (renamed?)" % title)
    return fails, notes


def selftest():
    before = (r"\begin{document}\section{Results}We find \pn{a} and \pn{a} \citep{k1,k2}; see \cref{fig:x}."
              r"\begin{equation}y = mx\label{eq:y}\end{equation}\section{Discussion}As \citet{k3} said."
              r"\end{document}")
    s0 = snapshot_of(before)
    cases = {
        "reworded only": (before.replace("We find", "We measure"), 0),
        "dropped a number": (before.replace(r" and \pn{a}", ""), 1),
        "added a citation": (before.replace("said.", r"said \citep{k9}."), 1),
        "changed an equation": (before.replace("y = mx", "y = m x + c"), 1),
        "equation reflowed": (before.replace("y = mx", "y =\n  mx"), 0),
        "moved between sections": (before.replace(r"\section{Discussion}As \citet{k3} said.",
                                                  r"As \citet{k3} said.\section{Discussion}Fine."), 0),
    }
    ok = True
    for name, (tex, want) in cases.items():
        fails, notes = compare(s0, snapshot_of(tex))
        got = 1 if fails else 0
        ok &= got == want
        print("%s  %-24s %s" % ("ok  " if got == want else "FAIL", name, fails or notes or "no change"))
    fails, _ = compare(s0, snapshot_of(before.replace(r"As \citet{k3}", "As k3")), section="Results")
    ok &= not fails
    print("%s  %-24s %s" % ("ok  " if not fails else "FAIL", "--section scope", fails or "no change"))
    return 0 if ok else 1


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("action", nargs="?", choices=("snapshot", "compare"))
    ap.add_argument("--paper")
    ap.add_argument("--section", help="compare only this section's keys")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        sys.exit(selftest())
    if not a.action:
        ap.error("snapshot or compare")
    paper = find_paper(a.paper)
    store = os.path.join(paper, ".write-paper", "invariants.json")
    snap = snapshot_of(flatten(os.path.join(paper, "main.tex")))
    if a.action == "snapshot":
        os.makedirs(os.path.dirname(store), exist_ok=True)
        with open(store, "w") as f:
            json.dump(snap, f, indent=1, sort_keys=True)
        counts = total(snap)
        print("snapshot: %d sections; %s -> %s" % (len(snap), ", ".join(
            "%d %s" % (sum(counts[k].values()), k) for k in KINDS), os.path.relpath(store)))
        return
    if not os.path.isfile(store):
        sys.exit("no snapshot at %s: run snapshot before the pass" % os.path.relpath(store))
    with open(store) as f:
        before = json.load(f)
    fails, notes = compare(before, snap, a.section)
    for n in notes:
        print("note: " + n)
    for x in fails:
        print("FAIL: " + x)
    print("check_invariants: %s" % ("FAIL (%d)" % len(fails) if fails else "ok"))
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()

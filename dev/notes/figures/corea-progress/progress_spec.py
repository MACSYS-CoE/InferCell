"""Fig 1d/2d spec -- build status, parsed from spec/spec.md rather than typed.

Fig 1c draws the *model*. These two figures draw *how much of it exists*, and a
figure that drifts from the task list is worse than no figure, so nothing here
states a phase's status. `phases()` reads spec/spec.md §11: the headings, the
`**PR:**` line and the `- [x]` ticks. Only the mapping from a phase to the thing
it puts on the state graph is authored, in ELEMENT_PHASE and DELIVERS.

The reduction layer -- what Core A' greys, crosses and relabels -- is inherited
whole from ../corea-coupling/corea_spec.py. This file changes no model content.
"""
import os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "corea-coupling"))
from corea_spec import (GREY_FILL, GREY_LINE, GREY_TEXT, GREY_EDGE,      # noqa: F401
                        CHANGED_FILL, CROSS, REMOVED_NODES, CHANGED_NODES,
                        CHARGING_LABEL, charging_pos, DEAD_ENDPOINT,
                        DEAD_COUPLING, STRIP_LABELS, RELABEL, NEW_EDGES, CHG)

SPEC = os.path.join(HERE, "..", "..", "..", "..", "spec", "spec.md")

# --- palette: a second axis, deliberately not the reduction's ---------------
# The reduction owns red ✗ and grey. Progress owns green and hollow slate, so a
# reader never has to ask whether a mark means "cut" or "not written yet".
BUILT_FILL, BUILT_LINE, BUILT_TEXT = "#d5f0dd", "#1e8449", "#145a32"
TODO_FILL, TODO_LINE, TODO_TEXT = "#ffffff", "#94a3ab", "#5d6d76"


# --- the parser -------------------------------------------------------------
PHASE_H = re.compile(r"^### Phase (\d+) — (.*)$")
PR_LINE = re.compile(r"^\*\*PR:\*\* (.*)$")
PR_NUM = re.compile(r"^#(\d+) \(merged ([0-9-]+)")


class Phase:
    def __init__(self, n, title):
        self.n, self.title = n, title
        self.pr, self.merged, self.done, self.total = None, None, 0, 0

    @property
    def built(self):
        return self.pr is not None

    @property
    def tag(self):
        """`P3 · #44` when merged, `P10` when not."""
        return f"P{self.n} · #{self.pr}" if self.built else f"P{self.n}"

    def __repr__(self):
        return f"<Phase {self.n} {self.tag} {self.done}/{self.total}>"


def phases(path=SPEC):
    """Every phase in §11, with its PR and its tick count. Raises if §11 moved."""
    out, cur, in_tasks = {}, None, False
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith("## "):
            in_tasks = line.startswith("## 11.")
            if not in_tasks:
                cur = None
            continue
        if not in_tasks:
            continue
        m = PHASE_H.match(line)
        if m:
            cur = Phase(int(m.group(1)), _plain(m.group(2)))
            out[cur.n] = cur
            continue
        if cur is None:
            continue
        m = PR_LINE.match(line)
        if m:
            assert cur.pr is None, f"phase {cur.n} has two PR lines"
            n = PR_NUM.match(m.group(1))
            if n:
                cur.pr, cur.merged = int(n.group(1)), n.group(2)
            else:
                assert m.group(1) in ("_open_", "_not started_"), \
                    f"phase {cur.n}: unrecognised PR field {m.group(1)!r}"
            continue
        if line.startswith("- [x]") or line.startswith("- [ ]"):
            cur.total += 1
            cur.done += line.startswith("- [x]")

    missing = [n for n in range(18) if n not in out]
    assert not missing, f"spec/spec.md §11 has no heading for phase(s) {missing}"
    # Phase 17 states its own tasks are written once phase 16 has run, so a
    # zero-task phase is legitimate; a wholesale format change is not.
    assert sum(p.total for p in out.values()) > 100, \
        "§11 parsed with almost no `- [ ]` tasks -- the checkbox format moved"
    for p in out.values():
        assert not (p.built and p.done < p.total), \
            f"phase {p.n} is merged as #{p.pr} but only {p.done}/{p.total} ticked"
        assert not (p.built and p.total == 0), \
            f"phase {p.n} is merged as #{p.pr} with no task list"
    return out


def _plain(s):
    """Strip the markdown emphasis a heading may carry, e.g. phase 3's."""
    return re.sub(r"\*\*|\*|~~", "", s).strip().rstrip(".")


# --- what each phase puts on the state graph --------------------------------
# Authored. The mapping is fig 1c's own node ids, so a renamed node fails the
# build rather than silently losing a badge.
ELEMENT_PHASE = {
    "HOOK": 3, "REBUILD": 4, "GROWTH": 5,
    "Central": 6, "Transport": 7, "Nucleotide": 8, CHG: 9,
    "CME:Transcription": 10, "CME:Translation": 11, "CME:mRNAdecay": 12,
}

# Which phase turned each edge kind from a declaration into something that runs.
# Rows are fig 1c's key rows, in fig 1c's order. None = still declare-only.
EDGE_KIND_PHASE = [
    ("mass", "#1b4f72", "———",
     "shared state, continuous; gradients cross inside the ODE block", 1),
    ("currency", "#95a5a6", "- - -", "routed via the pool node", 1),
    ("deferred counter", "#b9770e", "- - -",
     "CME accrues a cost, the hook debits it one step later", 3),
    ("catalytic", "#7d3c98", "· · · ·",
     "counts enter rate laws as parameters, no mass flows", 3),
    ("rate constant", "#1e8449", "- - -",
     "pools re-enter the CME at the 60 s rebuild, piecewise-constant", 4),
    ("volume", "#935116", "———",
     "counts set surface area, which sets V, which rescales concentrations", 5),
    ("clamped", "#c0392b", "- - ⊣", "a real dependence replaced by a constant", None),
]

# One line per phase, naming what it delivers *on the state graph*. Phases that
# put nothing there say so; that absence is the point of fig 2d.
DELIVERS = {
    0: "— nothing in the figure —",
    1: "the mass and currency arrows inside the ODE grid",
    2: "the CME block being three coexisting boxes at all",
    3: "the HOOK device · catalytic and deferred-counter arrows",
    4: "the REBUILD device · the rate-constant arrows",
    5: "the GROWTH device · the volume arrow",
    6: "box: Central (10 rxns)",
    7: "box: Transport (6 rxns)",
    8: "box: Nucleotide (5 rxns)",
    9: "box: tRNA charging (the box that moved)",
    10: "box: transcription (17 genes)",
    11: "box: translation (17, plus ptsG translocation)",
    12: "box: mRNA decay (17)",
    13: "the whole figure asserted complete, as an object in code",
    14: "conservation checks over the figure",
    15: "— nothing in the figure —",
    16: "— nothing in the figure —",
    17: "— nothing in the figure —",
}

# The dependency edges §11's parallelism paragraph states, and only those.
# (from, to, style, note)
DEPS = [
    (0, 1, "solid", ""), (1, 2, "solid", ""), (2, 3, "solid", ""),
    (3, 4, "solid", ""), (4, 5, "solid", ""),
    (1, 6, "solid", ""), (1, 7, "solid", ""), (1, 8, "solid", ""),
    (8, 9, "solid", "charging contributes to the pools recycling owns"),
    (2, 10, "solid", ""), (10, 11, "solid", ""), (10, 12, "solid", ""),
    (9, 11, "dashed", "can proceed on a double; R1 records the cost"),
    (5, 13, "solid", ""), (6, 13, "solid", ""), (7, 13, "solid", ""),
    (9, 13, "solid", ""), (11, 13, "solid", ""), (12, 13, "solid", ""),
    (13, 14, "solid", ""), (14, 15, "solid", ""), (15, 16, "solid", ""),
    (16, 17, "solid", ""),
]

# fig 2d's grid: (column, row). Row 0 is the framework spine, rows 1-5 the
# module fan-out it enables, row 6 the strictly serial tail. Two short rows
# rather than one long one, because eighteen chips in a line are unreadable.
GRID = {
    0: (0, 0), 1: (1, 0), 2: (2, 0), 3: (3, 0), 4: (4, 0), 5: (5, 0),
    6: (2, 1), 7: (2, 2), 8: (2, 3), 9: (4, 3),
    10: (3, 4), 11: (4, 4), 12: (4, 5),
    13: (1, 6), 14: (2, 6), 15: (3, 6), 16: (4, 6), 17: (5, 6),
}

# Which rows open a band, for the labels down the left-hand side.
ROW_BAND = {
    0: ("the framework", "no syn3A biology at all — D0's order, and every phase "
        "merged so far"),
    1: ("the fan-out — ODE block",
        "6, 7 and 8 are mutually independent; 9 waits on 8"),
    4: ("the fan-out — CME block",
        "10 before 11 and 12, which both read the transcripts 10 owns"),
    6: ("assembly, validation, inference",
        "strictly serial — each check catches what the next would mask"),
}

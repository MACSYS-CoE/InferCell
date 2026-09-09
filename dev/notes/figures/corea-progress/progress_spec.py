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
# Only what the derivation needs. GREY_EDGE, RELABEL, STRIP_LABELS and the two
# dead-edge lists are read by fig 1c's and fig 1r's own edge stages, which
# make_progress re-runs; the rest of fig 1c's node vocabulary is not, because
# fig 1d authors its own boxes.
from corea_spec import (GREY_EDGE, CHANGED_FILL, REMOVED_NODES, DEAD_ENDPOINT,
                        DEAD_COUPLING, STRIP_LABELS, RELABEL, NEW_EDGES, CHG)

SPEC = os.path.join(HERE, "..", "..", "..", "..", "spec", "spec.md")

# --- palette: a second axis, deliberately not the reduction's ---------------
# The reduction owns red ✗ and grey. Progress owns green and hollow slate, so a
# reader never has to ask whether a mark means "cut" or "not written yet".
BUILT_FILL, BUILT_LINE, BUILT_TEXT = "#d5f0dd", "#1e8449", "#145a32"
TODO_FILL, TODO_LINE, TODO_TEXT = "#ffffff", "#94a3ab", "#5d6d76"


# Every phase §11 must carry, written out. A literal rather than a `range` so a
# lettered phase is as load-bearing as a numbered one, and so adding one is a
# visible edit here instead of an off-by-one in an arithmetic expression.
EXPECTED_PHASES = ["0", "1", "2", "3", "4", "5", "5b",
                   "6", "7", "8", "9", "10", "11", "12",
                   "13", "14", "15", "16", "17"]

# --- the parser -------------------------------------------------------------
PHASE_H = re.compile(r"^### Phase (\d+)([a-z]?) — (.*)$")
PR_LINE = re.compile(r"^\*\*PR:\*\* (.*)$")
PR_NUM = re.compile(r"^#(\d+) \(merged ([0-9-]+)")


class Phase:
    def __init__(self, n, key, title):
        # `n` is the label as the spec writes it -- "5" or "5b" -- because it is
        # what every badge prints and what ELEMENT_PHASE, DELIVERS, GRID and
        # DEPS key on. Ordering goes through `key`, never through `n`: sorting
        # phase labels as strings would put "10" before "5". `key` is stored
        # rather than derived, since PHASE_H has already split the two parts.
        self.n, self.key, self.title = n, key, title
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
        # A heading the regex cannot read would be skipped in silence and its
        # phase would simply be absent from both figures -- or, worse, its PR
        # line would be attributed to the phase above it.
        assert m or not line.startswith("### Phase "), \
            f"§11 heading not parsed, so its phase would vanish: {line!r}"
        if m:
            cur = Phase(m.group(1) + m.group(2),
                        (int(m.group(1)), m.group(2)),
                        _plain(m.group(3)))
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

    missing = [n for n in EXPECTED_PHASES if n not in out]
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
    "HOOK": "3", "REBUILD": "4", "GROWTH": "5",
    "Central": "6", "Transport": "7", "Nucleotide": "8", CHG: "9",
    "CME:Transcription": "10", "CME:Translation": "11", "CME:mRNAdecay": "12",
}

# Which phase turned each edge kind from a declaration into something that runs.
# Rows are fig 1c's key rows, in fig 1c's order. None = still declare-only.
EDGE_KIND_PHASE = [
    ("mass", "#1b4f72", "———",
     "shared state, continuous; gradients cross inside the ODE block", "1"),
    ("currency", "#95a5a6", "- - -", "routed via the pool node", "1"),
    ("deferred counter", "#b9770e", "- - -",
     "CME accrues a cost, the hook debits it one step later", "3"),
    ("catalytic", "#7d3c98", "· · · ·",
     "counts enter rate laws as parameters, no mass flows", "3"),
    ("rate constant", "#1e8449", "- - -",
     "pools re-enter the CME at the 60 s rebuild, piecewise-constant", "4"),
    ("volume", "#935116", "———",
     "counts set surface area, which sets V, which rescales concentrations;"
     "<br/>the inbound half — geometry into an ODE rate law — arrives with 5b",
     "5"),
    ("clamped", "#c0392b", "- - ⊣", "a real dependence replaced by a constant", None),
]

# Drawn in the figure and explained in the key, but deliberately *not* one of
# the seven kinds and never counted among them: an intra-block coupling is not a
# CouplingEdge at all. Spec §12, 2026-09-04 settled that a jump module's peer
# writes are gated on `written_states` rather than on an edge.
EXTRA_KEY_ROWS = [
    ("intra-block", "#7d3c98", "———",
     "<i>not a CouplingEdge</i> — a jump module reads or writes a declared "
     "peer’s state,<br/>gated on <b>written_states</b>, which spec §12 "
     "(2026-09-04) settled", "2"),
]

# One line per phase, naming what it delivers *on the state graph*. Phases that
# put nothing there say so; that absence is the point of fig 2d.
DELIVERS = {
    "0": "— nothing in the figure —",
    "1": "the mass and currency arrows inside the ODE grid",
    "2": "the CME block being three coexisting boxes at all",
    "3": "the HOOK device · catalytic and deferred-counter arrows",
    "4": "the REBUILD device · the rate-constant arrows",
    "5": "the GROWTH device · the volume arrow",
    "5b": "the volume arrow's inbound half — geometry into a rate law",
    "6": "box: Central (10 rxns)",
    "7": "box: Transport (6 rxns)",
    "8": "box: Nucleotide (5 rxns)",
    "9": "box: tRNA charging (the box that moved)",
    "10": "box: transcription (17 genes)",
    "11": "box: translation (17, plus ptsG translocation)",
    "12": "box: mRNA decay (17)",
    "13": "the whole figure asserted complete, as an object in code",
    "14": "conservation checks over the figure",
    "15": "— nothing in the figure —",
    "16": "— nothing in the figure —",
    "17": "— nothing in the figure —",
}

# The dependency edges §11's parallelism paragraph states, and only those.
# (from, to, style, note)
DEPS = [
    ("0", "1", "solid", ""), ("1", "2", "solid", ""), ("2", "3", "solid", ""),
    ("3", "4", "solid", ""), ("4", "5", "solid", ""), ("5", "5b", "solid", ""),
    ("1", "6", "solid", ""), ("1", "7", "solid", ""), ("1", "8", "solid", ""),
    ("5b", "7", "solid", "the inbound volume channel task 7.2 needs"),
    ("8", "9", "solid", "charging contributes to the pools recycling owns"),
    ("2", "10", "solid", ""), ("10", "11", "solid", ""), ("10", "12", "solid", ""),
    ("9", "11", "dashed", "can proceed on a double; R1 records the cost"),
    ("5", "13", "solid", ""), ("6", "13", "solid", ""), ("7", "13", "solid", ""),
    ("9", "13", "solid", ""), ("11", "13", "solid", ""), ("12", "13", "solid", ""),
    ("13", "14", "solid", ""), ("14", "15", "solid", ""), ("15", "16", "solid", ""),
    ("16", "17", "solid", ""),
]

# ---------------------------------------------------------------------------
# fig 1d's layout. Fig 1c is fig 1's canvas with the reduction drawn *onto* it,
# so it carries every module Core A′ deletes and every coupling it cuts, greyed.
# That is the right figure for "what did the reduction do" and the wrong one for
# "what is built": half its ink is about a model we are not building. Fig 1d
# therefore keeps only what survives, and lays it out fresh.
#
# The shape is the one sentence the figure is trying to say: two blocks, three
# coupling devices, one cycle. CME on the left, ODE on the right, the 1 s hook
# between them; the 60 s rebuild and the volume chain below, closing the loop
# the other way.
# ---------------------------------------------------------------------------

ODE_C, CME_C, HOOK_C, RB_C, GR_C, CL_C = ("#1b4f72", "#7d3c98", "#117864",
                                          "#1e8449", "#935116", "#c0392b")
CUR_C = "#b3bfc4"
ODE_PANEL, CME_PANEL = "#f4f8fb", "#faf6fc"

POS = {
    # the CME block, a triangle so transcription reaches both its consumers
    "CME:Transcription": (190, 900), "CME:Translation": (190, 660),
    "CME:mRNAdecay": (425, 780),
    # the ODE block, a 2x2 grid over the shared currency pool. It sits this far
    # right because five labelled couplings leave the hook for it, and they need
    # the gap: closer together, their labels stack on top of one another.
    # Charging sits nearest the hook and Nucleotide furthest, so the hook's
    # single edge to charging is short and its label falls in open ground. The
    # other way round that edge grazes the top of whichever box it passes, and
    # every position for its label lands on one. Central over Nucleotide is a
    # free consequence: their two mass edges become a clean vertical pair.
    "Transport": (1385, 900), "Central": (1680, 900),
    "CME:tRNAcharging": (1385, 690), "Nucleotide": (1680, 690),
    "CURRENCY": (1532, 555),
    # the three devices, and the clamped medium above the block it feeds
    "HOOK": (740, 790), "REBUILD": (560, 420), "GROWTH": (1130, 400),
    "Medium": (1532, 1075),
    # panels, keys, title
    "CMEPANEL": (307, 795), "ODEPANEL": (1532, 760),
    "KEY": (430, 165), "PROGKEY": (1420, 172), "TITLE": (900, 1150),
}

# The two blocks drawn as ground rather than as a cluster, so an edge can
# terminate on a *block* where the coupling is into the block rather than into
# one module of it. Padding is (sides, bottom, top) in points, measured from the
# module *centres* the panel has to contain -- the top pad carries the block
# title, which is why it is the largest.
PANELS = {
    "CMEPANEL": (95, 60, 100, CME_C, CME_PANEL,
                 "CME block", "Gillespie / Lattice Microbes · stochastic · not "
                 "differentiable<br/>three modules · 17 × 3 reactions, plus one "
                 "translocation"),
    "ODEPANEL": (112, 58, 100, ODE_C, ODE_PANEL,
                 "ODE block", "odecell → LSODA · deterministic · differentiable"
                 "<br/>four modules · <b>190 → 22 reactions</b> · ~32 dynamic "
                 "states and 5 chemostats"),
}


def _mod(name, old_new, note=""):
    lab = f'<b>{name}</b><br/><font point-size="8.5"><b>{old_new}</b></font>'
    return lab + (f'<br/><font point-size="7.5">{note}</font>' if note else "")


# Module labels are fig 1c's, trimmed: with the greyed half of the figure gone
# there is no longer anything to disambiguate against.
NODE_LABEL = {
    "Transport": _mod("Transport", "87 → 6 rxns", "PTS cascade ×5 · lactate export"),
    "Central": _mod("Central", "26 → 10 rxns", "glycolysis through lactate · NOX dropped"),
    "Nucleotide": _mod("Nucleotide", "46 → 5 rxns", "PGK3, PYK3 · ADK1, PPA, GK1"),
    "CME:tRNAcharging": _mod("tRNA charging", "20 × 5 rxns → 1 lumped step",
                             "moved CME → ODE · owns both tRNA pools"),
    "CME:Transcription": _mod("transcription", "455 → 17 genes", ""),
    "CME:Translation": _mod("translation", "455 → 17 genes", "plus ptsG translocation"),
    "CME:mRNAdecay": _mod("mRNA decay", "455 → 17 mRNAs", ""),
}

DEVICE_LABEL = {
    "HOOK": ("#d4efdf", HOOK_C,
             '<b>the hook</b> — state<br/><font point-size="9.5">every <b>1.0 s</b> '
             'of cell time</font><br/> <br/>'
             '<font point-size="8.5">counts → mM · integrate 1 s (LSODA) · write '
             'counts back</font><br/>'
             '<font point-size="8"><i>~6,300 calls per cell cycle · gradients do '
             'not cross here</i></font>'),
    "REBUILD": (CHANGED_FILL, RB_C,
                '<b>the CME rebuild</b> — rate constants<br/>'
                '<font point-size="9.5">every <b>60 s</b>, MinCell_restart.py</font>'
                '<br/> <br/><font point-size="8.5">live pools → k_transcription, '
                'k_translation<br/>the whole CME is reconstructed</font><br/>'
                '<font point-size="8"><i>the only ODE→CME direction</i></font>'),
    "GROWTH": (CHANGED_FILL, GR_C,
               '<b>growth</b> — surface area<br/>'
               '<font point-size="9.5">every 1.0 s, in neither block</font><br/> <br/>'
               '<font point-size="8.5">CellSA = 28.0 nm² × membrane protein count'
               '<br/>r = √(CellSA/4π) &#160;&#160; V = (4/3)πr³</font><br/>'
               '<font point-size="8"><i>ptsG is the only one: 831 copies.</i> '
               '<b>V reaches ~1.07×</b>;<br/><i>report fractional growth, never a '
               'doubling time</i></font>'),
}

MEDIUM_LABEL = ('<b>growth medium</b> — <b>17 → 5 clamped concentrations</b><br/>'
                '<font point-size="8.5">glucose_e (40 mM) · CTP · UTP · aa pool · O₂'
                '</font>')
CURRENCY_LABEL = ('<b>currency pool</b> &#160;<font point-size="8.5">ATP · ADP · AMP · '
                  'P<sub>i</sub> · PP<sub>i</sub> · NAD(P)(H)</font><br/>'
                  '<font point-size="7.5"><i>not a module — a drawing device.</i> '
                  'Shared by all four ODE modules; NAD(P)(H) survives through LDH_L '
                  'alone.</font>')


# Every coupling Core A′ keeps, one row each, in fig 1c's colours and styles.
# Authored so the labels can be trimmed and the layout can breathe -- and then
# *checked*: make_progress asserts this table is exactly fig 1c's live edge set,
# so a coupling cannot be dropped by being forgotten. Additions are declared in
# ADDED_EDGES below and must be justified there.
#
# (src, dst, colour, style, penwidth, label, extra attributes)
E = lambda a, b, c, st, pw, lab="", x="": (a, b, c, st, pw, lab, x)


def LBL(colour, text, size=9):
    """An edge label on an opaque backing.

    Graphviz centres a label on its own edge, so the line runs through the
    text. Every label in this figure therefore sits on a white table cell with
    no border: invisible against the page, and it hides the couple of
    millimetres of line behind it."""
    return ('<<table border="0" cellborder="0" cellspacing="0" cellpadding="1" '
            f'bgcolor="white"><tr><td><font color="{colour}" point-size="{size}">'
            f'{text}</font></td></tr></table>>')


TEE = "arrowhead=tee"
CTR_C = "#b9770e"                      # deferred counter
DEBIT_NTP = "<b>debit NTP pools</b><br/>clamped at zero, deficit carried"
CREDIT_NMP = "<b>credit NMP</b><br/>decay recycling, unclamped"
DEBIT_TRNA = ("<b>debit charged tRNA · credit uncharged</b><br/>"
              "~553 residues/s against a pool of order 10³")
GLUCOSE_CLAMP = "glucose_e<br/><i>no feedback on the medium</i>"
NOHEAD = "arrowhead=none"

EDGES = [
    # --- mass, inside the ODE block -----------------------------------------
    E("Central", "Nucleotide", ODE_C, "solid", 2.1, "13dpg, pep"),
    E("Nucleotide", "Central", ODE_C, "solid", 1.9, "3pg, pyr"),
    E("Central", "Transport", ODE_C, "solid", 1.7, "lac, pep"),
    E("Transport", "Central", ODE_C, "solid", 1.3, "g6p, pyr"),
    # --- currency, routed via the pool device -------------------------------
    E("Transport", "CURRENCY", CUR_C, "dashed", 0.9, "", NOHEAD),
    E("Central", "CURRENCY", CUR_C, "dashed", 0.9, "", NOHEAD),
    E("Nucleotide", "CURRENCY", CUR_C, "dashed", 0.9, "", NOHEAD),
    E("CME:tRNAcharging", "CURRENCY", CUR_C, "dashed", 0.9, "", NOHEAD),
    # --- mass, inside the CME block -----------------------------------------
    E("CME:Transcription", "CME:Translation", CME_C, "solid", 1.4, "17 mRNAs"),
    E("CME:Transcription", "CME:mRNAdecay", CME_C, "solid", 1.4, "17 mRNAs"),
    # --- deferred counters: the CME accrues, the hook debits one step later --
    E("CME:Transcription", "HOOK", CTR_C, "dashed", 1.5, "atp, ctp, gtp, utp"),
    E("CME:Translation", "HOOK", CTR_C, "dashed", 1.5, "atp, gtp · charged tRNA"),
    E("CME:mRNAdecay", "HOOK", CTR_C, "dashed", 1.5, "amp, atp, cmp, gmp +1"),
    E("HOOK", "Central", CTR_C, "dashed", 2.0,
      "<b>debit ATP</b><br/>trsc, translat, mRNAdeg, transloc"),
    # Two counters run hook -> Nucleotide, so their labels are anchored at the
    # head at different angles; a mid-edge label would put both in one place.
    E("HOOK", "Nucleotide", CTR_C, "dashed", 2.4, "",
      f'headlabel={LBL(CTR_C, DEBIT_NTP)}, labeldistance=6.2, labelangle=-15'),
    E("HOOK", "Nucleotide", CTR_C, "dotted", 1.6, "",
      f'headlabel={LBL(CTR_C, CREDIT_NMP)}, labeldistance=5.5, labelangle=16'),
    E("HOOK", "CME:tRNAcharging", CTR_C, "dashed", 2.0, DEBIT_TRNA),
    # --- catalytic: counts enter rate laws, no mass flows --------------------
    E("CME:Translation", "HOOK", CME_C, "dotted", 2.0, "protein counts"),
    E("HOOK", "ODEPANEL", CME_C, "dotted", 2.0,
      "→ <b>enzyme concentrations</b><br/>17 single-gene GPRs, one-to-one<br/>"
      "<i>into the block, drawn once</i>"),
    # --- rate constant: the 60 s return path --------------------------------
    E("Nucleotide", "REBUILD", RB_C, "dashed", 1.7,
      "atp, gtp &#160;<i>ctp, utp chemostatted</i>"),
    E("CME:tRNAcharging", "REBUILD", RB_C, "dashed", 1.7, "charged-tRNA pool"),
    E("REBUILD", "CME:Transcription", RB_C, "dashed", 2.4,
      "<b>k_transcription</b><br/>per gene, from live NTP conc."),
    E("REBUILD", "CME:Translation", RB_C, "dashed", 2.4,
      "<b>k_translation</b><br/>per gene, from the live charged-tRNA pool"),
    # --- volume: on neither clock -------------------------------------------
    E("CME:Translation", "GROWTH", GR_C, "solid", 1.6, "1 membrane protein locus (ptsG)"),
    E("GROWTH", "HOOK", GR_C, "solid", 2.4,
      "<b>V rescales every count↔mM conversion</b><br/>and so dilutes all four ODE "
      "modules"),
    E("GROWTH", "REBUILD", GR_C, "solid", 1.6, "CellV → getConc"),
    # --- clamped: a real dependence replaced by a constant -------------------
    E("Medium", "Transport", CL_C, "dashed", 1.7, "",
      TEE + f', headlabel={LBL(CL_C, GLUCOSE_CLAMP)}, labeldistance=6.5, labelangle=20'),
    E("Medium", "CME:tRNAcharging", CL_C, "dashed", 1.7, "",
      TEE + f', headlabel={LBL(CL_C, "aa pool,<br/>chemostatted")}, labeldistance=5.5, labelangle=-19'),
    E("CME:Transcription", "CME:Transcription", CL_C, "dashed", 1.3,
      "RnaPconc = 187, frozen", TEE + ", tailport=w, headport=w"),
    E("CME:Translation", "CME:Translation", CL_C, "dashed", 1.3,
      "ribosomeConc = 503, frozen", TEE + ", tailport=w, headport=w"),
    E("GROWTH", "GROWTH", CL_C, "dashed", 1.3,
      "V clamped at 2× initial.<br/><b>There is no division event.</b>",
      TEE + ", tailport=e, headport=e"),
]

# Edges fig 1c does not draw between these two nodes, with the reason. Fig 1c
# draws the catalytic channel's ODE half into Cofactor, because fig 1 does and
# fig 1c holds fig 1's geometry; Core A′ deletes Cofactor, so fig 1c greys that
# edge and puts the real targets in a footnote. With the deleted modules gone
# there is nothing to grey and no footnote to hang it on, so the edge is drawn
# where the channel actually runs -- into the block, once, as fig 1 itself
# labels it.
ADDED_EDGES = {
    ("HOOK", "ODEPANEL"): "the catalytic channel's ODE half, which fig 1c can only "
                          "draw into a module Core A′ deletes (its CATNOTE)",
}


# fig 2d's grid: (column, row). Row 0 is the framework spine, rows 1-5 the
# module fan-out it enables, row 6 the strictly serial tail. Two short rows
# rather than one long one, because nineteen chips in a line are unreadable.
GRID = {
    "0": (0, 0), "1": (1, 0), "2": (2, 0), "3": (3, 0), "4": (4, 0),
    "5": (5, 0), "5b": (6, 0),
    "6": (2, 1), "7": (2, 2), "8": (2, 3), "9": (4, 3),
    "10": (3, 4), "11": (4, 4), "12": (4, 5),
    "13": (1, 6), "14": (2, 6), "15": (3, 6), "16": (4, 6), "17": (5, 6),
}

# Which rows open a band, for the labels down the left-hand side.
ROW_BAND = {
    0: ("the framework", "no syn3A biology at all — D0's order, and every phase "
        "merged so far except 5b, which is last in the band and is what "
        "unblocks the fan-out"),
    1: ("the fan-out — ODE block",
        "6, 7 and 8 are mutually independent; 9 waits on 8"),
    4: ("the fan-out — CME block",
        "10 before 11 and 12, which both read the transcripts 10 owns"),
    6: ("assembly, validation, inference",
        "strictly serial — each check catches what the next would mask"),
}

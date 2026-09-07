"""Fig 1d and fig 2d -- how much of Core A′ exists, as of the current spec.

    python3 make_progress.py     # -> fig1d_state_graph_progress.{pdf,png}
                                 #    fig2d_phase_roadmap.{pdf,png}

Fig 1d is fig 1c with a build-status overlay: it re-runs ../corea-coupling's
builder and adds badges, so no model content is restated here and the two
figures cannot disagree about the model. Fig 2d is drawn from scratch, because
the thing it shows -- the phase list in dependency order, including the five
phases that put nothing on a state graph -- has no representation in fig 1c.

Neither figure states a status. Both read spec/spec.md §11 through
progress_spec.phases().
"""
import datetime, os, re, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
COREA = os.path.join(HERE, "..", "corea-coupling")
BUILD = os.path.join(HERE, "build")
OUT = os.path.join(BUILD, "out")

sys.path.insert(0, HERE)
sys.path.insert(0, COREA)
import progress_spec as S
import make_corea as C

# Point fig 1c's builder, and through it fig 1r's stage-2/4 helpers, at this
# directory's build tree and this directory's spec. Same idiom make_corea.py
# uses on make_reduced (make_corea.py:31-35), one level further down the chain.
C.BUILD, C.OUT, C.S = BUILD, OUT, S
C.M.BUILD, C.M.OUT, C.M.S = BUILD, OUT, S

STEM1 = "fig1d_state_graph_progress"
STEM2 = "fig2d_phase_roadmap"


# --- helpers ----------------------------------------------------------------
def node_id(stmt):
    m = re.match(r'\s*"([^"]+)" \[', stmt)
    return m.group(1) if m and "->" not in stmt.split("[")[0] else None


def drop(stmts, ids):
    """Remove whole node statements by id, and assert each was there to remove."""
    kept = [s for s in stmts if node_id(s) not in ids]
    gone = {node_id(s) for s in stmts} & ids
    assert gone == ids, f"nothing to drop for {sorted(ids - gone)}"
    return kept


def measure(path, ref=None):
    """Every node's rendered centre and size at its pinned position.

    The badge has to sit on the corner of the box as *drawn*, and three boxes
    are not the size fig 1 drew them: the charging box is re-labelled and
    re-pinned, and the two keys are authored here. Asking neato is the only way
    to know, and it is what stage 2 already does for fig 1's own geometry.
    """
    dot = os.path.join(OUT, "measure.dot")
    subprocess.run(["neato", "-n", "-Tdot", path, "-o", dot], check=True)
    # -Tdot renormalises the origin, so a node overflowing the bounding box
    # shifts every coordinate it reports. Realigned against `ref` below.
    t = open(dot).read()
    g = {}
    for stmt in t.split("];"):
        head = stmt.split("[")[0]
        m = re.search(r'"?([A-Za-z_:][A-Za-z_:0-9]*)"?\s*$', head.strip())
        if not m or "->" in head:
            continue
        p = re.search(r'\bpos="([-\d.,e]+)"', stmt)
        w = re.search(r'\bwidth=([\d.]+)', stmt)
        h = re.search(r'\bheight=([\d.]+)', stmt)
        if p and w and h:
            x, y = (float(v) for v in p.group(1).split(","))
            g[m.group(1)] = (x, y, float(w.group(1)) * 72, float(h.group(1)) * 72)
    if ref:
        name, (rx, ry) = ref
        dx, dy = g[name][0] - rx, g[name][1] - ry
        g = {k: (x - dx, y - dy, w, h) for k, (x, y, w, h) in g.items()}
    return g


def render(dot_lines, stem, tag):
    path = os.path.join(BUILD, stem + ".dot")
    open(path, "w").write("\n".join(dot_lines))
    for fmt in ("pdf", "png"):
        subprocess.run(["neato", "-n", f"-T{fmt}", "-Gdpi=150", path,
                        "-o", os.path.join(HERE, f"{stem}.{fmt}")], check=True)
    size = subprocess.run(["pdfinfo", os.path.join(HERE, stem + ".pdf")],
                          capture_output=True, text=True
                          ).stdout.split("Page size:")[1].split("\n")[0].strip()
    print(f"ok  {tag}  canvas {size}  ({len(dot_lines)} statements)")
    return size


# --- fig 1d: the state graph, annotated -------------------------------------
def badge(nid, x, y, ph):
    """A pill on a box's top-right corner. Green and filled once the phase has
    merged, hollow while it has not. Never drawn on a removed module: those are
    on nobody's task list, and their ✗ already means something else."""
    fill, line, text = ((S.BUILT_FILL, S.BUILT_LINE, S.BUILT_TEXT) if ph.built
                        else (S.TODO_FILL, S.TODO_LINE, S.TODO_TEXT))
    mark = "&#10004;" if ph.built else "&#9744;"
    return (f'  "{nid}" [pos="{x:.1f},{y:.1f}!", shape=plaintext, style="", label=<'
            f'<table border="1" color="{line}" cellpadding="3" cellspacing="0" '
            f'bgcolor="{fill}"><tr><td><font color="{text}" point-size="10">'
            f'<b>{mark} {ph.tag}</b></font></td></tr></table>>];')


def channel_badge(nid, x, y, ph, what):
    """The two phases that built a *property* of a block rather than a box:
    phase 1's contribution channel and phase 2's jump composition. Neither has
    a box of its own to sit on, and badging one of the boxes would claim the
    box exists."""
    return (f'  "{nid}" [pos="{x:.1f},{y:.1f}!", shape=plaintext, style="", label=<'
            f'<table border="1" color="{S.BUILT_LINE}" cellpadding="4" cellspacing="0" '
            f'bgcolor="{S.BUILT_FILL}"><tr><td align="left">'
            f'<font color="{S.BUILT_TEXT}" point-size="10"><b>&#10004; {ph.tag}</b></font>'
            f'<font color="{S.BUILT_TEXT}" point-size="8.5"> &#160;{what}</font>'
            '</td></tr></table>>];')


def key_node(x, y, ps):
    """Fig 1c's edge-kind key, gaining the column the whole figure turns on.

    The seven rows and their descriptions are fig 1c's, word for word
    (make_corea.py:110-128); only the third column is new."""
    rows = []
    for name, colour, glyph, what, ph in S.EDGE_KIND_PHASE:
        if ph is None:
            since = f'<font color="{S.TODO_TEXT}"><b>&#9744;</b> declare-only</font>'
        elif ps[ph].built:
            since = f'<font color="{S.BUILT_TEXT}"><b>&#10004; {ps[ph].tag}</b></font>'
        else:
            since = f'<font color="{S.TODO_TEXT}"><b>&#9744; P{ph}</b></font>'
        rows.append(f'<tr><td align="right"><font color="{colour}"><b>{glyph}</b></font></td>'
                    f'<td align="left"><b>{name}</b> · {what}</td>'
                    f'<td align="center"><font point-size="9">{since}</font></td></tr>')
    live = sum(1 for *_, ph in S.EDGE_KIND_PHASE if ph is not None and ps[ph].built)
    return ('  "KEY" [pos="%.1f,%.1f!", shape=plaintext, style="", label=<'
            '<table border="1" color="#bdc3c7" cellpadding="5" cellspacing="0" bgcolor="white">'
            '<tr><td colspan="2" align="center">'
            '<b>edge kinds</b> — what actually crosses each interface</td>'
            '<td align="center"><b>executes since</b></td></tr>'
            + "".join(rows) +
            '<tr><td colspan="3" align="left"><font point-size="8.5">'
            f'<b>{live} of {len(S.EDGE_KIND_PHASE)} kinds run.</b> The seventh validates, '
            'resolves and reports, but nothing executes it.</font></td></tr>'
            '</table>>];') % (x, y)


def progress_key(x, y, ps):
    """The second axis, stated beside the first so no reader conflates them."""
    built = sum(p.built for p in ps.values())
    done = sum(p.done for p in ps.values())
    total = sum(p.total for p in ps.values())
    live = sum(1 for *_, ph in S.EDGE_KIND_PHASE if ph is not None and ps[ph].built)
    example = max((p for p in ps.values() if p.built), key=lambda p: p.n)
    todo = min((p for p in ps.values() if not p.built), key=lambda p: p.n)
    row = lambda mark, colour, fill, what: (
        f'<tr><td align="center" bgcolor="{fill}"><font color="{colour}" point-size="10">'
        f'<b>{mark}</b></font></td><td align="left"><font point-size="9">{what}</font>'
        '</td></tr>')
    return (f'  "PROGKEY" [pos="{x:.1f},{y:.1f}!", shape=plaintext, style="", label=<'
            f'<table border="1" color="{S.BUILT_LINE}" cellpadding="5" cellspacing="0" '
            'bgcolor="white">'
            '<tr><td colspan="2" align="center"><b>what is built</b> — a second axis, '
            'and not the reduction’s</td></tr>'
            + row(f"&#10004; {example.tag}", S.BUILT_TEXT, S.BUILT_FILL,
                  "<b>merged</b> · every task in the phase is ticked, and the PR is "
                  "the record of what it cost")
            + row(f"&#9744; P{todo.n}", S.TODO_TEXT, S.TODO_FILL,
                  "<b>not written</b> · the number is the spec §11 phase that writes it")
            + row("&#10007;", S.CROSS, S.GREY_FILL,
                  "<b>unchanged from fig 1c</b> · cut from the <i>model</i> by the "
                  "reduction. Says nothing about code.")
            + '<tr><td colspan="2" align="left"><font point-size="9">'
              f'<b>{built} of {len(ps)} phases have merged, and every one of them is '
              f'framework.</b> {live} of {len(S.EDGE_KIND_PHASE)} edge kinds execute;'
              '<br align="left"/><b>not one box in this figure exists in code.</b> '
              f'<font point-size="8.5">{done} of {total} tasks ticked; every badge here'
              '<br align="left"/>is parsed from spec §11, never typed.<br align="left"/>'
              '</font>'
              '</font></td></tr></table>>];')


def title_node(x, y, ps):
    built = [p for p in ps.values() if p.built]
    asof = max(p.merged for p in built)
    prs = f"#{min(p.pr for p in built)}–#{max(p.pr for p in built)}"
    return (f'  "TITLE" [pos="{x:.1f},{y:.1f}!", shape=plaintext, style="", label=<'
            '<font point-size="21"><b>State coupling in reduced syn3A — Core A′: '
            'what is built</b></font><br/>'
            '<font point-size="11">fig 1c, with every module and device badged by the '
            f'spec §11 phase that writes it · <b>{len(built)} of {len(ps)} phases merged</b> '
            f'({prs}) · as of {asof}</font>>];')


def fig1d(ps):
    src = C.baseline()
    g, bb = C.M.geometry()
    header = [l.replace("overlap=false", "overlap=true") for l in src[:5]]

    # Fig 1c's own statements, less the two nodes this figure re-authors and the
    # note whose slot it needs. MOVENOTE explains why charging sits where it
    # does -- settled, recorded in spec §12 amendment 1, and noise here.
    body = (C.strip_labels(C.M.edges(src)) + S.NEW_EDGES
            + drop(C.nodes(g), {"KEY", "TITLE"})
            + drop(C.overlays(g, bb), {"MOVENOTE"}))

    # Pass 1: place the keys roughly, then ask neato how big everything is.
    keys = [key_node(231, 104, ps), progress_key(1180, 373, ps), title_node(807, 1213, ps)]
    pass1 = os.path.join(BUILD, "_fig1d_pass1.dot")
    open(pass1, "w").write("\n".join(header + body + keys + ["}"]))
    geo = measure(pass1, ref=("HOOK", g["HOOK"][:2]))

    # Pass 2: pin the wide key by its lower-left corner, which is fig 1c's own
    # anchor for it, so widening it grows the table rather than the canvas.
    kx, ky, kw, kh = geo["KEY"]
    keys[0] = key_node(kw / 2, kh / 2, ps)

    dot = header + body + keys + badges(geo, ps) + ["}"]
    return render(dot, STEM1, "fig 1d")


# Badges straddle the box's top-right corner, pushed this far out along both
# axes so they clear the box's own text: the two device boxes carry a title
# that runs the full width, and a badge sitting square on the corner covers it.
OFF = 9.0


def badges(geo, ps):
    L = []
    for nid, n in S.ELEMENT_PHASE.items():
        assert nid in geo, f"fig 1c has no node {nid!r} to badge"
        x, y, w, h = geo[nid]
        L.append(badge(f"BADGE_{nid.replace(':', '_')}",
                       x + w / 2 + OFF, y + h / 2 + OFF, ps[n]))
    lx, ly, lw, lh = geo["LAB_CME"]
    L.append(channel_badge("BADGE_CMEBLOCK", lx + 46, ly - lh / 2 - 16, ps[2],
                           "three jump modules coexist and compose"))
    ox, oy, ow, oh = geo["LAB_ODE"]
    L.append(channel_badge("BADGE_ODEBLOCK", ox - 12, oy - oh / 2 - 16, ps[1],
                           "a module contributes to a state it does not own"))
    return L


# --- fig 2d: the phase list, in dependency order ----------------------------
COL, ROW = 268.0, 132.0          # grid pitch, in points
X0, Y0 = 210.0, 820.0            # top-left chip centre
CHIPW = 236                      # chip body width, so titles wrap predictably


def wrap(text, n=32, align="left"):
    """Break a title across lines without hyphenating, so chips stay one width."""
    lines, line = [], ""
    for word in text.split():
        if line and len(line) + 1 + len(word) > n:
            lines.append(line)
            line = word
        else:
            line = f"{line} {word}".strip()
    lines.append(line)
    return f'<br align="{align}"/>'.join(lines) + f'<br align="{align}"/>'



def pos(n):
    c, r = S.GRID[n]
    return X0 + c * COL, Y0 - r * ROW


def chip(ph):
    x, y = pos(ph.n)
    fill, line, text = ((S.BUILT_FILL, S.BUILT_LINE, S.BUILT_TEXT) if ph.built
                        else (S.TODO_FILL, S.TODO_LINE, S.TODO_TEXT))
    mark = "&#10004;" if ph.built else "&#9744;"
    head = (f'<b>{mark} phase {ph.n}</b>' +
            (f' &#160;·&#160; <b>#{ph.pr}</b>' if ph.built else ''))
    stamp = (f'merged {ph.merged} · {ph.done}/{ph.total} tasks' if ph.built else
             (f'{ph.total} tasks, none ticked' if ph.total else
              'tasks are written once phase 16 has run'))
    return (f'  "P{ph.n}" [pos="{x:.1f},{y:.1f}!", shape=plaintext, style="", label=<'
            f'<table border="1" color="{line}" cellpadding="0" cellspacing="0" '
            'bgcolor="white">'
            f'<tr><td bgcolor="{fill}" align="left" balign="left" cellpadding="4" width="{CHIPW}">'
            f'<font color="{text}" point-size="11">{head}</font></td></tr>'
            f'<tr><td align="left" balign="left" cellpadding="4" width="{CHIPW}"><font point-size="10">'
            f'<b>{wrap(ph.title)}</b></font></td></tr>'
            f'<tr><td align="left" balign="left" cellpadding="4" bgcolor="#fbfcfc" width="{CHIPW}">'
            f'<font point-size="8.5" color="#4d5b63">{wrap(S.DELIVERS[ph.n], 38)}</font>'
            f'<br align="left"/><font point-size="8" color="{S.TODO_TEXT}">{stamp}</font>'
            '<br align="left"/></td></tr></table>>];')


def dep_edges(ps):
    L = []
    for a, b, style, note in S.DEPS:
        live = ps[a].built and ps[b].built
        # The six edges converging on phase 13 span the whole grid; drawn
        # lighter so they read as "everything lands here" rather than as six
        # separate claims fighting the chips they pass.
        colour = S.BUILT_LINE if live else ("#c8d0d4" if b == 13 else "#a9b4ba")
        lab = (f', label=<<font point-size="8" color="{S.TODO_TEXT}"><i>{note}</i></font>>'
               if note else "")
        L.append(f'  "P{a}" -> "P{b}" [color="{colour}", style={style}, '
                 f'penwidth={1.6 if live else (0.9 if b == 13 else 1.1)}, arrowsize=0.7{lab}];')
    return L


def band_labels():
    L = []
    for r, (name, gloss) in S.ROW_BAND.items():
        x, y = X0 - COL * 0.72, Y0 - r * ROW
        L.append(f'  "BAND{r}" [pos="{x:.1f},{y:.1f}!", shape=plaintext, style="", label=<'
                 f'<font point-size="11" color="#4d5b63"><b>{wrap(name, 24, "right")}</b>'
                 f'</font><font point-size="8" color="{S.TODO_TEXT}">'
                 f'<i>{wrap(gloss, 30, "right")}</i></font>>];')
    return L


def fig2d(ps):
    built = [p for p in ps.values() if p.built]
    done = sum(p.done for p in ps.values())
    total = sum(p.total for p in ps.values())
    tx = X0 + 2.5 * COL
    head = (f'  "T2" [pos="{tx:.1f},{Y0 + 1.15 * ROW:.1f}!", shape=plaintext, style="", '
            'label=<<font point-size="20"><b>Core A′ — the eighteen phases, in the order '
            'spec §11 puts them</b></font><br/>'
            f'<font point-size="11">{len(built)} merged, {len(ps) - len(built)} to go · '
            f'{done} of {total} tasks ticked · as of {max(p.merged for p in built)} · '
            'the third line of each chip is what the phase puts on fig 1d</font>>];')
    fy = Y0 - 6.8 * ROW
    foot = (f'  "F2" [pos="{tx:.1f},{fy:.1f}!", shape=plaintext, style="", label=<'
            f'<table border="1" color="{S.BUILT_LINE}" cellpadding="6" cellspacing="0" '
            'bgcolor="white"><tr><td align="left"><font point-size="11">'
            '<b>Every merged phase is framework, and that is D0 on purpose</b> — the '
            'protocol changes first, then the kill risk on a two-module toy,<br align="left"/>'
            'before a line of syn3A biology. It is also why fig 1d looks untouched: the '
            'phases that draw <i>boxes</i> all lie ahead.<br align="left"/>'
            '<font point-size="9.5">Phases 15–17 put nothing on fig 1d at all. They are the '
            'layer the project exists for, and a state graph cannot show it.'
            '</font><br align="left"/></font></td></tr></table>>];')
    dot = ['digraph roadmap {',
           '  graph [bgcolor="white", margin=0.4, splines=true, overlap=true];',
           '  node [shape=plaintext, fontname="Helvetica"];',
           '  edge [fontname="Helvetica"];']
    dot += dep_edges(ps) + [chip(ps[n]) for n in sorted(ps)] + band_labels() + [head, foot]
    return render(dot + ['}'], STEM2, "fig 2d")


def main():
    os.makedirs(OUT, exist_ok=True)
    ps = S.phases()
    fig1d(ps)
    fig2d(ps)


if __name__ == "__main__":
    main()

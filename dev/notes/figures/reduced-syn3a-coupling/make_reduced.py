"""Fig 1r -- the Core A' reduction, drawn on fig 1's own canvas.

Builds fig 1's DOT from the sibling generator, freezes the geometry that
`overlap=false` produced, then rewrites nodes and edges under reduction_spec.py.
Freezing matters: neato's overlap removal is content-sensitive, so editing any
label without it shifts every node and the two figures no longer overlay.

    python3 make_reduced.py        # -> fig1r_state_graph_reduced.{pdf,png}
"""
import os, re, shutil, subprocess, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import reduction_spec as S

HERE = os.path.dirname(os.path.abspath(__file__))
FIG1 = os.path.join(HERE, "..", "minimal-cell-coupling")
BUILD = os.path.join(HERE, "build")
OUT = os.path.join(BUILD, "out")
STEM = "fig1r_state_graph_reduced"

C_ODE, C_CME, C_HOOK, C_CTR, C_CL = "#1b4f72", "#7d3c98", "#117864", "#b9770e", "#c0392b"
C_RB, C_GR = "#1e8449", "#935116"
_geom = {}


# --- stage 1: rebuild fig 1's DOT, and check it is the shipped figure -------
def baseline():
    os.makedirs(OUT, exist_ok=True)
    shutil.copy(os.path.join(FIG1, "state_interfaces.csv"), OUT)
    subprocess.run([sys.executable, os.path.join(FIG1, "draw_state.py")],
                   cwd=BUILD, check=True, stdout=subprocess.DEVNULL)
    import hashlib
    h = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()
    assert h(os.path.join(OUT, "fig1_state_graph.png")) == \
           h(os.path.join(FIG1, "fig1_state_graph.png")), \
           "regenerated fig 1 differs from the shipped one -- the layout has moved"
    return open(os.path.join(OUT, "fig1_state_graph.dot")).read().split("\n")


# --- stage 2: the geometry overlap removal actually produced ---------------
def geometry():
    dot = os.path.join(OUT, "fig1_geom.dot")
    subprocess.run(["neato", "-n", "-Tdot", os.path.join(OUT, "fig1_state_graph.dot"),
                    "-o", dot], check=True)
    t = open(dot).read()
    g = {}
    for stmt in t.split("];"):
        head = stmt.split("[")[0]
        m = re.search(r'"?([A-Za-z_:]+)"?\s*$', head.strip())
        if not m or "->" in head:
            continue
        p = re.search(r'\bpos="([-\d.,e]+)"', stmt)
        w = re.search(r'\bwidth=([\d.]+)', stmt)
        h = re.search(r'\bheight=([\d.]+)', stmt)
        if p and w and h:
            x, y = (float(v) for v in p.group(1).split(","))
            g[m.group(1)] = (x, y, float(w.group(1)) * 72, float(h.group(1)) * 72)
    bb = [float(v) for v in re.search(r'bb="([^"]*)"', t).group(1).split(",")]
    return g, bb


# --- stage 3: node statements, authored rather than patched ----------------
def nodes(g):
    P = lambda n: f'pos="{g[n][0]},{g[n][1]}!"'
    L = []
    for n, lab in S.REMOVED_NODES.items():
        L.append(f'  "{n}" [{P(n)}, label=<<font color="{S.GREY_TEXT}">{lab}</font>>, '
                 f'fillcolor="{S.GREY_FILL}", color="{S.GREY_LINE}"];')
    for n, lab in S.CHANGED_NODES.items():
        border = C_CME if n.startswith("CME:") else C_ODE
        L.append(f'  "{n}" [{P(n)}, label=<{lab}>, '
                 f'fillcolor="{S.CHANGED_FILL}", color="{border}"];')

    L.append(f'  "CURRENCY" [{P("CURRENCY")}, style="filled,dashed", fillcolor="#eef2f5", '
             f'color="#95a5a6", fontcolor="#566573", label=<<b>currency pool</b>  '
             '<font point-size="8.5">ATP · ADP · AMP · P<sub>i</sub> · PP<sub>i</sub> · '
             'NAD(P)(H)</font><br/>'
             '<font point-size="7.5"><i>not a module — a drawing device.</i> Now shared by '
             '<b>three</b> modules —<br/>Transport, Central, Nucleotide. NAD(P)(H) survives '
             'through LDH_L alone.</font>>];')
    L.append(f'  "HOOK" [{P("HOOK")}, fillcolor="#d4efdf", color="{C_HOOK}", penwidth=2.2, '
             'label=<<b>the hook</b> — state<br/><font point-size="9.5">every <b>1.0 s</b> '
             'of cell time</font><br/> <br/>'
             '<font point-size="8.5">counts → ODE initial conditions (partTomM)<br/>'
             'integrate 1 s (LSODA)<br/>write counts back (mMtoPart)</font><br/> <br/>'
             '<font point-size="8"><i>~6,300 calls per cell cycle<br/>'
             'gradients do not cross here</i></font>>];')
    L.append(f'  "REBUILD" [{P("REBUILD")}, fillcolor="{S.CHANGED_FILL}", color="{C_RB}", '
             'penwidth=2.2, label=<<b>the CME rebuild</b> — rate constants<br/>'
             '<font point-size="9.5">every <b>60 s</b>, MinCell_restart.py</font><br/> <br/>'
             '<font point-size="8.5">live pools → k_transcription, k_translation<br/>'
             'the whole CME is reconstructed</font><br/> <br/>'
             f'<font point-size="8" color="{S.GREY_TEXT}"><i>k_replication branch gone with '
             'replication</i></font><br/>'
             '<font point-size="8"><i>still the only ODE→CME direction</i></font>>];')
    L.append(f'  "GROWTH" [{P("GROWTH")}, fillcolor="{S.CHANGED_FILL}", color="{C_GR}", '
             'penwidth=2.2, label=<<b>growth</b> — surface-area accounting<br/>'
             '<font point-size="9.5">every 1.0 s, in neither block</font><br/> <br/>'
             f'<font point-size="8.5" color="{S.GREY_TEXT}">'
             '<s>CellSA_Lip = 0.513 × Σ (lipid count × headgroup area)</s></font><br/>'
             '<font point-size="8.5">CellSA_Prot = 28.0 nm² × (membrane protein count)<br/>'
             'r = √(CellSA/4π) &#160;&#160; V = (4/3)πr³</font><br/> <br/>'
             '<font point-size="8"><i>ptsG is the only membrane protein: 831 × 28 nm².<br/>'
             '<b>V reaches ~1.07× — the 2× clamp is never approached.</b><br/>'
             'Report fractional growth, never a doubling time.</i></font>>];')
    L.append(f'  "Medium" [{P("Medium")}, style="filled,dashed", fillcolor="#fdedec", '
             f'color="{C_CL}", fontcolor="{C_CL}", label=<<b>growth medium</b><br/>'
             '<font point-size="8.5"><b>17 → 5 clamped concentrations</b></font><br/>'
             '<font point-size="7.5">glucose_e (40 mM), CTP, UTP, aa pool, O₂</font>>];')

    L.append(f'  "LAB_CME" [{P("LAB_CME")}, shape=plaintext, style="", label=<'
             f'<font color="{C_CME}" point-size="15"><b>CME block</b></font><br/>'
             f'<font color="{C_CME}" point-size="10">Gillespie / Lattice Microbes · '
             'stochastic · not differentiable</font>>];')
    L.append(f'  "LAB_ODE" [{P("LAB_ODE")}, shape=plaintext, style="", label=<'
             f'<font color="{C_ODE}" point-size="15"><b>ODE block</b></font><br/>'
             f'<font color="{C_ODE}" point-size="10">odecell → LSODA · deterministic · '
             '<b>190 → 21 reactions</b></font><br/>'
             f'<font color="{C_ODE}" point-size="8.5"><i>three modules of six survive; '
             '~32 dynamic states and 5 chemostats</i></font>>];')
    L.append(f'  "TITLE" [{P("TITLE")}, shape=plaintext, style="", label=<'
             '<font point-size="21"><b>State coupling in reduced syn3A — Core A′</b></font><br/>'
             '<font point-size="11">the same graph as fig 1, with everything the Step 1 '
             'reduction deletes greyed and crossed · 21 metabolic reactions · 17 genes · '
             'positions identical to fig 1, so the two overlay</font>>];')
    L.append(f'  "KEY" [{P("KEY")}, shape=plaintext, style="", label=<'
             '<table border="1" color="#bdc3c7" cellpadding="5" cellspacing="0" bgcolor="white">'
             '<tr><td colspan="2" align="center">'
             '<b>edge kinds</b> — what actually crosses each interface</td></tr>'
             '<tr><td align="right"><font color="#1b4f72"><b>———</b></font></td><td align="left">'
             '<b>mass</b> · shared state, continuous; gradients cross inside the ODE block</td></tr>'
             '<tr><td align="right"><font color="#95a5a6"><b>- - -</b></font></td><td align="left">'
             '<b>currency</b> · routed via the pool node</td></tr>'
             '<tr><td align="right"><font color="#b9770e"><b>- - -</b></font></td><td align="left">'
             '<b>deferred counter</b> · CME accrues a cost, the hook debits it one step later</td></tr>'
             '<tr><td align="right"><font color="#7d3c98"><b>· · · ·</b></font></td><td align="left">'
             '<b>catalytic</b> · counts enter rate laws as parameters, no mass flows</td></tr>'
             '<tr><td align="right"><font color="#1e8449"><b>- - -</b></font></td><td align="left">'
             '<b>rate constant</b> · pools re-enter the CME at the 60 s rebuild, piecewise-constant</td></tr>'
             '<tr><td align="right"><font color="#935116"><b>———</b></font></td><td align="left">'
             '<b>volume</b> · counts set surface area, which sets V, which rescales concentrations</td></tr>'
             '<tr><td align="right"><font color="#c0392b"><b>- - ⊣</b></font></td><td align="left">'
             '<b>clamped</b> · a real dependence replaced by a constant</td></tr>'
             '</table>>];')
    return L


# --- stage 4: edge statements, taken from fig 1 and transformed -------------
HEX = re.compile(r'#[0-9a-fA-F]{6}')


def key_of(line):
    m = re.match(r'\s*"([^"]+)" -> "([^"]+)" \[', line)
    return (m.group(1), m.group(2)) if m else None


def match(line, spec_key):
    a, b, disc = spec_key
    k = key_of(line)
    return k == (a, b) and (disc in line if disc else True)


def edges(src):
    dead = {tuple(k) for k in S.DEAD_ENDPOINT} | {tuple(k) for k, _ in S.DEAD_COUPLING}
    hit = {k: 0 for k in dead} | {k: 0 for k in S.RELABEL}
    L = []
    for line in src:
        if not key_of(line):
            continue
        for k in S.RELABEL:
            if match(line, k):
                for old, new in S.RELABEL[k]:
                    assert old in line, f"relabel target missing: {k} {old!r}"
                    line = line.replace(old, new)
                hit[k] += 1
        for k in dead:
            if match(line, k):
                line = HEX.sub(S.GREY_EDGE, line)
                hit[k] += 1
                break
        L.append(line)
    bad = {k: v for k, v in hit.items() if v != 1}
    assert not bad, f"spec keys that did not match exactly one edge: {bad}"
    return L


# --- stage 5: overlays -----------------------------------------------------
def overlays(g, bb):
    L = []
    # The cross is a glyph node, not two edges: `outputorder=edgesfirst` paints
    # every edge beneath every node, so an edge-drawn cross would sit under the box.
    # Sized to the box it cancels, and slightly transparent so the name stays legible.
    for i, n in enumerate(S.REMOVED_NODES):
        x, y, _, h = g[n]
        L.append(f'  "XX{i}" [pos="{x},{y}!", shape=plaintext, style="", '
                 f'label=<<font color="{S.CROSS}d9" point-size="{h*1.30:.0f}">'
                 '<b>\u2717</b></font>>];')

    L.append('  "REDKEY" [pos="551,966!", shape=plaintext, style="", label=<'
             '<table border="1" color="#bdc3c7" cellpadding="5" cellspacing="0" bgcolor="white">'
             '<tr><td colspan="2" align="center"><b>the Core A′ reduction</b></td></tr>'
             f'<tr><td align="right"><font color="{S.CROSS}" point-size="13"><b>✗</b></font>'
             f'</td><td align="left"><font point-size="9"><font color="{S.GREY_TEXT}"><b>removed</b></font> · the module is cut, and every edge<br/>touching it is greyed</font></td></tr>'
             f'<tr><td align="right" bgcolor="{S.CHANGED_FILL}">&#160;&#160;&#160;</td>'
             '<td align="left"><font point-size="9"><b>changed</b> · the module survives in reduced form;<br/>the box carries <b>old → new</b></font></td></tr>'
             f'<tr><td align="right"><font color="{S.GREY_TEXT}"><b>———</b></font></td>'
             '<td align="left"><font point-size="9"><b>greyed edge</b> · the coupling is gone.<br/>Everything still coloured is live in Core A′</font></td></tr>'
             '</table>>];')

    L.append('  "CATNOTE" [pos="757,516!", shape=plaintext, style="", label=<'
             '<table border="1" color="#d5b8e0" cellpadding="4" cellspacing="0" bgcolor="white">'
             f'<tr><td align="left"><font color="{C_CME}" point-size="8">'
             '<b>the catalytic channel is not cut.</b> Fig 1 draws it once, into<br/>'
             'Cofactor, which Core A′ deletes — hence the grey edge above.<br/>'
             'In Core A′ it runs into <b>Transport, Central and Nucleotide</b>:<br/>'
             '17 single-gene GPRs, one-to-one, no AND/OR logic needed.</font></td></tr>'
             '</table>>];')

    L.append(f'  "ANCHOR_LL" [pos="{bb[0]+0.36},{bb[1]+0.36}!", shape=point, style=invis, '
             'width=0.01, height=0.01, label=""];')
    L.append(f'  "ANCHOR_UR" [pos="{bb[2]-0.36},{bb[3]-0.36}!", shape=point, style=invis, '
             'width=0.01, height=0.01, label=""];')
    return L


# --- stage 6: prove the two figures overlay ---------------------------------
def verify():
    """Rasterise both PDFs and check they land on the same pixels.

    Node centres are frozen, so everything fig 1 and fig 1r share must be
    pixel-identical. Boxes that changed grow about their fixed centre, and the
    edges attached to them shift by a few points -- that is the whole of the
    difference, and it is confined to the modules the reduction touches.
    """
    ref = os.path.join(FIG1, "fig1_state_graph.pdf")
    a, b = os.path.join(BUILD, "_fig1.png"), os.path.join(BUILD, "_fig1r.png")
    for pdf, png in ((ref, a), (os.path.join(HERE, STEM + ".pdf"), b)):
        subprocess.run(["pdftoppm", "-r", "150", "-png", "-singlefile", pdf,
                        png[:-4]], check=True)
    dim = lambda p: subprocess.run(["magick", "identify", "-format", "%wx%h", p],
                                   capture_output=True, text=True, check=True).stdout
    size = dim(a)
    assert size == dim(b), f"rasters differ in size: {size} vs {dim(b)}"

    # the hook is untouched by the reduction: it must be identical, pixel for pixel
    px = int(size.split("x")[0]) / (1349.6 + 2 * 28.8)     # px per graphviz point
    gx, gy, gw, gh = _geom["HOOK"]
    x = int((gx + 28.8 - gw / 2 + 10) * px); w = int((gw - 20) * px)
    y = int(int(size.split("x")[1]) - (gy + 28.8 + gh / 2 - 10) * px); h = int((gh - 20) * px)
    crop = f"{w}x{h}+{x}+{y}"
    d = subprocess.run(["magick", "(", a, "-crop", crop, "+repage", "-colorspace", "Gray", ")",
                        "(", b, "-crop", crop, "+repage", "-colorspace", "Gray", ")",
                        "-compose", "difference", "-composite",
                        "-format", "%[fx:maxima]", "info:"],
                       capture_output=True, text=True, check=True).stdout.strip()
    assert float(d) == 0.0, f"the hook box is not pixel-identical (max diff {d})"

    # red = fig 1, cyan = fig 1r; anything black sits in both
    ov = os.path.join(BUILD, "overlay_check.png")
    subprocess.run(["magick", "-size", size, "xc:white",
                    "(", a, "-colorspace", "Gray", "-level", "0%,60%", ")",
                    "-compose", "CopyRed", "-composite",
                    "(", b, "-colorspace", "Gray", "-level", "0%,60%", ")",
                    "-compose", "CopyGreen", "-composite",
                    "(", b, "-colorspace", "Gray", "-level", "0%,60%", ")",
                    "-compose", "CopyBlue", "-composite",
                    "-resize", "2400x", ov], check=True)
    for f in (a, b):
        os.remove(f)
    print(f"ok  rasters {size} identical; hook pixel-exact; wrote {ov}")


def main():
    src = baseline()
    g, bb = geometry()
    _geom.update(g)
    header = [l.replace("overlap=false", "overlap=true") for l in src[:5]]
    dot = header + edges(src) + nodes(g) + overlays(g, bb) + ["}"]
    path = os.path.join(BUILD, STEM + ".dot")
    open(path, "w").write("\n".join(dot))
    for fmt in ("pdf", "png"):
        subprocess.run(["neato", "-n", f"-T{fmt}", "-Gdpi=150", path,
                        "-o", os.path.join(HERE, f"{STEM}.{fmt}")], check=True)
    size = lambda p: subprocess.run(["pdfinfo", p], capture_output=True, text=True
                                    ).stdout.split("Page size:")[1].split("\n")[0].strip()
    ref, new = size(os.path.join(FIG1, "fig1_state_graph.pdf")), size(os.path.join(HERE, STEM + ".pdf"))
    assert ref == new, f"canvas differs -- fig 1 {ref}, reduced {new}"
    print(f"ok  canvas {new}  ({len(dot)} statements)")
    if "--verify" in sys.argv:
        verify()


if __name__ == "__main__":
    main()

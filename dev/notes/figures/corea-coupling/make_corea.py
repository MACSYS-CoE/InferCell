"""Fig 1c -- Core A′ with tRNA charging in the ODE block, on fig 1's canvas.

Fig 1r (../reduced-syn3a-coupling) draws the lumped charging step inside the
CME, where the wave plan put it. The scoping note, the frozen registry and
spec/spec.md §12 amendment 1 all place it in the ODE block. This figure is fig 1r
with that one module moved and its edges redrawn by kind. Every other node sits
where fig 1 puts it.

    python3 make_corea.py        # -> fig1c_state_graph_corea.{pdf,png}

Stages 2 and 4 -- freezing the geometry overlap removal produced, and
transforming fig 1's edge statements -- are fig 1r's own functions, imported and
pointed at this directory's build/ and this directory's spec. Stage 1 is
re-implemented here without fig 1r's byte-identity assertion: the shipped fig 1
raster was made with a different Graphviz and font stack, and on the cluster
(Graphviz 2.44) the rebuilt fig 1 renders at a different canvas size, so the
assertion cannot pass and fig 1r itself does not build here. The DOT this figure
freezes is fig 1's pinned layout after *this* Graphviz's overlap removal. A
mismatch is reported, not fatal. No raster verification: charging is meant to
move, so fig 1c does not overlay fig 1 and no check claims it does.
"""
import os, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
FIG1 = os.path.join(HERE, "..", "minimal-cell-coupling")
FIG1R = os.path.join(HERE, "..", "reduced-syn3a-coupling")
BUILD = os.path.join(HERE, "build")
OUT = os.path.join(BUILD, "out")
STEM = "fig1c_state_graph_corea"

sys.path.insert(0, HERE)
sys.path.insert(0, FIG1R)
import corea_spec as S
import make_reduced as M
M.BUILD, M.OUT, M.S = BUILD, OUT, S      # reuse baseline(), geometry(), edges()

C_ODE, C_CME, C_HOOK, C_CL = "#1b4f72", "#7d3c98", "#117864", "#c0392b"
C_RB, C_GR = "#1e8449", "#935116"


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
    cx, cy = S.charging_pos(g)
    L.append(f'  "{S.CHG}" [pos="{cx},{cy}!", label=<{S.CHARGING_LABEL}>, '
             f'fillcolor="{S.CHANGED_FILL}", color="{C_ODE}", penwidth=2.0];')

    L.append(f'  "CURRENCY" [{P("CURRENCY")}, style="filled,dashed", fillcolor="#eef2f5", '
             f'color="#95a5a6", fontcolor="#566573", label=<<b>currency pool</b>  '
             '<font point-size="8.5">ATP · ADP · AMP · P<sub>i</sub> · PP<sub>i</sub> · '
             'NAD(P)(H)</font><br/>'
             '<font point-size="7.5"><i>not a module — a drawing device.</i> Now shared by '
             '<b>four</b> modules —<br/>Transport, Central, Nucleotide, tRNA charging. '
             'NAD(P)(H) survives through LDH_L alone.</font>>];')
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
             'stochastic · not differentiable</font><br/>'
             f'<font color="{C_CME}" point-size="8.5"><i>three modules: 17 × 3 reactions '
             'plus one translocation</i></font>>];')
    L.append(f'  "LAB_ODE" [{P("LAB_ODE")}, shape=plaintext, style="", label=<'
             f'<font color="{C_ODE}" point-size="15"><b>ODE block</b></font><br/>'
             f'<font color="{C_ODE}" point-size="10">odecell → LSODA · deterministic · '
             '<b>190 → 22 reactions</b></font><br/>'
             f'<font color="{C_ODE}" point-size="8.5"><i>four modules: three of six survive, '
             'plus tRNA charging moved in from the CME; ~32 dynamic states and 5 chemostats</i></font>>];')
    L.append(f'  "TITLE" [{P("TITLE")}, shape=plaintext, style="", label=<'
             '<font point-size="21"><b>State coupling in reduced syn3A — Core A′, '
             'with tRNA charging in the ODE block</b></font><br/>'
             '<font point-size="11">fig 1r redrawn with the lumped charging step where the '
             'scoping note places it · 22 metabolic reactions · 17 genes · '
             'every other node sits where fig 1 puts it</font>>];')
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


def overlays(g, bb):
    L = []
    for i, n in enumerate(S.REMOVED_NODES):
        x, y, _, h = g[n]
        L.append(f'  "XX{i}" [pos="{x},{y}!", shape=plaintext, style="", '
                 f'label=<<font color="{S.CROSS}d9" point-size="{h*1.30:.0f}">'
                 '<b>✗</b></font>>];')

    rx, ry = g["CME:Replication"][0] + 160, g["Medium"][1] - 15   # the gap left of the medium
    L.append(f'  "REDKEY" [pos="{rx},{ry}!", shape=plaintext, style="", label=<'
             '<table border="1" color="#bdc3c7" cellpadding="5" cellspacing="0" bgcolor="white">'
             '<tr><td colspan="2" align="center"><b>the Core A′ reduction</b></td></tr>'
             f'<tr><td align="right"><font color="{S.CROSS}" point-size="13"><b>✗</b></font>'
             f'</td><td align="left"><font point-size="9"><font color="{S.GREY_TEXT}"><b>removed</b></font> · the module is cut, and every edge<br/>touching it is greyed</font></td></tr>'
             f'<tr><td align="right" bgcolor="{S.CHANGED_FILL}">&#160;&#160;&#160;</td>'
             '<td align="left"><font point-size="9"><b>changed</b> · the module survives in reduced form;<br/>the box carries <b>old → new</b></font></td></tr>'
             f'<tr><td align="right"><font color="{S.GREY_TEXT}"><b>———</b></font></td>'
             '<td align="left"><font point-size="9"><b>greyed edge</b> · the coupling is gone.<br/>Everything still coloured is live in Core A′</font></td></tr>'
             f'<tr><td align="right"><font color="{C_ODE}" point-size="11"><b>▭</b></font></td>'
             '<td align="left"><font point-size="9"><b>moved</b> · one module crosses from the CME to the ODE block;<br/>its box is the only one not at fig 1’s position</font></td></tr>'
             '</table>>];')

    hx, hy = g["HOOK"][0] + 50, g["HOOK"][1] - 150
    L.append(f'  "CATNOTE" [pos="{hx},{hy}!", shape=plaintext, style="", label=<'
             '<table border="1" color="#d5b8e0" cellpadding="4" cellspacing="0" bgcolor="white">'
             f'<tr><td align="left"><font color="{C_CME}" point-size="8">'
             '<b>the catalytic channel is not cut.</b> Fig 1 draws it once, into<br/>'
             'Cofactor, which Core A′ deletes — hence the grey edge above.<br/>'
             'In Core A′ it runs into <b>Transport, Central and Nucleotide</b>:<br/>'
             '17 single-gene GPRs, one-to-one, no AND/OR logic needed.</font></td></tr>'
             '</table>>];')

    cx, cy = S.charging_pos(g)
    L.append(f'  "MOVENOTE" [pos="{cx},{cy-105}!", shape=plaintext, style="", label=<'
             '<table border="1" color="#aec6d8" cellpadding="4" cellspacing="0" bgcolor="white">'
             f'<tr><td align="left"><font color="{C_ODE}" point-size="8">'
             '<b>why charging moved.</b> Fig 1 draws it in the CME because the published<br/>'
             'model runs it there. The scoping note always placed the lumped step in the<br/>'
             'ODE block, and the frozen registry marks both tRNA species metabolite-scale;<br/>'
             'the wave plan mis-assigned it and fig 1r inherited that. Spec §12 amendment 1.<br/>'
             'Its 553 events/s were 353× every other stochastic event combined.</font></td></tr>'
             '</table>>];')

    L.append(f'  "ANCHOR_LL" [pos="{bb[0]+0.36},{bb[1]+0.36}!", shape=point, style=invis, '
             'width=0.01, height=0.01, label=""];')
    L.append(f'  "ANCHOR_UR" [pos="{bb[2]-0.36},{bb[3]-0.36}!", shape=point, style=invis, '
             'width=0.01, height=0.01, label=""];')
    return L


def baseline():
    """Fig 1r's stage 1, with the shipped-raster assertion downgraded to a report."""
    import hashlib, shutil
    os.makedirs(OUT, exist_ok=True)
    shutil.copy(os.path.join(FIG1, "state_interfaces.csv"), OUT)
    subprocess.run([sys.executable, os.path.join(FIG1, "draw_state.py")],
                   cwd=BUILD, check=True, stdout=subprocess.DEVNULL)
    h = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()
    if h(os.path.join(OUT, "fig1_state_graph.png")) != h(os.path.join(FIG1, "fig1_state_graph.png")):
        print("note  rebuilt fig 1 is not byte-identical to the shipped raster; "
              "positions are fig 1's pins after this Graphviz's overlap removal")
    return open(os.path.join(OUT, "fig1_state_graph.dot")).read().split("\n")


def strip_labels(lines):
    """Dead edges whose fig 1 labels would pile up where charging used to sit."""
    import re
    out = []
    for line in lines:
        if any(M.match(line, k) for k in S.STRIP_LABELS):
            line = re.sub(r', label=<.*>\]', "]", line)
        out.append(line)
    return out


def main():
    src = baseline()
    g, bb = M.geometry()
    header = [l.replace("overlap=false", "overlap=true") for l in src[:5]]
    dot = header + strip_labels(M.edges(src)) + S.NEW_EDGES + nodes(g) + overlays(g, bb) + ["}"]
    path = os.path.join(BUILD, STEM + ".dot")
    open(path, "w").write("\n".join(dot))
    for fmt in ("pdf", "png"):
        subprocess.run(["neato", "-n", f"-T{fmt}", "-Gdpi=150", path,
                        "-o", os.path.join(HERE, f"{STEM}.{fmt}")], check=True)
    size = lambda p: subprocess.run(["pdfinfo", p], capture_output=True, text=True
                                    ).stdout.split("Page size:")[1].split("\n")[0].strip()
    ref, new = size(os.path.join(FIG1, "fig1_state_graph.pdf")), size(os.path.join(HERE, STEM + ".pdf"))
    print(f"ok  canvas {new} (shipped fig 1: {ref})  ({len(dot)} statements)")


if __name__ == "__main__":
    main()

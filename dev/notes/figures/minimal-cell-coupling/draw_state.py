"""Fig 1 -- state coupling graph of the CME-ODE well-stirred minimal cell.
Edges derived from out/state_interfaces.csv; node positions hand-set (neato -n)."""
import csv, collections, subprocess, html

rows = list(csv.DictReader(open("out/state_interfaces.csv")))
ODE = ["Transport", "Central", "Nucleotide", "Lipid", "Cofactor", "AminoAcid"]
CME = ["CME:Replication", "CME:Transcription", "CME:mRNAdecay",
       "CME:Translation", "CME:tRNAcharging"]
NRXN = {"Central": 26, "Nucleotide": 46, "Transport": 87, "Lipid": 17,
        "Cofactor": 13, "AminoAcid": 1}
NAME = {"CME:Replication": "DNA replication", "CME:Transcription": "transcription",
        "CME:mRNAdecay": "mRNA decay", "CME:Translation": "translation",
        "CME:tRNAcharging": "tRNA charging", "AminoAcid": "Amino acid"}
POS = {                                                     # inches
    "CME:Replication": (5.9, 13.7), "CME:mRNAdecay":    (5.9, 11.6),
    "CME:Transcription": (2.2, 11.0), "CME:Translation": (2.2, 8.2),
    "CME:tRNAcharging": (4.9, 4.9),
    "HOOK": (9.4, 9.6), "REBUILD": (6.4, 3.4), "GROWTH": (12.6, 1.4),
    "Transport": (13.8, 12.8), "Central": (19.1, 12.8),
    "Nucleotide": (13.8, 9.6), "Lipid": (19.1, 9.6),
    "Cofactor": (13.8, 6.4),   "AminoAcid": (19.1, 6.4),
    "CURRENCY": (18.8, 3.2),   "Medium": (12.4, 15.1),
    "KEY": (3.0, 1.1), "TITLE": (11.0, 16.5),
    "LAB_CME": (2.6, 15.4), "LAB_ODE": (16.4, 14.8),
}
SHORT = lambda s: s.replace("M_", "").replace("_c", "").replace("__L", "")
C_ODE, C_CME, C_HOOK, C_CTR, C_CAT, C_CL = ("#1b4f72", "#7d3c98", "#117864",
                                            "#b9770e", "#7d3c98", "#c0392b")
C_RB, C_GR = "#1e8449", "#935116"          # 60 s rebuild channel; growth block
def agg(k):
    d = collections.defaultdict(list)
    for r in rows:
        if r["kind"] == k: d[(r["source"], r["target"])].append(r["species"])
    return d
mass, currency = agg("mass"), agg("currency-mass")
counter = {**agg("counter-debit"), **agg("counter-credit")}
ratek, volume, derived = agg("rate-constant"), agg("volume"), agg("derived")
def lab(sps, cap=2):
    s = sorted({x for g in sps for x in g.split(",")})
    return html.escape(", ".join(SHORT(x) for x in s[:cap])
                       + (f" +{len(s)-cap}" if len(s) > cap else ""))
def P(n): return f'{POS[n][0]*72:.0f},{POS[n][1]*72:.0f}!'

L = []; A = L.append
A('digraph G {')
A('  splines=curved; overlap=false; outputorder=edgesfirst; pad=0.4; bgcolor="white";')
A('  graph [fontname="Helvetica"];')
A('  node [fontname="Helvetica", fontsize=11, shape=box, style="rounded,filled", penwidth=1.3, margin="0.17,0.10"];')
A('  edge [fontname="Helvetica", fontsize=8.5, arrowsize=0.7];')

for m in CME:
    A(f'  "{m}" [pos="{P(m)}", label=<<b>{NAME[m]}</b>>, fillcolor="#ead6f0", color="{C_CME}"];')
for m in ODE:
    A(f'  "{m}" [pos="{P(m)}", label=<<b>{NAME.get(m, m)}</b><br/>'
      f'<font point-size="8.5">{NRXN[m]} rxns</font>>, fillcolor="#d6e4f0", color="{C_ODE}"];')
A(f'  "CURRENCY" [pos="{P("CURRENCY")}", style="filled,dashed", fillcolor="#eef2f5", '
  f'color="#95a5a6", fontcolor="#566573", label=<<b>currency pool</b>  '
  '<font point-size="8.5">ATP · ADP · AMP · P<sub>i</sub> · PP<sub>i</sub> · NAD(P)(H) · CoA</font><br/>'
  '<font point-size="7.5"><i>not a module — a drawing device.</i> Shared by <b>five</b> modules —<br/>'
  'Amino acid has none, its one reaction (FMETTRS) uses no currency.</font>>];')
A(f'  "HOOK" [pos="{P("HOOK")}", fillcolor="#d4efdf", color="{C_HOOK}", penwidth=2.2, '
  f'label=<<b>the hook</b> — state<br/><font point-size="9.5">every <b>1.0 s</b> of cell time</font><br/> <br/>'
  f'<font point-size="8.5">counts → ODE initial conditions (partTomM)<br/>integrate 1 s (LSODA)<br/>'
  f'write counts back (mMtoPart)</font><br/> <br/>'
  f'<font point-size="8"><i>~6,300 calls per cell cycle<br/>gradients do not cross here</i></font>>];')
A(f'  "REBUILD" [pos="{P("REBUILD")}", fillcolor="#d5f5e3", color="{C_RB}", penwidth=2.2, '
  f'label=<<b>the CME rebuild</b> — rate constants<br/>'
  f'<font point-size="9.5">every <b>60 s</b>, MinCell_restart.py</font><br/> <br/>'
  f'<font point-size="8.5">live pools → k_transcription, k_translation<br/>'
  f'the whole CME is reconstructed</font><br/> <br/>'
  f'<font point-size="8"><i>the ODE→CME direction.<br/>'
  f'piecewise-constant, not a continuous propensity</i></font>>];')
A(f'  "GROWTH" [pos="{P("GROWTH")}", fillcolor="#fdf2e9", color="{C_GR}", penwidth=2.2, '
  f'label=<<b>growth</b> — surface-area accounting<br/>'
  f'<font point-size="9.5">every 1.0 s, in neither block</font><br/> <br/>'
  f'<font point-size="8.5">CellSA_Lip = 0.513 × Σ (lipid count × headgroup area)<br/>'
  f'CellSA_Prot = 28.0 nm² × (membrane protein count)<br/>'
  f'r = √(CellSA/4π) &#160;&#160; V = (4/3)πr³</font><br/> <br/>'
  f'<font point-size="8"><i>54% protein is the calibrated <b>initial</b> area split.<br/>'
  f'Which term dominates the <b>increment</b> is not set by the code.</i></font>>];')
A(f'  "Medium" [pos="{P("Medium")}", style="filled,dashed", fillcolor="#fdedec", '
  f'color="{C_CL}", fontcolor="{C_CL}", label=<<b>growth medium</b><br/>'
  f'<font point-size="8.5">17 clamped external concentrations</font>>];')

# --- edges -----------------------------------------------------------------
for (a, b), sp in mass.items():
    if a in ODE and b in ODE:
        A(f'  "{a}" -> "{b}" [color="{C_ODE}", penwidth={0.9+0.20*len(set(sp)):.2f}, '
          f'label=<<font color="{C_ODE}">{lab(sp)}</font>>];')
for m in ODE:
    if any(m in e for e in currency):
        A(f'  "{m}" -> "CURRENCY" [color="#b3bfc4", penwidth=0.8, arrowhead=none, style=dashed];')
A(f'  "Medium" -> "Transport" [color="{C_CL}", style=dashed, penwidth=1.7, '
  f'label=<<font color="{C_CL}">no feedback on the medium</font>>];')

# 1 s: deferred cost counters, CME -> hook -> ODE pools
for a in sorted({a for a, _ in counter}):
    sp = sum((v for (x, _), v in counter.items() if x == a), [])
    A(f'  "{a}" -> "HOOK" [color="{C_CTR}", style=dashed, penwidth=1.5, '
      f'label=<<font color="{C_CTR}">{lab(sp, 4)}</font>>];')
A(f'  "HOOK" -> "Nucleotide" [color="{C_CTR}", style=dashed, penwidth=2.4, '
  f'headlabel=<<font color="{C_CTR}"><b>debit NTP / dNTP pools</b><br/>'
  f'clamped at zero, deficit carried<br/><i>a max(0,·) on the interface</i></font>>, '
  f'labeldistance=8.5, labelangle=-13];')
A(f'  "HOOK" -> "Nucleotide" [color="{C_CTR}", style=dotted, penwidth=1.6, '
  f'headlabel=<<font color="{C_CTR}"><b>credit NMP</b><br/>mRNA decay recycling,<br/>'
  f'<i>unclamped</i></font>>, labeldistance=5.0, labelangle=17];')
A(f'  "HOOK" -> "Central" [color="{C_CTR}", style=dashed, penwidth=2.0, '
  f'headlabel=<<font color="{C_CTR}"><b>debit ATP</b><br/>hydrolysis: trsc, translat,<br/>'
  f'mRNAdeg, DNArep, transloc</font>>, labeldistance=7.0, labelangle=-15];')
A(f'  "HOOK" -> "AminoAcid" [color="{C_CTR}", style=dashed, penwidth=1.6, '
  f'headlabel=<<font color="{C_CTR}"><b>debit fmet-tRNA</b><br/>the one aa-tRNA<br/>'
  f'reconciled at the hook</font>>, labeldistance=6.5, labelangle=14];')
# 1 s: protein counts -> enzyme concentrations (drawn once; all six modules)
A(f'  "CME:Translation" -> "HOOK" [color="{C_CAT}", style=dotted, penwidth=2.0];')
A(f'  "HOOK" -> "Cofactor" [color="{C_CAT}", style=dotted, penwidth=2.4, '
  f'label=<<font color="{C_CAT}"><b>protein counts → enzyme conc.</b><br/>'
  f'GPR: 87 single-gene · 4 OR (sum) · 5 AND (<b>min</b>)<br/>'
  f'<i>all six ODE modules; drawn once</i></font>>];')
# 1 s: live ODE -> CME mass
for (a, b), sp in mass.items():                      # every reactant of the charging chain
    if b != "CME:tRNAcharging": continue
    col = C_ODE if a in ODE else C_CME
    note = "<br/><i>synced counts, every 1 s</i>" if a in ODE else ""
    A(f'  "{a}" -> "CME:tRNAcharging" [color="{col}", penwidth=1.6, '
      f'label=<<font color="{col}"><b>{html.escape(sorted(sp)[0]) if a in ODE and "aa pools" in sp[0] else lab(sp, 3)}</b>{note}</font>>];')
A(f'  "AminoAcid" -> "CME:Translation" [color="{C_ODE}", penwidth=1.4, '
  f'label=<<font color="{C_ODE}">fmet-tRNA</font>>];')

# 60 s: the CME rebuild -- ODE pools re-enter as rate constants
for (a, b), sp in ratek.items():
    A(f'  "{a}" -> "REBUILD" [color="{C_RB}", style="dashed", penwidth=1.7, '
      f'label=<<font color="{C_RB}">{lab(sp, 4)}</font>>];')
A(f'  "REBUILD" -> "CME:Transcription" [color="{C_RB}", style="dashed", penwidth=2.4, '
  f'label=<<font color="{C_RB}"><b>k_transcription</b><br/>'
  f'per gene, from live NTP conc.</font>>];')
A(f'  "REBUILD" -> "CME:Translation" [color="{C_RB}", style="dashed", penwidth=2.4, '
  f'label=<<font color="{C_RB}"><b>k_translation</b><br/>'
  f'per gene, from live charged-tRNA counts</font>>];')
A(f'  "REBUILD" -> "CME:Replication" [color="{C_RB}", style="dashed", penwidth=2.0, '
  f'headlabel=<<font color="{C_RB}"><b>k_replication</b><br/>'
  f'elongation, from live dNTP conc.</font>>, labeldistance=7.5, labelangle=13];')

# growth: surface area in, volume out
for (a, b), sp in derived.items():
    A(f'  "{a}" -> "GROWTH" [color="{C_GR}", penwidth=1.6, '
      f'label=<<font color="{C_GR}">{html.escape(sp[0])}</font>>];')
A(f'  "GROWTH" -> "HOOK" [color="{C_GR}", penwidth=2.4, '
  f'label=<<font color="{C_GR}"><b>V rescales every count↔mM conversion</b><br/>'
  f'growth dilutes all six ODE modules; drawn once</font>>];')
A(f'  "GROWTH" -> "REBUILD" [color="{C_GR}", penwidth=1.6, '
  f'label=<<font color="{C_GR}">CellV → getConc</font>>];')
A(f'  "GROWTH" -> "GROWTH" [color="{C_CL}", style=dashed, penwidth=1.3, arrowhead=tee, '
  f'label=<<font color="{C_CL}">V clamped at 6.70e-17 L = 2× initial.<br/>'
  f'<b>There is no division event.</b></font>>];')

# intra-CME mass
for (a, b), sp in mass.items():
    if a in CME and b in CME:
        if b == "CME:tRNAcharging": continue          # already drawn above, with the ODE reactants
        if (a, b) == ("CME:Translation", "CME:Replication"):
            A(f'  "{a}" -> "{b}" [color="{C_CME}", penwidth=1.3, '
              f'taillabel=<<font color="{C_CME}"><b>DnaA</b><br/>'
              f'<i>initiates replication</i></font>>, labeldistance=6.0, labelangle=-18];')
        else:
            A(f'  "{a}" -> "{b}" [color="{C_CME}", penwidth=1.3, '
              f'label=<<font color="{C_CME}">{lab(sp)}</font>>];')
# what remains genuinely frozen through the rebuild
for n, t in (("CME:Transcription", "RnaPconc = 187, frozen"),
             ("CME:Translation", "ribosomeConc = 503, frozen"),
             ("CME:Replication", "DNApol3 = 35, frozen")):
    A(f'  "{n}" -> "{n}" [color="{C_CL}", style=dashed, penwidth=1.3, arrowhead=tee, '
      f'tailport=w, headport=w, label=<<font color="{C_CL}">{t}</font>>];')

# --- text ------------------------------------------------------------------
A(f'  "LAB_CME" [pos="{P("LAB_CME")}", shape=plaintext, style="", label=<'
  f'<font color="{C_CME}" point-size="15"><b>CME block</b></font><br/>'
  f'<font color="{C_CME}" point-size="10">Gillespie / Lattice Microbes · stochastic · not differentiable</font>>];')
A(f'  "LAB_ODE" [pos="{P("LAB_ODE")}", shape=plaintext, style="", label=<'
  f'<font color="{C_ODE}" point-size="15"><b>ODE block</b></font><br/>'
  f'<font color="{C_ODE}" point-size="10">odecell → LSODA · deterministic · 190 reactions</font><br/>'
  f'<font color="{C_ODE}" point-size="8.5"><i>differentiable in principle within the ODE subsystem; the shipped '
  f'SciPy/LSODA implementation is not AD-enabled</i></font>>];')
A(f'  "TITLE" [pos="{P("TITLE")}", shape=plaintext, style="", label=<'
  '<font point-size="21"><b>State coupling in the well-stirred minimal cell</b></font><br/>'
  '<font point-size="11">JCVI-syn3A CME–ODE model · Minimal_Cell@db048ac · '
  'every edge derived from SBtab stoichiometry, defMetRxns.py rate forms, and the production driver</font>>];')
A(f'  "KEY" [pos="{P("KEY")}", shape=plaintext, style="", label=<'
  '<table border="1" color="#bdc3c7" cellpadding="5" cellspacing="0" bgcolor="white"><tr><td colspan="2" align="center">'
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
A('}')
open("out/fig1_state_graph.dot", "w").write("\n".join(L))
for fmt in ("pdf", "png"):
    subprocess.run(["neato", "-n", f"-T{fmt}", "-Gdpi=150",
                    "out/fig1_state_graph.dot", "-o", f"out/fig1_state_graph.{fmt}"], check=True)
print("ok")

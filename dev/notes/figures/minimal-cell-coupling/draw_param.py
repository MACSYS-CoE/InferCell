"""Fig 2 -- parameter coupling. Counts from out/param_prior_dominance.csv and
out/mu0_cross_file.csv; structure from each SBtab file's `Parameter prior`
MatrixInfo column (the balancing dependency map)."""
import subprocess, pandas as pd, numpy as np

dom = pd.read_csv("out/param_prior_dominance.csv")
mu  = pd.read_csv("out/mu0_cross_file.csv")
RT  = 8.314e-3 * 310.15

def stat(mod, qt):
    r = dom[(dom.module == mod) & (dom.quantity == qt)]
    if r.empty: return None, None
    n = int(r.n.iloc[0]); a = r.at_prior.iloc[0]
    return n, (None if pd.isna(a) else int(a))
def pct(mod, qt):
    n, a = stat(mod, qt)
    return "—" if n is None or a is None else f"{n} · {100*a/n:.0f}% at prior"

MODS = ["Central+AA+Cof", "Nucleotide", "Lipid"]
QT = {"mu":"standard chemical potential", "KM":"Michaelis constant",
      "KV":"catalytic rate constant geometric mean", "c":"concentration",
      "u":"concentration of enzyme"}

BASIC, DERIV, FACTOR = "#1b4f72", "#7d3c98", "#117864"
BROKEN, GREY = "#c0392b", "#7f8c8d"
POS = {}
def N(name, x, y): POS[name] = (x, y); return name
def P(n): return f'{POS[n][0]*72:.0f},{POS[n][1]*72:.0f}!'

L = []; A = L.append
A('digraph G {')
A('  splines=curved; overlap=false; outputorder=edgesfirst; pad=0.45; bgcolor="white";')
A('  graph [fontname="Helvetica"];')
A('  node [fontname="Helvetica", fontsize=11, shape=box, style="rounded,filled", penwidth=1.3, margin="0.16,0.09"];')
A('  edge [fontname="Helvetica", fontsize=9, arrowsize=0.7];')

def basic(name, x, y, label):
    N(name, x, y)
    A(f'  "{name}" [pos="{P(name)}", fillcolor="#d6e4f0", color="{BASIC}", label=<{label}>];')
def derived(name, x, y, label):
    N(name, x, y)
    A(f'  "{name}" [pos="{P(name)}", style="rounded,filled,dashed", fillcolor="#f4e8f7", '
      f'color="{DERIV}", label=<{label}>];')
def factor(name, x, y, label):
    N(name, x, y)
    A(f'  "{name}" [pos="{P(name)}", shape=box, style="filled", fillcolor="#d4efdf", '
      f'color="{FACTOR}", penwidth=1.8, label=<{label}>];')
def plain(name, x, y, label, **kw):
    N(name, x, y)
    extra = "".join(f', {k}={v}' for k, v in kw.items())
    A(f'  "{name}" [pos="{P(name)}", shape=plaintext, style="", label=<{label}>{extra}];')

# ============================ PANEL A: within a module =====================
plain("A_TITLE", 6.2, 15.4,
      '<font point-size="15"><b>A · Inside one balanced module</b></font><br/>'
      '<font point-size="10">the dependency map is declared in each SBtab file\'s '
      '<i>Parameter prior</i> MatrixInfo column</font>')

basic("MU", 2.0, 14.4, '<b>μ°<sub>i</sub></b>  standard chemical potential<br/>'
      f'<font point-size="8.5">one per species · {stat("Central+AA+Cof","standard chemical potential")[0]}'
      f' / {stat("Nucleotide","standard chemical potential")[0]}'
      f' / {stat("Lipid","standard chemical potential")[0]} per module</font>')
factor("F1", 5.9, 14.4, '<b>N<sup>T</sup> · thermodynamics</b><br/>'
       '<font point-size="9">K<sub>eq,r</sub> = exp(−Σ<sub>i</sub> n<sub>ir</sub> μ°<sub>i</sub> / RT)</font>')
derived("KEQ", 9.9, 14.4, '<b>K<sub>eq,r</sub></b>  equilibrium constant<br/>'
        '<font point-size="8.5">derived · 0% at prior width, always</font>')
basic("KV", 2.0, 13.0, '<b>KV<sub>r</sub></b>  catalytic rate constant (geometric mean)<br/>'
      f'<font point-size="8.5">{pct("Central+AA+Cof", QT["KV"])} · {pct("Nucleotide", QT["KV"])} · {pct("Lipid", QT["KV"])}</font>')
basic("KM", 2.0, 11.6, '<b>KM<sub>r,i</sub></b>  Michaelis constants<br/>'
      f'<font point-size="8.5">{pct("Central+AA+Cof", QT["KM"])} · {pct("Nucleotide", QT["KM"])} · {pct("Lipid", QT["KM"])}</font>')
factor("F2", 5.9, 12.4, '<b>Haldane ± ½</b><br/>'
       '<font point-size="9">k<sub>cat</sub><sup>±</sup> = KV<sub>r</sub> · K<sub>eq,r</sub><sup>±½</sup> · Π KM<sup>∓½</sup></font>')
derived("KCF", 9.9, 13.0, '<b>k<sub>cat</sub><sup>f</sup></b>  substrate catalytic rate<br/>'
        '<font point-size="8.5">derived · 0% at prior width</font>')
derived("KCR", 9.9, 11.6, '<b>k<sub>cat</sub><sup>r</sup></b>  product catalytic rate<br/>'
        '<font point-size="8.5">derived · 0% at prior width</font>')
basic("CONC", 2.0, 9.0, '<b>c<sub>i</sub></b>  metabolite concentrations (initial state)<br/>'
      f'<font point-size="8.5">{pct("Central+AA+Cof", QT["c"])} · {pct("Nucleotide", QT["c"])} · {pct("Lipid", QT["c"])}</font>')
N("U", 2.0, 10.2)
A(f'  "U" [pos="{P("U")}", fillcolor="#fdf0e3", color="#b9770e", label=<'
  '<b>u<sub>r</sub></b>  enzyme concentration<br/>'
  '<font point-size="8.5">100% at prior width in every module — <b>and it never matters:</b><br/>'
  'overwritten every 1.0 s from CME protein counts via the GPR rule</font>>];')

derived("RATE", 9.9, 9.6, '<b>v<sub>r</sub>(t)</b>  reaction rate law<br/>'
        '<font point-size="8.5">reversible convenience kinetics — what the integrator sees</font>')
for a in ("KCF", "KCR", "KM", "CONC", "U"):
    A(f'  "{a}" -> "RATE" [color="#95a5a6", penwidth=1.1, style=dashed];')
for a, b in (("MU","F1"), ("F1","KEQ"), ("KV","F2"), ("KM","F2"), ("KEQ","F2")):
    col = FACTOR if b.startswith("F") else FACTOR
    A(f'  "{a}" -> "{b}" [color="{col}", penwidth=1.4];')
A(f'  "F2" -> "KCF" [color="{FACTOR}", penwidth=1.4];')
A(f'  "F2" -> "KCR" [color="{FACTOR}", penwidth=1.4];')

# ============================ PANEL B: between modules =====================
plain("B_TITLE", 18.4, 15.4,
      '<font point-size="15"><b>B · Between modules: the constraint that isn\'t there</b></font><br/>'
      '<font point-size="10">each file was balanced separately, so μ° for a shared compound '
      'was estimated more than once</font>')

def blockbox(name, x, y, title, mod, nrxn, nsp):
    N(name, x, y)
    A(f'  "{name}" [pos="{P(name)}", fillcolor="#eaf2f8", color="{BASIC}", penwidth=1.7, label=<'
      f'<table border="0" cellpadding="2" cellspacing="0">'
      f'<tr><td colspan="2"><b>{title}</b></td></tr>'
      f'<tr><td colspan="2"><font point-size="8.5">{nrxn} reactions · {nsp} species</font></td></tr>'
      f'<tr><td align="left"><font point-size="8.5">μ°</font></td><td align="left"><font point-size="8.5">{stat(mod, QT["mu"])[0]}</font></td></tr>'
      f'<tr><td align="left"><font point-size="8.5">KM</font></td><td align="left"><font point-size="8.5">{pct(mod, QT["KM"])}</font></td></tr>'
      f'<tr><td align="left"><font point-size="8.5">KV</font></td><td align="left"><font point-size="8.5">{pct(mod, QT["KV"])}</font></td></tr>'
      f'<tr><td align="left"><font point-size="8.5">c</font></td><td align="left"><font point-size="8.5">{pct(mod, QT["c"])}</font></td></tr>'
      f'</table>>];')

blockbox("B1", 15.0, 13.9, "Central + Amino acid + Cofactor", "Central+AA+Cof", 37, 55)
blockbox("B2", 21.8, 13.9, "Nucleotide", "Nucleotide", 45, 53)
blockbox("B3", 15.0, 10.6, "Lipid", "Lipid", 17, 32)
N("B4", 21.8, 10.6)
A(f'  "B4" [pos="{P("B4")}", style="rounded,filled,dashed", fillcolor="#fdedec", color="{BROKEN}", '
  f'fontcolor="{BROKEN}", label=<<b>Transport</b>  <font point-size="8.5">87 reactions</font><br/>'
  f'<font point-size="8.5"><b>no balanced parameter table at all</b><br/>'
  f'k<sub>cat</sub> and K<sub>M</sub> written straight into defMetRxns.py;<br/>'
  f'external concentrations are literals</font>>];')

pairs = {("Central+AA+Cof","Nucleotide"): ("B1","B2"),
         ("Central+AA+Cof","Lipid"):      ("B1","B3"),
         ("Nucleotide","Lipid"):          ("B2","B3")}
for (a, b), (na, nb) in pairs.items():
    g = mu[(mu.fileA == a) & (mu.fileB == b)]
    if g.empty:
        g = mu[(mu.fileA == b) & (mu.fileB == a)]
    med = g.delta.abs().median(); n = len(g)
    fac = np.exp(med / RT)
    fs = f"{fac:,.0f}×" if fac >= 10 else f"{fac:.1f}×"
    A(f'  "{na}" -> "{nb}" [color="{BROKEN}", style=dashed, penwidth=1.8, dir=both, '
      f'arrowhead=tee, arrowtail=tee, label=<<font color="{BROKEN}">'
      f'{n} shared compounds<br/>median |Δμ°| = {med:.1f} kJ/mol<br/>'
      f'<b>K<sub>eq</sub> differs by ~{fs}</b></font>>];')

plain("B_NOTE", 18.4, 9.0,
      f'<table border="1" color="{BROKEN}" cellpadding="6" cellspacing="0" bgcolor="#fdedec">'
      f'<tr><td align="left"><font color="{BROKEN}" point-size="10">'
      f'Wegscheider constraints tie K<sub>eq</sub> to μ° <b>within</b> each block — which is why the derived rows<br/>'
      f'never sit at prior width. Across blocks there is no such tie: the same compound carries a different<br/>'
      f'μ° in each file, so the assembled model is not thermodynamically consistent as a whole.<br/>'
      f'RT = {RT:.2f} kJ/mol at 37 °C, so a ~19 kJ/mol gap is three orders of magnitude in K<sub>eq</sub>.'
      f'</font></td></tr></table>')

# ============================ PANEL C: gene expression =====================
plain("C_TITLE", 11.6, 7.3,
      '<font point-size="15"><b>C · Gene expression: no factor structure at all</b></font><br/>'
      '<font point-size="10">metabolism at least went through balancing; this block did not</font>')

N("SCALARS", 3.6, 5.8)
A(f'  "SCALARS" [pos="{P("SCALARS")}", fillcolor="#fdedec", color="{BROKEN}", penwidth=1.7, label=<'
  '<b>~19 hard-coded global scalars</b><br/>'
  '<font point-size="8.5">rnaPolKcat · rnaPolKd · rnaPolK0 · rrnaPolKcat · trnaPolKcat<br/>'
  'krnadeg · ptnDegRate · riboKcat · riboK0 · riboKd<br/>'
  'RnaPconc=187 · ribosomeConc=503 · ctRNAconc=150<br/>'
  'ATPconc · UTPconc · CTPconc · GTPconc<br/><br/>'
  '<b>no prior, no posterior, no provenance</b></font>>];')
N("DATA", 3.6, 3.0)
A(f'  "DATA" [pos="{P("DATA")}", fillcolor="#eaf2f8", color="{BASIC}", label=<'
  '<b>fixed data inputs</b><br/><font point-size="8.5">'
  'syn3A.gb genome sequence<br/>proteomics counts (Breuer et al. 2019)</font>>];')
N("PLATE", 11.6, 4.6)
A(f'  "PLATE" [pos="{P("PLATE")}", fillcolor="#ead6f0", color="{DERIV}", penwidth=1.7, label=<'
  '<table border="0" cellpadding="3" cellspacing="0">'
  '<tr><td colspan="2"><b>455 mRNA-coding loci</b>  '
  '<font point-size="8.5">(+ rRNA and tRNA operons)</font></td></tr>'
  '<tr><td align="right"><font point-size="9">k<sub>transcription</sub></font></td>'
  '<td align="left"><font point-size="9">= f(seq) · min(rnaPolKcat · ptnCount<sub>g</sub>/180, 180)</font></td></tr>'
  '<tr><td align="right"><font point-size="9">k<sub>degradation</sub></font></td>'
  '<td align="left"><font point-size="9">= (18/452)·88 / n<sub>tot</sub>(g)  '
  '<font color="#c0392b">— length only</font></font></td></tr>'
  '<tr><td align="right"><font point-size="9">k<sub>translation</sub></font></td>'
  '<td align="left"><font point-size="9">= f(aa sequence, riboKcat, riboKd, ribosomeConc)</font></td></tr>'
  '</table>>];')
N("CNOTE", 19.4, 4.8)
A(f'  "CNOTE" [pos="{P("CNOTE")}", shape=plaintext, style="", label=<'
  f'<table border="1" color="{GREY}" cellpadding="6" cellspacing="0" bgcolor="white">'
  f'<tr><td align="left"><font point-size="10">'
  f'Every one of the 455 rate constants is a <b>deterministic</b> function of<br/>'
  f'sequence and a handful of shared scalars. The graph is a star, not a chain:<br/>'
  f'perturb one scalar and all 455 rates move together — so the parameters are<br/>'
  f'maximally coupled and individually unidentifiable by construction.<br/><br/>'
  f'<font color="{BROKEN}">riboKcat is defined twice (10 in MinCell_CMEODE.py, 12 in<br/>'
  f'translation_rate_start.py); the imported module wins, so 10 is dead code.</font>'
  f'</font></td></tr></table>>];')
A(f'  "SCALARS" -> "PLATE" [color="{BROKEN}", penwidth=2.6, '
  f'label=<<font color="{BROKEN}">shared by all 455</font>>];')
A(f'  "DATA" -> "PLATE" [color="{BASIC}", penwidth=1.6, '
  f'label=<<font color="{BASIC}">per-gene sequence and proteomics count</font>>];')

# ---------------------------------------------------------------- title/key
plain("TITLE", 11.8, 16.8,
      '<font point-size="21"><b>Parameter coupling in the well-stirred minimal cell</b></font><br/>'
      '<font point-size="11">Which parameters constrain which — a different graph from the state '
      'coupling, over the same modules.<br/>Counts computed over the reactions actually active in '
      'the model, not the whole reconstruction.</font>')
plain("KEY", 8.8, 8.9,
      '<table border="1" color="#bdc3c7" cellpadding="5" cellspacing="0" bgcolor="white">'
      '<tr><td colspan="2" align="center"><b>node kinds</b></td></tr>'
      f'<tr><td align="right"><font color="{BASIC}">▬</font></td><td align="left">'
      '<b>basic</b> · independently specified, has its own prior</td></tr>'
      f'<tr><td align="right"><font color="{FACTOR}">▬</font></td><td align="left">'
      '<b>factor</b> · deterministic map declared by the balancing matrices</td></tr>'
      f'<tr><td align="right"><font color="{DERIV}">▭</font></td><td align="left">'
      '<b>derived</b> · a function of basic quantities, never free</td></tr>'
      f'<tr><td align="right"><font color="{BROKEN}">▭</font></td><td align="left">'
      '<b>unconstrained</b> · no prior, no balancing, or a broken link</td></tr>'
      '</table>')
A('}')
open("out/fig2_param_graph.dot", "w").write("\n".join(L))
for fmt in ("pdf", "png"):
    subprocess.run(["neato", "-n", f"-T{fmt}", "-Gdpi=150",
                    "out/fig2_param_graph.dot", "-o", f"out/fig2_param_graph.{fmt}"], check=True)
print("ok")

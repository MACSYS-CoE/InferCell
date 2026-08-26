"""Parameter-coupling analysis restricted to each module's ACTIVE reactions."""
import sys, collections
sys.path.insert(0, ".")
import pandas as pd, numpy as np
from sbtab import read_sbtab, parse_formula
from modules import *

OUT = "out/"
mods = module_reactions()
CENT = mods["Central"] + mods["AminoAcid"] + mods["Cofactor"]
GROUPS = [("Central+AA+Cof", CENTRAL_FILE, CENT),
          ("Nucleotide",     NUCL_FILE,    mods["Nucleotide"]),
          ("Lipid",          LIPID_FILE,   mods["Lipid"])]

def species_of(path, rxnids):
    rt = read_sbtab(path)["Reaction"]
    sp, kept = set(), []
    for r in rt.itertuples():
        if r.Reaction in rxnids:
            s, p = parse_formula(r.ReactionFormula)
            sp |= {x[1] for x in s} | {x[1] for x in p}
            kept.append(r.Reaction)
    return sp, kept

def prior_map(tab):
    return {r.QuantityType: (r.PriorGeometricStd, r.Dependence)
            for r in tab["Parameter prior"].itertuples()}

rows, mu, spec = [], {}, {}
for lab, path, rxnids in GROUPS:
    tab = read_sbtab(path)
    sp, kept = species_of(path, rxnids)
    spec[lab] = sp
    par = tab["Parameter"].copy()
    par["gstd"] = pd.to_numeric(par["UnconstrainedGeometricStd"], errors="coerce")
    par["mean"] = pd.to_numeric(par["UnconstrainedMean"], errors="coerce")
    rid, cid = "Reaction:SBML:reaction:id", "Compound:SBML:species:id"
    keep = par[rid].isin(kept) | (par[rid].isin(["", "nan"]) & par[cid].isin(sp)) \
           | ((par[rid] == "") & par[cid].isin(sp))
    act = par[keep]
    pm = prior_map(tab)
    print(f"\n--- {lab}: {len(kept)} active reactions, {len(sp)} species, "
          f"{len(act)} active parameters ---")
    for qt, grp in act.groupby("QuantityType"):
        pg, dep = pm.get(qt, ("", "?"))
        try: pgv = float(pg)
        except ValueError: pgv = np.nan
        n = len(grp)
        at = int((abs(grp["gstd"] - pgv)/pgv < 1e-3).sum()) if pgv == pgv else -1
        f = "  --" if at < 0 else f"{100*at/n:5.0f}%"
        print(f"  {qt:<40} {dep:<8} n={n:<5} at prior={max(at,0):<5} {f}")
        rows.append(dict(module=lab, quantity=qt, dependence=dep, n=n,
                         at_prior=(at if at >= 0 else None)))
    m = act[act.QuantityType == "standard chemical potential"]
    mu[lab] = dict(zip(m[cid], m["mean"]))

pd.DataFrame(rows).to_csv(OUT + "param_prior_dominance.csv", index=False)

print("\n=== mu0 disagreement, restricted to compounds ACTIVE in both files ===")
out = []
labs = list(mu)
for i in range(len(labs)):
    for j in range(i+1, len(labs)):
        a, b = labs[i], labs[j]
        shared = sorted((spec[a] & spec[b]) & set(mu[a]) & set(mu[b]))
        d = [(s, mu[a][s], mu[b][s]) for s in shared
             if mu[a][s] == mu[a][s] and mu[b][s] == mu[b][s]]
        if not d: continue
        dv = np.array([abs(x-y) for _, x, y in d])
        print(f"  {a:>14} vs {b:<12} shared={len(d):<4} median|dmu0|={np.median(dv):6.2f} "
              f"max={dv.max():7.2f}  >1 kJ/mol: {(dv>1).sum()}/{len(d)}")
        for s, x, y in sorted(d, key=lambda t: -abs(t[1]-t[2]))[:5]:
            print(f"      {s:<16} {x:9.2f}  vs {y:9.2f}   d={x-y:8.2f}")
        out += [dict(fileA=a, fileB=b, compound=s, mu0_A=x, mu0_B=y, delta=x-y)
                for s, x, y in d]
pd.DataFrame(out).to_csv(OUT + "mu0_cross_file.csv", index=False)

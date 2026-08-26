"""Derive the state-coupling (species-flow) module graph of the CME-ODE
well-stirred minimal cell, Minimal_Cell@db048ac, CME_ODE/.

Output: state_interfaces.csv -- one row per (source, target, species, kind).
"""
import sys, csv, collections
sys.path.insert(0, ".")
from sbtab import read_sbtab, parse_formula
from modules import *
from handcoded import handcoded_reactions

OUT = "out/"

# ---------------------------------------------------------------- ODE modules
def ode_reactions():
    """rxn id -> (module, substrates, products, reversible)"""
    mods = module_reactions()
    rxn2mod = {r: m for m, rs in mods.items() for r in rs}
    R = {}
    for path in (CENTRAL_FILE, NUCL_FILE, LIPID_FILE):
        for row in read_sbtab(path)["Reaction"].itertuples():
            if row.Reaction in rxn2mod and row.Reaction not in R:
                s, p = parse_formula(row.ReactionFormula)
                R[row.Reaction] = (rxn2mod[row.Reaction], [x[1] for x in s],
                                   [x[1] for x in p],
                                   str(row.IsReversible).upper() == "TRUE")
    for row in read_sbtab(TRANS_FILE)["Reaction"].itertuples():
        if row.Reaction == "R_Ht":
            continue                                     # defMetRxns.py:1117
        s, p = parse_formula(row.ReactionFormula)
        R[row.Reaction] = ("Transport", [x[1] for x in s], [x[1] for x in p],
                           str(row.IsReversible).upper() == "TRUE")
    # hand-coded rate forms (defMetRxns.py:1869-2582)
    HC_MOD = {"PUNP5": "Nucleotide", "ATPase": "Central",
              "GHMT2": "Cofactor", "NADHK": "Cofactor"}
    HC_REV = {"PUNP5", "ATPase", "GHMT2", "NADHK"}
    for r in handcoded_reactions():
        mod = HC_MOD.get(r["id"], "Transport")
        R[r["id"]] = (mod, r["subs"], r["prods"], r["id"] in HC_REV)
        for _, const in r["clamped"]:                    # medium is a clamped source
            R[r["id"]] = (mod, r["subs"] + ["MEDIUM:" + const], r["prods"],
                          r["id"] in HC_REV)
    return R

CURRENCY = {"M_atp_c","M_adp_c","M_amp_c","M_pi_c","M_ppi_c","M_h2o_c","M_h_c",
            "M_nad_c","M_nadh_c","M_nadp_c","M_nadph_c","M_coa_c","M_accoa_c",
            "M_gtp_c","M_gdp_c"}

def ode_edges(R):
    prod = collections.defaultdict(set)   # species -> modules producing it
    cons = collections.defaultdict(set)   # species -> modules consuming it
    for rid, (mod, s, p, rev) in R.items():
        for x in s:
            cons[x].add(mod)
            if rev: prod[x].add(mod)
        for x in p:
            prod[x].add(mod)
            if rev: cons[x].add(mod)
    rows = []
    for sp in sorted(set(prod) | set(cons)):
        if sp.startswith("MEDIUM:"):
            for tgt in sorted(cons[sp]):
                rows.append(dict(source="Medium", target=tgt,
                                 species=sp.split(":")[1], kind="clamped-source",
                                 protocol="constant (growth medium)",
                                 differentiable="n/a"))
            continue
        for a in sorted(prod[sp]):
            for b in sorted(cons[sp]):
                if a == b: continue
                rows.append(dict(
                    source=a, target=b, species=sp,
                    kind="currency-mass" if sp in CURRENCY else "mass",
                    protocol="shared ODE state (continuous)",
                    differentiable="yes"))
    return rows

# ----------------------------------------------------------------- CME blocks
# Hand-encoded from source; every row carries its file:line provenance.
CME_ROWS = [
 # --- deferred cost counters: CME debits ODE pools at the hook write-back ---
 ("CME:Transcription","Nucleotide","M_atp_c,M_utp_c,M_ctp_c,M_gtp_c","counter-debit",
  "deferred 1 s; zero-clamped, deficit carried","no","in_out.py:180-259"),
 ("CME:Transcription","Central","M_atp_c","counter-debit",
  "ATP_trsc hydrolysis at hook","no","in_out.py:178,195"),
 ("CME:mRNAdecay","Nucleotide","M_amp_c,M_ump_c,M_cmp_c,M_gmp_c","counter-credit",
  "NMP recycling at hook (unclamped)","no","in_out.py:185,286"),
 ("CME:mRNAdecay","Central","M_atp_c","counter-debit",
  "ATP_mRNAdeg hydrolysis at hook","no","in_out.py:178"),
 ("CME:Translation","Nucleotide","M_gtp_c","counter-debit",
  "ATP_translat counter is paid in GTP","no","in_out.py:197-217"),
 ("CME:Translation","Central","M_atp_c","counter-debit",
  "ATP_transloc hydrolysis (membrane ptns)","no","in_out.py:178"),
 ("CME:Replication","Nucleotide","M_datp_c,M_dttp_c,M_dctp_c,M_dgtp_c","counter-debit",
  "deferred 1 s; zero-clamped","no","in_out.py:189,263"),
 ("CME:Replication","Central","M_atp_c","counter-debit",
  "ATP_DNArep hydrolysis at hook","no","in_out.py:178"),
 ("CME:Translation","AminoAcid","M_fmettrna_c","counter-debit",
  "only aa-tRNA reconciled at the hook","no","in_out.py:192,296"),
 # --- live ODE -> CME edges (mass action on synced counts) ---
 ("Transport","CME:tRNAcharging","20 aa pools","mass",
  "CME reactant, counts synced each 1 s","no","MinCell_CMEODE.py:1308"),
 ("Central","CME:tRNAcharging","M_atp_c","mass",
  "CME reactant, counts synced each 1 s","no","MinCell_CMEODE.py:1307"),
 ("AminoAcid","CME:Translation","M_fmettrna_c","mass",
  "R_FMETTRS is the one aa reaction left in the ODE","no","defMetRxns.py:288"),
 # --- intra-CME ---
 ("CME:tRNAcharging","CME:Translation","20 charged tRNAs","mass",
  "aa cost counters paid intra-CME, k=100/s","no","MinCell_CMEODE.py:1315"),
 ("CME:Transcription","CME:Translation","455 mRNAs","mass",
  "mRNA is the translation reactant","no","MinCell_CMEODE.py:730"),
 ("CME:Transcription","CME:mRNAdecay","455 mRNAs","mass",
  "mRNA is the decay reactant","no","MinCell_CMEODE.py:729"),
 ("CME:Replication","CME:Transcription","gene copy number","mass",
  "gene locus is the transcription reactant","no","rep_start.py"),
 # --- catalytic / parametric: protein counts become ODE enzyme concentrations ---
 *[("CME:Translation", m, "M_PTN_* counts", "catalytic",
    "enzyme conc. rebuilt from counts every 1 s", "no",
    "Simp.py:initModel / defMetRxns.py:1900+")
   for m in ("Central","Nucleotide","Transport","Lipid","Cofactor","AminoAcid")],
 # --- 60 s CME rebuild: ODE pools re-enter the CME as rate CONSTANTS ---
 # MinCell_restart.py loops once per biological minute and rebuilds the whole
 # CME. The hard-coded first-minute scalars of MinCell_CMEODE.py are commented
 # out there and replaced by live pool reads. This is the ODE -> CME direction.
 ("Nucleotide","CME:Transcription","M_atp_c,M_utp_c,M_ctp_c,M_gtp_c","rate-constant",
  "60 s rebuild; NTP conc. enters k_transcription","no","MinCell_restart.py:442-445,468,564"),
 ("CME:tRNAcharging","CME:Translation","20 charged tRNAs","rate-constant",
  "60 s rebuild; live counts enter NMonoSum","no","translation_rate_restart.py:89"),
 ("AminoAcid","CME:Translation","M_fmettrna_c","rate-constant",
  "60 s rebuild; live count enters k_translation","no","translation_rate_restart.py:113"),
 # --- growth: surface-area accounting, its own block ---
 ("Lipid","Growth","9 lipid headgroups","derived",
  "CellSA_Lip = 0.513 x sum(count x headgroup area)","yes","in_out.py:73-95"),
 ("CME:Translation","Growth","93 membrane protein loci","derived",
  "CellSA_Prot = 28.0 nm^2 x count; 94 species terms (0779 -> ptsg, ptsg_P)","no","in_out.py:31,48,58"),
 *[("Growth", m, "CellV", "volume",
    "V rescales every count<->mM conversion at the hook", "yes",
    "Rxns.partTomM:42 / in_out.mMtoPart:143")
   for m in ("Central","Nucleotide","Transport","Lipid","Cofactor","AminoAcid")],
 ("Growth","CME:Transcription","CellV","volume",
  "CellV feeds getConc at the 60 s rebuild","no","MinCell_restart.py:173-176"),
 ("Growth","CME:Translation","CellV","volume",
  "CellV feeds partTomM at the 60 s rebuild","no","translation_rate_restart.py:89,113"),
 # --- replication elongation is rebuilt from live dNTP pools, same as transcription ---
 ("Nucleotide","CME:Replication","M_datp_c,M_dttp_c,M_dctp_c,M_dgtp_c","rate-constant",
  "60 s rebuild; dNTP conc. enters replication elongation","no","rep_restart.py:239-243"),
 ("Growth","CME:Replication","CellV","volume",
  "CellV feeds getConc at the 60 s rebuild","no","rep_restart.py:239-243"),
 # --- intra-CME reactants of the charging chain ---
 ("CME:Transcription","CME:tRNAcharging","20 uncharged tRNAs","mass",
  "tRNA transcription produces the charging reactant; 29 loci","no",
  "MinCell_restart.py:1533,1566"),
 ("CME:Translation","CME:tRNAcharging","20 M_PTN_* synthetases","mass",
  "synthetase is an explicit charging reactant/product","no",
  "MinCell_restart.py:1502-1509,1519-1525"),
 # --- DnaA links translation to replication initiation ---
 ("CME:Translation","CME:Replication","M_DnaA_c","mass",
  "DnaA filament on oriC initiates replication","no",
  "protein_metabolites_frac.csv:7 / rep_restart.py:55-95"),
 # --- clamped: edges that exist physically but are constants in the code ---
 ("CME:Translation","CME:Translation","ribosomeConc","clamped",
  "503 ribosomes, frozen through the rebuild","severed","MinCell_restart.py:506"),
 ("CME:Transcription","CME:Transcription","RnaPconc","clamped",
  "187 RNAP, frozen through the rebuild","severed","MinCell_restart.py:462"),
 ("CME:Replication","CME:Replication","DNApol3","clamped",
  "35 DNA pol III, frozen through the rebuild","severed","rep_restart.py:244"),
 ("Growth","Growth","CellV","clamped",
  "V clamped at 6.70e-17 L = 2x initial; no division event","no","in_out.py:118-120"),
]

def main():
    R = ode_reactions()
    rows = ode_edges(R)
    for s, t, sp, k, pr, df, prov in CME_ROWS:
        rows.append(dict(source=s, target=t, species=sp, kind=k,
                         protocol=pr, differentiable=df, provenance=prov))
    for r in rows:
        r.setdefault("provenance", "SBtab ReactionFormula")
    with open(OUT + "state_interfaces.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["source","target","species","kind",
                                           "protocol","differentiable","provenance"])
        w.writeheader(); w.writerows(rows)

    print(f"ODE reactions in model: {len(R)}")
    bymod = collections.Counter(v[0] for v in R.values())
    for m, n in bymod.most_common(): print(f"  {m:<12} {n}")
    print(f"\ninterface rows: {len(rows)}")
    print(collections.Counter(r['kind'] for r in rows))
    agg = collections.Counter((r["source"], r["target"]) for r in rows)
    print(f"\ndistinct module->module edges: {len(agg)}")

main()

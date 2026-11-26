"""KIM-177: the non-closure inventory of the full syn3A cell, counted from the tables.

Reads the exported reference model (syn3a-ensemble `data/model/`) and two SpeedyCell
source files, and prints every count the note `README.md` quotes. Standard library only.

    /fred/oz022/tkimpson/miniconda3/bin/python3 -I nonclosure_inventory.py \
        [--model DIR] [--speedycell DIR] > nonclosure_inventory.out.txt

The hard asserts are the issue's verify-by clauses: the second-order families sum to
363, the replication reactions number 493, the count-dependent rate rules 1,434, and
every reaction falls in exactly one family.
"""
import argparse
import collections
import csv
import fractions
import hashlib
import json
import os
import re
import subprocess
import sys

ENSEMBLE = "/fred/oz022/tkimpson/syn3a-ensemble"
SPEEDYCELL_JL = "/fred/oz022/tkimpson/SpeedyCell.jl"
HERE = os.path.dirname(os.path.abspath(__file__))
REDUCED_GENES = os.path.join(HERE, "..", "..", "..", "src", "organisms", "coreA", "data", "transcription_genes.tsv")
NONSMOOTH = ("min", "max", "round", "select_gt", "floor", "trunc", "floordiv", "mod", "abs")


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_csv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def parse_stoich(s):
    out = []
    for item in s.split(";"):
        i, d = item.split(":")
        out.append((int(i), int(d)))
    return out


# ---------------------------------------------------------------------------
# Families. Every CME reaction gets exactly one, by an explicit rule.
# ---------------------------------------------------------------------------

AA_RE = r"M_(ala|arg|asn|asp|cys|gln|glu|gly|his|ile|leu|lys|met|phe|pro|ser|thr|trp|tyr|val)(__L)?_c"


def second_order_family(reac, prods):
    """(family, subfamily) of a type-2 reaction from its reactant and product names."""
    a, b = reac
    pos = {n for n, d in prods if d > 0}
    if "M_DnaA_c" in (a, b):
        other = b if a == "M_DnaA_c" else a
        if re.fullmatch(r"High_Affinity_Site\d?", other):
            return "DnaA at the origins", "DnaA binds a high-affinity site off oriC"
        if re.fullmatch(r"High_Affinity_Site_oriC\d?", other):
            return "DnaA at the origins", "DnaA binds the high-affinity site at oriC"
        if re.fullmatch(r"Low_Affinity_Site\d?_\d", other):
            return "DnaA at the origins", "DnaA binds a low-affinity site at oriC"
        if re.fullmatch(r"ssDNAunboundSite\d?_\d+", other):
            if "Initiator_C" in pos:
                return "DnaA at the origins", "DnaA binds the last ssDNA site (filament complete, makes Initiator_C and Initiator_CC)"
            return "DnaA at the origins", "DnaA extends the ssDNA filament"
    if "M_DnaA_c" in pos:
        if any(re.fullmatch(r"ssDNA_Unbinding\d?_\d+", n) for n in (a, b)):
            return "DnaA at the origins", "DnaA leaves the filament (unbinding cascade, k = 600)"
        if any(re.fullmatch(r"ssDNAunboundSite\d?_\d+", n) for n in (a, b)):
            return "DnaA at the origins", "DnaA leaves the filament tip (k = 0.55)"
    if {re.sub(r"\d+$", "", a), re.sub(r"\d+$", "", b)} in ({"Initiator_C", "chromosome_C"},
                                                           {"Initiator_CC", "chromosome_CC"}):
        return "Replication initiation", "Initiator + chromosome: the fork starts and copies its first gene"
    if re.fullmatch(r"[A-Z]{3}_cost", a) and re.fullmatch(r"M_[a-z]{3}trna_c|M_glntrna_c", b):
        return "Translation cost payment", "residue counter + charged tRNA -> uncharged tRNA"
    if "glutrnagln" in a + b or "glutrnagln" in " ".join(pos):
        return "tRNA charging", "Gln-tRNA(Gln) amidotransferase (GatCAB, JCVISYN3A_0689) steps"
    if "M_atp_c" in (a, b) and any(n.startswith("M_PTN_JCVISYN3A_") for n in (a, b)):
        return "tRNA charging", "synthetase binds ATP"
    if any(re.fullmatch(r"JCVISYN3A_\d+_ATP", n) for n in (a, b)) and any(re.fullmatch(AA_RE, n) for n in (a, b)):
        return "tRNA charging", "synthetase-ATP binds its amino acid"
    if any(re.fullmatch(r"JCVISYN3A_\d+_ATP_AA", n) for n in (a, b)) and any(re.fullmatch(r"M_trna[a-z]{3}_c", n) for n in (a, b)):
        return "tRNA charging", "synthetase-ATP-amino acid binds its tRNA"
    return None, None


def first_order_family(reac, prods):
    pos = {n for n, d in prods if d > 0}
    if "ATP_trsc" in pos:
        return "Transcription"
    if "ATP_translat" in pos:
        return "Translation"
    if "ATP_mRNAdeg" in pos:
        return "mRNA degradation"
    if "ATP_DNArep" in pos:
        return "DNA replication (fork step)"
    if "ATP_transloc" in pos:
        return "Translocation"
    (a,) = reac
    if re.fullmatch(r"JCVISYN3A_\d+_ATP_AA_tRNA(gln)?|glutrnagln_enz_atp_(gln|asn)", a):
        return "tRNA charging (release, first order)"
    return None


# ---------------------------------------------------------------------------
# The node DAG (README.md "nodes.csv"), evaluated with Python semantics.
# ---------------------------------------------------------------------------

class Dag:
    def __init__(self, path):
        self.op, self.args = {}, {}
        for r in read_csv(path):
            i = int(r["id"])
            self.op[i] = r["op"]
            self.args[i] = r["args"]
        self.n = len(self.op)
        assert sorted(self.op) == list(range(1, self.n + 1))
        self.kids = {i: ([] if self.op[i] in ("const", "count") else [int(x) for x in self.args[i].split(";")])
                     for i in self.op}
        # Ops and count leaves reachable from each node (ids are topological).
        self.ops_below, self.counts_below = {}, {}
        for i in range(1, self.n + 1):
            ops = {self.op[i]}
            cnt = {self.args[i]} if self.op[i] == "count" else set()
            for k in self.kids[i]:
                ops |= self.ops_below[k]
                cnt |= self.counts_below[k]
            self.ops_below[i] = frozenset(ops)
            self.counts_below[i] = frozenset(cnt)

    def reachable(self, roots):
        seen, stack = set(), list(roots)
        while stack:
            i = stack.pop()
            if i in seen:
                continue
            seen.add(i)
            stack.extend(self.kids[i])
        return seen

    def evaluate(self, counts):
        """Every node's value at `counts` (name -> int), README semantics."""
        v = {}
        for i in range(1, self.n + 1):
            op, k = self.op[i], self.kids[i]
            if op == "const":
                s = self.args[i]
                v[i] = float(s) if any(t in s for t in (".", "e", "inf", "nan")) else int(s)
            elif op == "count":
                v[i] = counts[self.args[i]]
            elif op == "add":
                v[i] = v[k[0]] + v[k[1]]
            elif op == "sub":
                v[i] = v[k[0]] - v[k[1]]
            elif op == "mul":
                v[i] = v[k[0]] * v[k[1]]
            elif op == "div":
                v[i] = v[k[0]] / v[k[1]]
            elif op == "pow":
                v[i] = v[k[0]] ** v[k[1]]
            elif op == "round":
                v[i] = round(v[k[0]])
            elif op == "min":
                v[i] = min(v[j] for j in k)
            elif op == "max":
                v[i] = max(v[j] for j in k)
            elif op == "select_gt":
                a, b, x, y = (v[j] for j in k)
                v[i] = x if a > b else y
            else:
                raise ValueError("op %s is not in this export's histogram" % op)
        return v

    def branch(self, i, v):
        """Which argument a min/max/select_gt node takes at values v, and whether it is a floor/ceiling."""
        op, k = self.op[i], self.kids[i]
        if op in ("min", "max"):
            vals = [v[j] for j in k]
            best = min(vals) if op == "min" else max(vals)
            return vals.index(best)
        if op == "select_gt":
            return 0 if v[k[0]] > v[k[1]] else 1
        return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", default=os.path.join(ENSEMBLE, "data", "model"))
    ap.add_argument("--speedycell", default=os.path.join(ENSEMBLE, "code", "SpeedyCell", "src"))
    ap.add_argument("--speedycell-jl", default=SPEEDYCELL_JL,
                    help="optional: a SpeedyCell.jl checkout, to compare its model tables")
    a = ap.parse_args()
    M, S = a.model, a.speedycell
    P = lambda *p: os.path.join(M, *p)

    # ---------------------------------------------------------------- provenance
    print("== Provenance")
    head = subprocess.run(["git", "-C", M, "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    print("syn3a-ensemble HEAD:", head or "(not a git checkout)")
    with open(P("export_meta.json")) as f:
        meta = json.load(f)
    print("export repo_git_sha:", meta["repo_git_sha"], " Minimal_Cell:", meta["minimal_cell_commit"])
    tables = ["cme/reactions.csv", "cme/rate_rules.csv", "cme/species.csv", "cme/initial_counts.csv", "nodes.csv",
              "ode/species.csv", "ode/reactions.csv", "ode/stoichiometry.csv", "export_meta.json", "inventory.md"]
    for t in tables:
        print("sha256 %s  %s" % (sha256(P(t)), t))
    for t in ["Cell/writeback.jl", "Cell/hook.jl"]:
        print("sha256 %s  SpeedyCell/src/%s" % (sha256(os.path.join(S, t)), t))
    J = a.speedycell_jl
    if os.path.isdir(os.path.join(J, "data", "model")):
        jhead = subprocess.run(["git", "-C", J, "rev-parse", "--short", "HEAD"], capture_output=True, text=True).stdout.strip()
        same = [t for t in ("cme/reactions.csv", "cme/rate_rules.csv", "cme/species.csv", "cme/initial_counts.csv")
                if sha256(os.path.join(J, "data", "model", t)) == sha256(P(t))]
        print("SpeedyCell.jl HEAD %s: CME tables byte-identical to the above: %d of 4 (%s); its inventory.md sha256 %s"
              % (jhead, len(same), ", ".join(same), sha256(os.path.join(J, "data", "model", "inventory.md"))))
    else:
        print("SpeedyCell.jl checkout not found at", J, "(comparison skipped)")

    sp_rows = read_csv(P("cme", "species.csv"))
    sp = {int(r["index"]): r["name"] for r in sp_rows}
    x0 = {r["name"]: int(r["initial_count"]) for r in sp_rows}
    rx = read_csv(P("cme", "reactions.csv"))
    rr = {int(r["index"]): int(r["rate_node"]) for r in read_csv(P("cme", "rate_rules.csv"))}
    ode_sp = read_csv(P("ode", "species.csv"))
    ode_ids = {r["id"] for r in ode_sp}
    ode_rx = read_csv(P("ode", "reactions.csv"))
    dag = Dag(P("nodes.csv"))
    assert len(sp) == meta["cme_species"] and len(rx) == meta["cme_reactions"]

    R = []  # parsed reactions
    for r in rx:
        i = int(r["index"])
        reac = [sp[int(x)] for x in r["reactants"].split(";")]
        prods = [(sp[j], d) for j, d in parse_stoich(r["stoichiometry"])]
        R.append(dict(i=i, t=int(r["lm_type"]), k=float(r["rate"]), reac=reac, prods=prods, node=rr[i]))

    # ---------------------------------------------------------------- order
    print("\n== Reaction order (README: lm_type 0-3 = propensities k, k*x, k*x*y, k*x(x-1)/2)")
    by_t = collections.Counter(x["t"] for x in R)
    print("lm_type counts:", dict(sorted(by_t.items())))
    for x in R:
        assert len(x["reac"]) == x["t"] or (x["t"] == 0 and len(x["reac"]) == 0), x
        if x["t"] == 2:
            assert x["reac"][0] != x["reac"][1]
    assert by_t[2] == 363, by_t
    print("second-order (type 2, two distinct reactants):", by_t[2])

    # ---------------------------------------------------------------- families
    for x in R:
        if x["t"] == 2:
            x["fam"], x["sub"] = second_order_family(x["reac"], x["prods"])
        else:
            x["fam"], x["sub"] = first_order_family(x["reac"], x["prods"]), None
        assert x["fam"] is not None, ("unclassified", x["i"], x["reac"], x["prods"][:4])

    print("\n== The 363 second-order reactions by family")
    fam2 = collections.Counter(x["fam"] for x in R if x["t"] == 2)
    sub2 = collections.Counter((x["fam"], x["sub"]) for x in R if x["t"] == 2)
    for f, n in fam2.most_common():
        print("%4d  %s" % (n, f))
        for (ff, s), m in sorted(sub2.items(), key=lambda z: -z[1]):
            if ff == f:
                print("        %4d  %s" % (m, s))
    assert sum(fam2.values()) == 363
    print("sum:", sum(fam2.values()))

    # The inputs of CellSA, and so of doubling time (t_SA2): writeback.jl's lipids and membrane proteins.
    with open(os.path.join(S, "Cell/writeback.jl")) as f:
        wbsrc = f.read()
    lipids = re.findall(r'\("(M_\w+_c)", [0-9.]+\)', re.search(r"const SA_DICT = \[(.*?)\]\n", wbsrc, re.S).group(1))
    memids = re.findall(r'"(\d{4})"', re.search(r"const MEM_PROT_IDS = \[(.*?)\]\n", wbsrc, re.S).group(1))
    cellsa_inputs = set(lipids) | {"M_PTN_JCVISYN3A_%s_c" % i for i in memids if i != "0779"} | {"ptsg", "ptsg_P"}
    assert cellsa_inputs <= set(x0), sorted(cellsa_inputs - set(x0))
    print("\nCellSA inputs (writeback.jl lipid_sa!): %d lipids, %d membrane-protein ids, %d species" % (len(lipids), len(memids), len(cellsa_inputs)))

    print("\n-- per family: rate node constant?, reactant initial-count range, and what the reactants/products touch")
    rate_rule_inputs = set()
    for x in R:
        rate_rule_inputs |= dag.counts_below[x["node"]]
    for f in fam2:
        xs = [x for x in R if x["t"] == 2 and x["fam"] == f]
        const = sum(dag.op[x["node"]] == "const" for x in xs)
        lo = [min(x0[n] for n in x["reac"]) for x in xs]
        hi = [max(x0[n] for n in x["reac"]) for x in xs]
        ks = sorted({x["k"] for x in xs})
        names = set()
        for x in xs:
            names |= set(x["reac"]) | {n for n, d in x["prods"] if d != 0}
        touch = dict(
            ode_species=sorted(n for n in names if n in ode_ids),
            proteins=sorted(n for n in names if n.startswith("M_PTN_")) + (["M_DnaA_c"] if "M_DnaA_c" in names else []),
            genes=sorted(n for n in names if n.endswith("_gene")),
            transcripts=sorted(n for n in names if n.startswith("M_RNA_")),
            read_by_rate_rules=sorted(n for n in names if n in rate_rule_inputs),
        )
        touch["cellsa_inputs"] = sorted(n for n in names if n in cellsa_inputs)
        print("* %s (%d)" % (f, len(xs)))
        print("    rate node is a const: %d of %d; first-minute k values: %s" % (const, len(xs), ks if len(ks) <= 8 else "%d distinct, %g to %g" % (len(ks), ks[0], ks[-1])))
        print("    smaller reactant's initial count: %d to %d; larger: %d to %d" % (min(lo), max(lo), min(hi), max(hi)))
        for key, v in touch.items():
            print("    %s: %d%s" % (key, len(v), (" " + ", ".join(v[:12]) + (" ..." if len(v) > 12 else "")) if v else ""))
        nonconst = [x for x in xs if dag.op[x["node"]] != "const"]
        for x in nonconst:
            print("    count-dependent rate: reaction %d reads %s" % (x["i"], sorted(dag.counts_below[x["node"]])))

    print("\n-- cost payment: which *_cost_paid counter each credits")
    paid = collections.Counter(n for x in R if x["fam"] == "Translation cost payment" for n, d in x["prods"] if d > 0 and n.endswith("_cost_paid"))
    print("   ", dict(sorted(paid.items())))
    writers = collections.Counter(n for x in R for n, d in x["prods"] if d != 0 and n.endswith("_cost_paid"))
    cost_paid_species = sorted(n for n in x0 if n.endswith("_cost_paid"))
    print("    *_cost_paid species: %d; written by any reaction: %s" % (len(cost_paid_species), dict(sorted(writers.items()))))
    all_leaves = set().union(*(dag.counts_below[i] for i in dag.op))
    print("    *_cost_paid read as a reactant: %d; read as a count by any node: %d" % (
        sum(n.endswith("_cost_paid") for x in R for n in x["reac"]), sum(n.endswith("_cost_paid") for n in all_leaves)))

    print("\n== All 2,365 reactions by family (inventory.md's table, refined)")
    famall = collections.Counter((x["fam"], x["t"]) for x in R)
    for (f, t), n in sorted(famall.items(), key=lambda z: (-z[1])):
        print("%5d  type %d  %s" % (n, t, f))
    print("total:", sum(famall.values()))
    inv_other = sum(n for (f, t), n in famall.items()
                    if f not in ("Transcription", "Translation", "mRNA degradation", "Translocation")
                    and not f.startswith("DNA replication") and f != "Replication initiation")
    print("inventory.md 'Other' (383) recomputed:", inv_other)
    rep_total = famall[("DNA replication (fork step)", 1)] + famall[("Replication initiation", 2)]
    print("inventory.md 'DNA replication' (493) recomputed: %d = %d first-order fork steps + %d second-order initiations"
          % (rep_total, famall[("DNA replication (fork step)", 1)], famall[("Replication initiation", 2)]))
    assert rep_total == 493 and inv_other == 383
    assert sum(famall.values()) == 2365

    # ---------------------------------------------------------------- rate rules
    print("\n== Rate rules (rate_rules.csv roots into nodes.csv)")
    dep = [x for x in R if dag.op[x["node"]] != "const"]
    print("count-dependent rate rules:", len(dep), " constant:", len(R) - len(dep))
    assert len(dep) == meta["cme_rates_depending_on_counts"] == 1434
    print("-- by family: count-dependent / total")
    for f in sorted({x["fam"] for x in R}):
        xs = [x for x in R if x["fam"] == f]
        print("  %4d / %4d  %s" % (sum(dag.op[x["node"]] != "const" for x in xs), len(xs), f))

    cme_nodes = dag.reachable([x["node"] for x in R])
    ode_roots = [int(m) for r in ode_rx for m in re.findall(r"N__(\d+)", r["rate"])]
    ode_init_roots = [int(r["init_node"]) for r in ode_sp]
    ode_rate_nodes = dag.reachable(ode_roots)
    ode_init_nodes = dag.reachable(ode_init_roots)
    print("-- nodes per op: whole DAG (export_meta) | reachable from CME rate roots | from ODE rate N__ refs | from ODE init_node")
    for op in NONSMOOTH:
        tot = sum(1 for i in dag.op if dag.op[i] == op)
        if tot == 0:
            continue
        assert tot == meta["node_ops"].get(op, 0)
        c = sum(1 for i in cme_nodes if dag.op[i] == op)
        o = sum(1 for i in ode_rate_nodes if dag.op[i] == op)
        oi = sum(1 for i in ode_init_nodes if dag.op[i] == op)
        print("  %-10s %4d | %4d | %4d | %4d" % (op, tot, c, o, oi))
    absent = [op for op in NONSMOOTH if meta["node_ops"].get(op, 0) == 0]
    print("  absent from the export:", ", ".join(absent))

    # Classify each min/max/round/select_gt node by what it is.
    v0 = dag.evaluate(x0)

    def node_kind(i):
        op, k = dag.op[i], dag.kids[i]
        if op == "max" and len(k) == 2 and dag.op[k[0]] == "const" and v0[k[0]] == 1 and dag.op[k[1]] == "count":
            return "max(1, count) floor on a pool before division"
        if op == "min" and len(k) == 2 and dag.op[k[1]] == "const" and v0[k[1]] == 180:
            return "min(c * protein/180, 180): promoter strength from a live protein count, capped at 180"
        if op == "min" and all(dag.op[j] == "count" or dag.op[j] == "min" for j in k):
            return "min over subunit counts: complex AND rule (ODE enzyme)"
        if op == "round":
            return "round(V * 1e19): CellV"
        if op == "select_gt":
            return "select_gt(V, 6.7e-17, ...): the volume cap"
        raise ValueError("node %d (%s) fits no known kind" % (i, op))

    kinds = collections.defaultdict(list)
    for i in dag.op:
        if dag.op[i] in ("min", "max", "round", "select_gt"):
            kinds[node_kind(i)].append(i)
    print("-- the nonsmooth nodes by kind (node count; where reachable; active branch at procid 1's initial counts)")
    for kd, ids in sorted(kinds.items(), key=lambda z: -len(z[1])):
        where = []
        if any(i in cme_nodes for i in ids):
            where.append("CME rules %d" % sum(i in cme_nodes for i in ids))
        if any(i in ode_rate_nodes for i in ids):
            where.append("ODE rates %d" % sum(i in ode_rate_nodes for i in ids))
        if any(i in ode_init_nodes for i in ids):
            where.append("ODE init %d" % sum(i in ode_init_nodes for i in ids))
        if kd.startswith("max(1"):
            act = "floor active (pool < 1) at t0: %d" % sum(dag.branch(i, v0) == 0 for i in ids)
        elif kd.startswith("min(c"):
            act = "ceiling binding at t0: %d" % sum(dag.branch(i, v0) == 1 for i in ids)
        elif kd.startswith("select_gt"):
            act = "at the cap at t0: %d" % sum(dag.branch(i, v0) == 0 for i in ids)
        elif kd.startswith("min over"):
            act = "arguments per node: %s" % dict(collections.Counter(len(dag.kids[i]) for i in ids))
        else:
            act = "CellV at t0 = %s" % v0[ids[0]]
        print("  %4d  %s  [%s]  %s" % (len(ids), kd, "; ".join(where) or "unreached", act))

    floors = sorted(dag.args[dag.kids[i][1]] for i in kinds["max(1, count) floor on a pool before division"])
    print("  the pools under a max(1, .) floor (%d): %s" % (len(floors), " ".join(floors)))
    print("-- count-dependent CME rules reaching each nonsmooth op, by family")
    rows = collections.defaultdict(collections.Counter)
    for x in dep:
        for op in ("min", "max", "round", "select_gt"):
            if op in dag.ops_below[x["node"]]:
                rows[x["fam"]][op] += 1
        if not ({"min", "max", "round", "select_gt"} & dag.ops_below[x["node"]]):
            rows[x["fam"]]["smooth only"] += 1
    any_ns = sum(1 for x in dep if {"min", "max", "round", "select_gt"} & dag.ops_below[x["node"]])
    print("   family: min max round select_gt smooth-only")
    for f, c in sorted(rows.items()):
        print("   %-40s %4d %4d %4d %4d %4d" % (f, c["min"], c["max"], c["round"], c["select_gt"], c["smooth only"]))
    print("  count-dependent rules reaching at least one of min/max/round/select_gt: %d of %d" % (any_ns, len(dep)))
    with_min = [x for x in dep if "min" in dag.ops_below[x["node"]]]
    print("  rules reaching min (the protein-count promoter ceiling): %d, all %s" % (len(with_min), sorted({x["fam"] for x in with_min})))
    own = 0
    for x in with_min:
        prot = sorted(n for n in dag.counts_below[x["node"]] if n.startswith("M_PTN_"))
        gene = x["reac"][0]
        tag = re.search(r"JCVISYN3A_(\d+)", gene).group(1)
        own += prot == ["M_PTN_JCVISYN3A_%s_c" % tag]
    print("  ... of which the protein read is the gene's own product: %d of %d" % (own, len(with_min)))
    print("  those genes:", " ".join(sorted(re.search(r"JCVISYN3A_(\d+)", x["reac"][0]).group(1) for x in with_min)))
    reads = collections.Counter()
    for x in dep:
        for n in dag.counts_below[x["node"]]:
            reads[n] += 1
    print("  distinct count species read by the CME rules: %d" % len(reads))
    trl = [x for x in R if x["fam"] == "Translation"]
    per = collections.Counter(sum("trna" in n for n in dag.counts_below[x["node"]]) for x in trl)
    print("  charged tRNA pools read per translation rule: %s" % dict(per))
    groups = collections.Counter()
    for n in reads:
        if n.startswith("M_PTN_") or n == "ptsg":
            groups["protein counts"] += 1
        elif n in ("CellSA", "CellV"):
            groups[n] += 1
        elif re.fullmatch(r"M_d[atgc]tp_c", n):
            groups["dNTP pools"] += 1
        elif re.fullmatch(r"M_[aucg]tp_c", n):
            groups["NTP pools"] += 1
        elif "trna" in n:
            groups["charged tRNA pools"] += 1
        else:
            groups["other: " + n] += 1
    print("   ", dict(sorted(groups.items())))
    print("  count species read by the ODE (rates + init): %d" % len(set().union(*(dag.counts_below[i] for i in ode_roots + ode_init_roots))))
    ode_with = collections.Counter()
    for r in ode_rx:
        ns = [int(m) for m in re.findall(r"N__(\d+)", r["rate"])]
        ops = set().union(*(dag.ops_below[i] for i in ns)) if ns else set()
        for op in ("min", "max", "round", "select_gt"):
            ode_with[op] += op in ops
        ode_with["reads a count node"] += bool(ns)
    print("  ODE reactions (of %d) whose rate reaches: %s" % (len(ode_rx), dict(sorted(ode_with.items()))))

    # ---------------------------------------------------------------- code rules
    print("\n== Count-dependent rules that are code, not data (SpeedyCell/src/Cell)")
    for fn in ("Cell/writeback.jl", "Cell/hook.jl"):
        with open(os.path.join(S, fn)) as f:
            lines = f.read().split("\n")
        for ln, line in enumerate(lines, 1):
            code = line.split("#")[0]
            if re.search(r"\b(min|max|round|trunc|floor)\(", code) or re.search(r"V > V_CAP", code):
                print("  %s:%d: %s" % (fn, ln, line.strip()))
    with open(os.path.join(S, "Cell/writeback.jl")) as f:
        wb = f.read()

    def count_list(name, width):
        """Entries of `const name = [...]`: bare strings (width 1) or tuples of `width` strings."""
        body = re.search(r"const %s = \[(.*?)\]\n" % name, wb, re.S).group(1)
        strings = re.findall(r'"\w+"', body)
        if width == 1:
            assert "(" not in body, name
            return len(strings)
        tuples = re.findall(r"\(([^()]*)\)", body)
        assert all(len(re.findall(r'"\w+"', t)) == width for t in tuples), name
        assert len(strings) == width * len(tuples), name
        return len(tuples)

    capped = {n: count_list(n, w) for n, w in (("ATP_HYDRO_COUNTERS", 1), ("NTP_COUNTERS", 2), ("DNTP_COUNTERS", 2), ("AATRNA_COUNTERS", 3))}
    assert sum(capped.values()) == 22, capped
    print("  cost counters charged with min(counter, pool):", capped, "total", sum(capped.values()))
    print("  NMP recycling counters (credited whole, no cap):", count_list("NMP_COUNTERS", 2))
    print("  ODE species rounded to counts at every 1 s hook:", len(ode_sp))

    # ---------------------------------------------------------------- genes
    print("\n== Gene copy number")
    genes = sorted(n for n in x0 if n.endswith("_gene"))
    print("*_gene species: %d; initial counts: %s" % (len(genes), dict(collections.Counter(x0[g] for g in genes))))
    tx = [x for x in R if x["fam"] == "Transcription"]
    assert all(x["t"] == 1 and x["reac"][0].endswith("_gene") for x in tx)
    assert all(all(n != x["reac"][0] for n, d in x["prods"]) for x in tx)
    txg = collections.Counter(x["reac"][0] for x in tx)
    print("transcription reactions: %d, all type 1 with one *_gene reactant that the reaction leaves unchanged (propensity k*G)" % len(tx))
    print("genes with a transcription reaction: %d (each has %s); genes without: %d %s"
          % (len(txg), sorted(set(txg.values())), len(genes) - len(txg), sorted(set(genes) - set(txg))))
    prods_of_tx = collections.Counter()
    for x in tx:
        made = sorted(n for n, d in x["prods"] if d > 0 and not re.fullmatch(r"[ACGUT]TP_\w+", n))
        kind = ("mRNA (M_RNA_*)" if all(n.startswith("M_RNA_") for n in made) else
                "tRNA" if all("trna" in n for n in made) else
                "rRNA operon (%s)" % " ".join(made) if all("rRNA" in n for n in made) else " ".join(made))
        prods_of_tx[kind] += 1
    print("transcription reactions by product:", dict(sorted(prods_of_tx.items())))
    gene_readers = collections.Counter(x["fam"] for x in R for g in x["reac"] if g.endswith("_gene"))
    print("reactions with a *_gene reactant, by family:", dict(sorted(gene_readers.items())))
    print("rate rules that read a *_gene count:", sum(any(n.endswith("_gene") for n in dag.counts_below[x["node"]]) for x in R))
    tx_reads = collections.Counter()
    for x in tx:
        for n in dag.counts_below[x["node"]]:
            tx_reads["protein" if n.startswith("M_PTN_") else n] += 1
    print("what the transcription rate rules read (rules per input):", dict(sorted(tx_reads.items())))
    print("transcription rules reaching min / round / select_gt / max: %d / %d / %d / %d" % tuple(
        sum(op in dag.ops_below[x["node"]] for x in tx) for op in ("min", "round", "select_gt", "max")))

    rep = [x for x in R if x["fam"] in ("DNA replication (fork step)", "Replication initiation")]
    inc = collections.Counter(n for x in rep for n, d in x["prods"] if n.endswith("_gene") and d > 0)
    dec = [n for x in R for n, d in x["prods"] if n.endswith("_gene") and d < 0]
    print("replication reactions: %d; each increments one gene by +1: %s; genes incremented: %d; reactions that decrement a gene: %d"
          % (len(rep), all(sum(1 for n, d in x["prods"] if n.endswith("_gene")) == 1 for x in rep), len(inc), len(dec)))
    assert len(rep) == 493 and len(inc) == len(genes) and set(inc.values()) == {1}
    # The two forks: follow *_inter species from each initiation.
    nxt = {x["reac"][0]: x for x in rep if x["t"] == 1}
    for start in [x for x in rep if x["t"] == 2]:
        cur = [n for n, d in start["prods"] if n.endswith("_inter") and d > 0][0]
        steps, last = 0, start
        while cur in nxt:
            last = nxt[cur]
            steps += 1
            outs = [n for n, d in last["prods"] if n.endswith("_inter") and d > 0]
            if not outs:
                break
            cur = outs[0]
        first_gene = [n for n, d in start["prods"] if n.endswith("_gene")][0]
        print("fork from reaction %d (%s + %s, k = %g): copies %s first, then %d first-order steps, ending at reaction %d"
              % (start["i"], start["reac"][0], start["reac"][1], start["k"], first_gene, steps, last["i"]))
        # Where M0's two genes sit on this fork (step 0 is the initiation).
        cur, pos = start, 0
        while True:
            g = [n for n, d in cur["prods"] if n.endswith("_gene")][0]
            if g in ("JCVISYN3A_0607_gene", "JCVISYN3A_0779_gene"):
                print("    M0 gene %s is copied at step %d of this fork (reaction %d, first-minute k = %g)" % (g, pos, cur["i"], cur["k"]))
            outs = [n for n, d in cur["prods"] if n.endswith("_inter") and d > 0]
            if not outs or outs[0] not in nxt:
                break
            cur, pos = nxt[outs[0]], pos + 1
    for x in rep:
        extra = sorted("%s(%+d)" % (n, d) for n, d in x["prods"] if not re.search(r"_(gene|inter)$|DNArep", n))
        if extra:
            print("replication reaction %d (%s, copies %s) also makes: %s" % (
                x["i"], " + ".join(x["reac"]), [n for n, d in x["prods"] if n.endswith("_gene")][0], " ".join(extra)))
    rep_reads = set().union(*(dag.counts_below[x["node"]] for x in rep))
    print("what the replication rate rules read:", sorted(rep_reads))
    print("replication rules count-dependent: %d of %d; reaching max / round / select_gt: %d / %d / %d" % (
        sum(dag.op[x["node"]] != "const" for x in rep), len(rep),
        *(sum(op in dag.ops_below[x["node"]] for x in rep) for op in ("max", "round", "select_gt"))))
    dnaa_makers = [x["i"] for x in R if x["fam"] == "Translation" and any(n == "M_DnaA_c" and d > 0 for n, d in x["prods"])]
    print("M_DnaA_c is made by translation reaction(s): %s (gene %s)" % (dnaa_makers, [x["reac"][0] for x in R if x["i"] in dnaa_makers]))
    print("initial counts: M_DnaA_c %d, High_Affinity_Site_oriC %d, chromosome_C %d, chromosome_CC %d"
          % (x0["M_DnaA_c"], x0["High_Affinity_Site_oriC"], x0["chromosome_C"], x0["chromosome_CC"]))

    # ---------------------------------------------------------------- checks the review asked for
    print("\n== Rate constants: the first minute against the rule from minute 2 on")
    for x in R:
        if dag.op[x["node"]] == "const" and v0[x["node"]] != x["k"]:
            print("  reaction %d (%s, %s): first-minute k = %g, rule from minute 2 = const %g"
                  % (x["i"], x["fam"], " + ".join(x["reac"]), x["k"], v0[x["node"]]))
    made_by = collections.defaultdict(list)
    for x in R:
        for n, d in x["prods"]:
            if d > 0:
                made_by[n].append(x["i"])
    dead = [x for x in R if any(x0[n] == 0 and not made_by[n] for n in x["reac"])]
    print("reactions that can never fire (a reactant starts at 0 and nothing makes it): %s"
          % ["%d (%s)" % (x["i"], " + ".join(x["reac"])) for x in dead])
    print("reactions that make Initiator_C: %s; Initiator_CC: %s" % (made_by["Initiator_C"], made_by["Initiator_CC"]))
    rrna = [x for x in tx if any("rRNA" in n for n, d in x["prods"] if d > 0)]
    print("rRNA operon reactions: %s" % ["%d (%s)" % (x["i"], x["reac"][0]) for x in rrna])
    if os.path.exists(REDUCED_GENES):
        with open(REDUCED_GENES) as f:
            red = [l.split("\t")[0] for l in f if l.startswith("JCVISYN3A_")]
        own29 = {re.search(r"JCVISYN3A_\d+", x["reac"][0]).group(0) for x in with_min}
        print("reduced model genes (%s): %d; shared with the 29 protein-count promoters: %d"
              % (os.path.relpath(REDUCED_GENES, HERE), len(red), len(own29 & set(red))))

    # ---------------------------------------------------------------- exact closure
    print("\n== Exact closure of the DnaA cascade, by its reachable states")
    # One origin's cascade, from its oriC site at 1 and every other cascade species at 0.
    # DnaA is treated as always available, so every binding the sites allow can fire.
    # The off-oriC high-affinity sites (multi-copy) are not part of an origin's state.
    dnaa = [x for x in R if x["fam"] == "DnaA at the origins"]
    offsite = [x for x in dnaa if any(re.fullmatch(r"High_Affinity_Site\d?", n) for n in x["reac"])]
    cascade = [x for x in dnaa if x not in offsite]
    outside = {"M_DnaA_c", "Initiator_C", "Initiator_CC"}

    def reachable_states(rxns, start):
        seen, todo = {start}, [start]
        while todo:
            st = todo.pop()
            c = dict(st)
            for x in rxns:
                if all(n in outside or c.get(n, 0) >= 1 for n in x["reac"]):
                    nc = dict(c)
                    for n, d in x["prods"]:
                        if n not in outside:
                            nc[n] = nc.get(n, 0) + d
                    assert all(v >= 0 for v in nc.values())
                    ns = tuple(sorted((n, v) for n, v in nc.items() if v))
                    if ns not in seen:
                        seen.add(ns)
                        todo.append(ns)
        return [dict(st) for st in seen]

    closes, bilinear, fired = [], [], set()
    for oric in ("High_Affinity_Site_oriC", "High_Affinity_Site_oriC2", "High_Affinity_Site_oriC3"):
        states = reachable_states(cascade, ((oric, 1),))
        mine = [x for x in cascade if any(any(st.get(n, 0) for st in states) for n in x["reac"] if n not in outside)
                and all(n in outside or any(st.get(n, 0) for st in states) for n in x["reac"])]
        fired |= {x["i"] for x in mine}
        sites = [x for x in mine if "M_DnaA_c" not in x["reac"]]
        binds = [x for x in mine if "M_DnaA_c" in x["reac"]]
        for x in sites:
            a_, b_ = x["reac"]
            prod = [st.get(a_, 0) * st.get(b_, 0) for st in states]
            same_a = prod == [st.get(a_, 0) for st in states]
            same_b = prod == [st.get(b_, 0) for st in states]
            assert same_a or same_b, ("site x site product is not one of its factors", x["i"])
            closes.append(x["i"])
        for x in binds:
            site = [n for n in x["reac"] if n != "M_DnaA_c"][0]
            assert {st.get(site, 0) for st in states} <= {0, 1}, (x["i"], site)
            bilinear.append(x["i"])
        print("  origin %s: %d reachable states; %d site x site reactions, each equal to one of its factors in every state; "
              "%d DnaA x site bindings, each site 0 or 1 in every state" % (oric, len(states), len(sites), len(binds)))
    assert len(closes) == 177 and len(bilinear) == 99 and fired == {x["i"] for x in cascade}
    offlive = [x for x in offsite if x not in dead]
    print("  off-oriC high-affinity bindings (DnaA x a multi-copy site): %d, of which can fire: %s"
          % (len(offsite), ["%d (%s, initial %d)" % (x["i"], x["reac"][1], x0[x["reac"][1]]) for x in offlive]))
    print("  so of the 279: %d close exactly (site x site), %d are DnaA x a 0/1 site, %d DnaA x a multi-copy site, %d dead"
          % (len(closes), len(bilinear), len(offlive), len(offsite) - len(offlive)))
    assert len(closes) + len(bilinear) + len(offsite) == 279

    def affine_on(states, keys, f):
        """True if f(state) is an affine function of the state over `states` (exact rational elimination)."""
        rows = [[fractions.Fraction(st.get(k, 0)) for k in keys] + [fractions.Fraction(1), fractions.Fraction(f(st))] for st in states]
        ncol = len(keys) + 1
        r = 0
        for c in range(ncol):
            piv = next((i for i in range(r, len(rows)) if rows[i][c] != 0), None)
            if piv is None:
                continue
            rows[r], rows[piv] = rows[piv], rows[r]
            for i in range(len(rows)):
                if i != r and rows[i][c] != 0:
                    q = rows[i][c] / rows[r][c]
                    rows[i] = [u - q * w for u, w in zip(rows[i], rows[r])]
            r += 1
        return all(any(v != 0 for v in row[:ncol]) or row[ncol] == 0 for row in rows)

    def toy(rxns, start):
        """Reachable states of a small hand-written model: rxns = [(reactants, {species: change})]."""
        seen, todo = {start}, [start]
        while todo:
            st = dict(todo.pop())
            for reac, ch in rxns:
                if all(st.get(n, 0) >= 1 for n in reac):
                    nc = dict(st)
                    for n, d in ch.items():
                        nc[n] = nc.get(n, 0) + d
                    ns = tuple(sorted(nc.items()))
                    if ns not in seen:
                        seen.add(ns)
                        todo.append(ns)
        return [dict(s_) for s_ in seen]

    print("-- initiation, Initiator x chromosome (hand-written abstractions of the tables)")
    # First round: one origin makes one Initiator (reaction 1656); one chromosome copy (initial 1).
    first = toy([((), {}), (("O",), {"O": -1, "I": 1}), (("I", "C"), {"I": -1, "C": -1})], (("C", 1), ("I", 0), ("O", 1)))
    a1 = affine_on(first, ["O", "I", "C"], lambda st: st["I"] * st["C"])
    print("  first round (O=1, I=0, C=1): %d states %s; I*C affine in the state: %s"
          % (len(first), sorted(tuple(st[k] for k in "OIC") for st in first), a1))
    # After reaction 2364: chromosome +2, two daughter origins, each making one Initiator.
    second = toy([(("Oa",), {"Oa": -1, "I": 1}), (("Ob",), {"Ob": -1, "I": 1}), (("I", "C"), {"I": -1, "C": -1})],
                 (("C", 2), ("I", 0), ("Oa", 1), ("Ob", 1)))
    a2 = affine_on(second, ["Oa", "Ob", "I", "C"], lambda st: st["I"] * st["C"])
    print("  second round (Oa=Ob=1, I=0, C=2): %d states; I*C affine in the state: %s" % (len(second), a2))
    # The M13 model's R2 without and with a DnaA source (D capped at 3 extra for the enumeration).
    noR1 = toy([(("D", "O"), {"D": -1, "O": -1, "F": 1})], (("D", 3), ("F", 0), ("O", 1)))
    print("  M13 R2 without R1 (D0=3): %d states; D*O affine in the state: %s"
          % (len(noR1), affine_on(noR1, ["D", "O", "F"], lambda st: st["D"] * st["O"])))
    withR1 = [st for st in toy([(("S",), {"S": -1, "D": 1}), (("D", "O"), {"D": -1, "O": -1, "F": 1})],
                                (("D", 3), ("F", 0), ("O", 1), ("S", 3)))]
    print("  M13 R2 with R1 (D0=3, three births): %d states; D*O affine in the state: %s"
          % (len(withR1), affine_on(withR1, ["D", "O", "F", "S"], lambda st: st["D"] * st["O"])))
    assert a1 and not a2

    # ---------------------------------------------------------------- hops to the observables
    print("\n== Hops from each kind to the nearest observable, in the dependency graph")
    # Edge a -> b, labelled by what links them: a reaction or code rule that reads a and changes b.
    edges = collections.defaultdict(dict)

    def link(a_, b_, lab):
        if a_ != b_ and b_ not in edges[a_]:
            edges[a_][b_] = lab

    for x in R:
        reads = set(x["reac"]) | dag.counts_below[x["node"]]
        for n, d in x["prods"]:
            if d != 0:
                for r_ in reads:
                    link(r_, n, "CME %d" % x["i"])
    ode_st = collections.defaultdict(dict)
    for r in read_csv(P("ode", "stoichiometry.csv")):
        ode_st[int(r["reaction"])][r["species"]] = float(r["coefficient"])
    for r in ode_rx:
        num = int(r["number"])
        reads = {t for t in re.findall(r"[A-Za-z_][A-Za-z0-9_]*", r["rate"]) if t in ode_ids}
        for m in re.findall(r"N__(\d+)", r["rate"]):
            reads |= dag.counts_below[int(m)]
        for n, c in ode_st[num].items():
            if c != 0:
                for r_ in reads:
                    link(r_, n, "ODE %s" % r["id"])
    # The hook's code rules (writeback.jl): counters paid from pools, the surface area, the volume.
    for nm, w in (("NTP_COUNTERS", 2), ("DNTP_COUNTERS", 2), ("NMP_COUNTERS", 2)):
        for c, m in re.findall(r'\("(\w+)", "(\w+)"\)', re.search(r"const %s = \[(.*?)\]\n" % nm, wb, re.S).group(1)):
            link(c, m, "writeback.jl counter %s" % c)
    for c in re.findall(r'"(\w+)"', re.search(r"const ATP_HYDRO_COUNTERS = \[(.*?)\]\n", wb, re.S).group(1)):
        for m in (("M_gtp_c", "M_gdp_c") if "translat" in c else ("M_atp_c", "M_adp_c")) + ("M_pi_c",):
            link(c, m, "writeback.jl counter %s" % c)
    for c, ch, un in re.findall(r'\("(\w+)", "(\w+)", "(\w+)"\)', re.search(r"const AATRNA_COUNTERS = \[(.*?)\]\n", wb, re.S).group(1)):
        link(c, ch, "writeback.jl counter %s" % c)
        link(c, un, "writeback.jl counter %s" % c)
    for n in cellsa_inputs:
        link(n, "CellSA", "writeback.jl lipid_sa!")
    link("CellSA", "CellV", "writeback.jl cell_volume")
    obs = {"transcripts": {n for n in x0 if n.startswith("M_RNA_")}, "bulk metabolites": set(ode_ids),
           "doubling time (CellSA)": {"CellSA"}}

    def dist_to(target):
        """Hops from every species to the target set, and the next species on a shortest route."""
        inv = collections.defaultdict(list)
        for a_, bs in edges.items():
            for b_ in bs:
                inv[b_].append(a_)
        dist, nxt_ = {t: 0 for t in target}, {}
        frontier = sorted(target)
        while frontier:
            new = []
            for b_ in frontier:
                for a_ in sorted(inv[b_]):
                    if a_ not in dist:
                        dist[a_] = dist[b_] + 1
                        nxt_[a_] = b_
                        new.append(a_)
            frontier = new
        return dist, nxt_

    dists = {k: dist_to(v) for k, v in obs.items()}

    def kind_hops(writes, label, offset=1):
        """offset + the hops from the species the kind changes (offset 1: the kind's own reaction is a hop;
        0 for a kind that is a species); the route to the nearest observable of each type."""
        out = []
        for ob, (dist, nxt_) in dists.items():
            best = min(((dist[n], n) for n in writes if n in dist), default=None)
            if best is None:
                out.append("%s: unreached" % ob)
                continue
            route, n = [best[1]], best[1]
            while n in nxt_:
                lab = edges[n][nxt_[n]]
                n = nxt_[n]
                route.append("[%s] %s" % (lab, n))
            out.append("%s: %d (%s)" % (ob, best[0] + offset, " -> ".join(route)))
        print("* %s\n    %s" % (label, "\n    ".join(out)))

    def writes_of(xs):
        return sorted({n for x in xs for n, d in x["prods"] if d != 0})

    kind_hops(writes_of([x for x in R if x["i"] in bilinear]), "DnaA x 0/1 site bindings (99)")
    kind_hops(writes_of(offlive), "DnaA x multi-copy site bindings (%d live)" % len(offlive))
    kind_hops(writes_of([x for x in R if x["i"] in closes]), "cascade site x site (177)")
    kind_hops(writes_of([x for x in R if x["fam"] == "Replication initiation"]), "replication initiations (2)")
    kind_hops(writes_of([x for x in R if x["fam"] == "tRNA charging"]), "tRNA charging, second order (62)")
    kind_hops(writes_of([x for x in R if x["fam"] == "Translation cost payment"]), "translation cost payment (20)")
    kind_hops(writes_of([x for x in R if x["fam"] == "DNA replication (fork step)"]), "fork steps (491)")
    kind_hops(genes, "a gene copy number (the *_gene species themselves)", offset=0)
    kind_hops(writes_of(with_min), "transcription with an own-protein promoter (29)")
    kind_hops(writes_of(dep), "every count-dependent CME rule (1,434), so every floor, CellV and cap reader")
    and_rx = [r for r in ode_rx if any("min" in dag.ops_below[int(m)] for m in re.findall(r"N__(\d+)", r["rate"]))]
    kind_hops(sorted({n for r in and_rx for n, c in ode_st[int(r["number"])].items() if c != 0}), "ODE reactions with an AND rule (%d)" % len(and_rx))

    print("\nOK: 363 second-order, 493 replication, 1434 count-dependent rules; every reaction classified.")


if __name__ == "__main__":
    sys.exit(main())

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
import hashlib
import os
import re
import subprocess
import sys

ENSEMBLE = "/fred/oz022/tkimpson/syn3a-ensemble"
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
    a = ap.parse_args()
    M, S = a.model, a.speedycell
    P = lambda *p: os.path.join(M, *p)

    # ---------------------------------------------------------------- provenance
    print("== Provenance")
    head = subprocess.run(["git", "-C", M, "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    print("syn3a-ensemble HEAD:", head or "(not a git checkout)")
    meta = {}
    with open(P("export_meta.json")) as f:
        import json
        meta = json.load(f)
    print("export repo_git_sha:", meta["repo_git_sha"], " Minimal_Cell:", meta["minimal_cell_commit"])
    tables = ["cme/reactions.csv", "cme/rate_rules.csv", "cme/species.csv", "nodes.csv",
              "ode/species.csv", "ode/reactions.csv", "export_meta.json", "inventory.md"]
    for t in tables:
        print("sha256 %s  %s" % (sha256(P(t)), t))
    for t in ["Cell/writeback.jl", "Cell/hook.jl"]:
        print("sha256 %s  SpeedyCell/src/%s" % (sha256(os.path.join(S, t)), t))

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
    print("   ", dict(paid))
    writers = collections.Counter(n for x in R for n, d in x["prods"] if d != 0 and n.endswith("_cost_paid"))
    cost_paid_species = sorted(n for n in x0 if n.endswith("_cost_paid"))
    print("    *_cost_paid species: %d; written by any reaction: %s" % (len(cost_paid_species), dict(writers)))
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
        if op == "min":
            return "min over subunit/mixed terms: complex AND rule (ODE enzyme)"
        if op == "round":
            return "round(V * 1e19): CellV"
        if op == "select_gt":
            return "select_gt(V, 6.7e-17, ...): the volume cap"
        return op

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
    print("   ", dict(groups))
    print("  count species read by the ODE (rates + init): %d" % len(set().union(*(dag.counts_below[i] for i in ode_roots + ode_init_roots))))
    ode_with = collections.Counter()
    for r in ode_rx:
        ns = [int(m) for m in re.findall(r"N__(\d+)", r["rate"])]
        ops = set().union(*(dag.ops_below[i] for i in ns)) if ns else set()
        for op in ("min", "max", "round", "select_gt"):
            ode_with[op] += op in ops
        ode_with["reads a count node"] += bool(ns)
    print("  ODE reactions (of %d) whose rate reaches: %s" % (len(ode_rx), dict(ode_with)))

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

    def count_list(name):
        m = re.search(r"const %s = \[(.*?)\]\n" % name, wb, re.S)
        return len(re.findall(r'"\w+"', m.group(1))) // (1 if name == "ATP_HYDRO_COUNTERS" else (3 if name == "AATRNA_COUNTERS" else 2))

    capped = {n: count_list(n) for n in ("ATP_HYDRO_COUNTERS", "NTP_COUNTERS", "DNTP_COUNTERS", "AATRNA_COUNTERS")}
    print("  cost counters charged with min(counter, pool):", capped, "total", sum(capped.values()))
    print("  NMP recycling counters (credited whole, no cap):", count_list("NMP_COUNTERS"))
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
        for n, d in x["prods"]:
            if d > 0 and not re.fullmatch(r"[ACGUT]TP_\w+", n):
                prods_of_tx["mRNA (M_RNA_*)" if n.startswith("M_RNA_") else ("tRNA" if "trna" in n else ("rRNA" if "rRNA" in n else n))] += 1
    print("transcription products:", dict(prods_of_tx))
    gene_readers = collections.Counter(x["fam"] for x in R for g in x["reac"] if g.endswith("_gene"))
    print("reactions with a *_gene reactant, by family:", dict(gene_readers))
    print("rate rules that read a *_gene count:", sum(any(n.endswith("_gene") for n in dag.counts_below[x["node"]]) for x in R))
    tx_reads = collections.Counter()
    for x in tx:
        for n in dag.counts_below[x["node"]]:
            tx_reads["protein" if n.startswith("M_PTN_") else n] += 1
    print("what the transcription rate rules read (rules per input):", dict(tx_reads))
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

    print("\nOK: 363 second-order, 493 replication, 1434 count-dependent rules; every reaction classified.")


if __name__ == "__main__":
    sys.exit(main())

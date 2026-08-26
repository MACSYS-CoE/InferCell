"""Extract hand-coded model.addReaction blocks from defMetRxns.py (db048ac).

Sequential scan: attribute addSubstrate/addProduct/addParameter to the most
recent addReaction unless an explicit quoted reaction ID is given.
"""
import re

SRC = "mc/CME_ODE/program/defMetRxns.py"

RX_RXN   = re.compile(r"model\.addReaction\(\s*([^,]+?)\s*,")
RX_SUB   = re.compile(r"model\.addSubstrate\(\s*([^,]+?)\s*,\s*'(Sub\d)'\s*,\s*'([^']+)'")
RX_PROD  = re.compile(r"model\.addProduct\(\s*([^,]+?)\s*,\s*'(Prod\d)'\s*,\s*'([^']+)'")
RX_CLAMP = re.compile(r"model\.addParameter\(\s*([^,]+?)\s*,\s*'(Sub\d|Prod\d)'\s*,\s*([A-Za-z_]\w*|[\d.]+)\s*\)")

AA_TRANSPORT = [  # defMetRxns.py:2455 -- (symport, abc, cytoplasmic metabolite)
    ('R_ARGt2r','R_ARG4abc','M_arg__L_c'), ('R_ASPt2pr','R_ASP4abc','M_asp__L_c'),
    ('R_GLYt2r','R_GLY4abc','M_gly_c'),    ('R_ISOt2r','R_ISO4abc','M_ile__L_c'),
    ('R_ALAt2r','R_ALA4abc','M_ala__L_c'), ('R_ASNt2r','R_ASN4abc','M_asn__L_c'),
    ('R_LEUt2r','R_LEU4abc','M_leu__L_c'), ('R_HISt2r','R_HIS4abc','M_his__L_c'),
    ('R_LYSt2r','R_LYS4abc','M_lys__L_c'), ('R_PROt2r','R_PRO4abc','M_pro__L_c'),
    ('R_PHEt2r','R_PHE4abc','M_phe__L_c'), ('R_THRt2r','R_THR4abc','M_thr__L_c'),
    ('R_TRPt2r','R_TRP4abc','M_trp__L_c'), ('R_TYRt2r','R_TYR4abc','M_tyr__L_c'),
    ('R_VALt2r','R_VAL4abc','M_val__L_c'), ('R_SERt2r','R_SER4abc','M_ser__L_c'),
    ('R_METt2r','R_MET4abc','M_met__L_c'), ('R_CYSt2r','R_CYS4abc','M_cys__L_c'),
    ('R_GLUt2pr','R_GLU4abc','M_glu__L_c'),('R_GLNt2r','R_GLN4abc','M_gln__L_c'),
]

def _unq(tok):
    return tok.strip().strip("'\"") if tok.strip().startswith(("'", '"')) else None

def handcoded_reactions():
    lines = open(SRC).read().split("\n")
    rxns, order, cur = {}, [], None
    for i, l in enumerate(lines):
        if l.lstrip().startswith("#"):
            continue
        m = RX_RXN.search(l)
        if m:
            rid = _unq(m.group(1))
            if rid is None:                       # model.addReaction(rxnID, ...) in a loop
                cur = None
                continue
            cur = rid
            if rid not in rxns:
                rxns[rid] = dict(id=rid, line=i+1, subs=[], prods=[], clamped=[])
                order.append(rid)
            continue
        for rx, key in ((RX_SUB, "subs"), (RX_PROD, "prods")):
            m = rx.search(l)
            if m:
                rid = _unq(m.group(1)) or cur
                if rid in rxns:
                    rxns[rid][key].append(m.group(3))
        m = RX_CLAMP.search(l)
        if m:
            rid = _unq(m.group(1)) or cur
            if rid in rxns:
                rxns[rid]["clamped"].append((m.group(2), m.group(3)))
    out = [rxns[r] for r in order]
    # the aaTransport loop, expanded (external conc clamped, ABC form costs ATP)
    for symp, abc, met in AA_TRANSPORT:
        out.append(dict(id=symp, line=2497, subs=[], prods=[met],
                        clamped=[("Sub1", "medium_aa")]))
        out.append(dict(id=abc, line=2522, subs=["M_atp_c"],
                        prods=[met, "M_adp_c", "M_pi_c"],
                        clamped=[("Sub1", "medium_aa")]))
    return out

if __name__ == "__main__":
    rs = handcoded_reactions()
    print(f"{len(rs)} hand-coded ODE reactions")
    for r in rs[:24]:
        print(f"  {r['id']:<12} L{r['line']:<5} S={r['subs']} P={r['prods']} clamp={[c[1] for c in r['clamped']]}")

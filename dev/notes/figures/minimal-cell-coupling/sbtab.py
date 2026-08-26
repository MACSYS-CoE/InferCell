"""Minimal SBtab reader: split a multi-table TSV into named DataFrames."""
import pandas as pd, io, re

def read_sbtab(path):
    tables = {}
    with open(path) as fh:
        lines = fh.read().split("\n")
    starts = [i for i, l in enumerate(lines) if l.startswith("!!SBtab")]
    starts.append(len(lines))
    for a, b in zip(starts, starts[1:]):
        hdr = lines[a]
        m = re.search(r"TableName='([^']+)'", hdr)
        name = m.group(1) if m else re.search(r"TableType='([^']+)'", hdr).group(1)
        block = [l for l in lines[a+1:b] if l.strip()]
        if not block:
            continue
        cols = [c.strip().lstrip("!") for c in block[0].split("\t")]
        rows = [r.split("\t") for r in block[1:]]
        rows = [r + [""]*(len(cols)-len(r)) for r in rows]
        rows = [r[:len(cols)] for r in rows]
        df = pd.DataFrame(rows, columns=cols)
        tables[name] = df
    return tables

_ARROW = re.compile(r"\s*(<=>|=>|<=|->|<->)\s*")

def parse_formula(f):
    """ReactionFormula -> (substrates, products) as [(coeff, species)]."""
    parts = _ARROW.split(f.strip())
    lhs, rhs = parts[0], parts[2] if len(parts) > 2 else ""
    def side(s):
        out = []
        for term in s.split("+"):
            term = term.strip()
            if not term:
                continue
            bits = term.split()
            if len(bits) == 2 and re.fullmatch(r"[\d.]+", bits[0]):
                out.append((float(bits[0]), bits[1]))
            else:
                out.append((1.0, bits[-1]))
        return out
    return side(lhs), side(rhs)

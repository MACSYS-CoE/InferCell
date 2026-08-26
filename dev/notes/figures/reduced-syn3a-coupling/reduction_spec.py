"""Core A' reduction spec -- what fig 1 keeps, changes and deletes.

Every entry is traceable to dev/notes/reduced-syn3a-scoping.md. Nothing here is
derived from the model files: fig 1 is a view of state_interfaces.csv, this is a
hand-written overlay on top of it. See README.md, "What is derived and what is not".
"""

# --- palette ---------------------------------------------------------------
GREY_FILL = "#eceff1"          # removed module body
GREY_LINE = "#aab2b7"          # removed module border
GREY_TEXT = "#8a9399"          # removed module text
GREY_EDGE = "#9aa5ab99"        # removed edge: grey at 60% alpha, so it recedes
                               # behind the solid grey (#b3bfc4) currency edges
CHANGED_FILL = "#ffe9a8"       # survives in altered form; border colour unchanged
CROSS = "#c0392b"              # the red cross; fig 1's "clamped" red

# --- nodes -----------------------------------------------------------------
# Removed: greyed and crossed. Label is fig 1's, in grey.
REMOVED_NODES = {
    "CME:Replication": "<b>DNA replication</b>",
    "Lipid":     '<b>Lipid</b><br/><font point-size="8.5">17 rxns</font>',
    "Cofactor":  '<b>Cofactor</b><br/><font point-size="8.5">13 rxns</font>',
    "AminoAcid": '<b>Amino acid</b><br/><font point-size="8.5">1 rxns</font>',
}

# Changed: block border kept, fill goes amber, old -> new noted.
def _mod(name, old_new, note):
    return (f'<b>{name}</b><br/><font point-size="8.5"><b>{old_new}</b></font>'
            f'<br/><font point-size="7.5">{note}</font>')

CHANGED_NODES = {
    "Transport":  _mod("Transport",  "87 rxns → 6 rxns",
                       "PTS cascade ×5 · lactate export"),
    "Central":    _mod("Central",    "26 rxns → 10 rxns",
                       "glycolysis through lactate · NOX dropped"),
    "Nucleotide": _mod("Nucleotide", "46 rxns → 5 rxns",
                       "PGK3, PYK3 · ADK1, PPA, GK1"),
    "CME:Transcription":
        '<b>transcription</b><br/><font point-size="8.5"><b>455 → 17 genes</b></font>',
    "CME:mRNAdecay":
        '<b>mRNA decay</b><br/><font point-size="8.5"><b>455 → 17 mRNAs</b></font>',
    "CME:Translation":
        '<b>translation</b><br/><font point-size="8.5"><b>455 → 17 genes</b></font>',
    "CME:tRNAcharging":
        '<b>tRNA charging</b><br/><font point-size="8.5"><b>20 × 5 rxns → 1 lumped step</b>'
        '</font><br/><font point-size="7.5"><i>ours, not the model’s</i></font>',
}

# --- edges -----------------------------------------------------------------
# Greyed because an endpoint is removed. Keyed (source, target, discriminator);
# discriminator is a substring that picks one of a set of parallel edges, or "".
DEAD_ENDPOINT = [
    ("Cofactor", "AminoAcid", ""), ("AminoAcid", "Cofactor", ""),
    ("Transport", "Cofactor", ""), ("Central", "Cofactor", ""),
    ("Nucleotide", "Cofactor", ""), ("Transport", "Lipid", ""),
    ("Central", "Lipid", ""), ("Lipid", "Central", ""),
    ("Lipid", "Nucleotide", ""), ("Nucleotide", "Lipid", ""),
    ("Lipid", "CURRENCY", ""), ("Cofactor", "CURRENCY", ""),
    ("CME:Replication", "HOOK", ""), ("HOOK", "AminoAcid", ""),
    ("AminoAcid", "CME:Translation", ""), ("AminoAcid", "REBUILD", ""),
    ("REBUILD", "CME:Replication", ""),
    ("Nucleotide", "REBUILD", "datp"),          # the dNTP rebuild channel
    ("Lipid", "GROWTH", ""),
    ("CME:Replication", "CME:Transcription", ""),
    ("CME:Translation", "CME:Replication", ""),
    ("CME:Replication", "CME:Replication", ""),  # DNApol3 = 35, frozen
]

# Greyed although both endpoints survive: the coupling itself is gone.
DEAD_COUPLING = [
    (("Transport", "Nucleotide", ""),
     "the salvage network and its transport are cut entirely"),
    (("Nucleotide", "Transport", ""),
     "no retained transport reaction consumes pyruvate; PTS produces it"),
    (("CME:Transcription", "CME:tRNAcharging", ""),
     "the 17 genes are protein-coding; trna becomes a conserved pool"),
    (("CME:Translation", "CME:tRNAcharging", ""),
     "the lumped charging step carries no synthetase species"),
    (("HOOK", "Cofactor", ""),
     "fig 1's 'drawn once' catalytic edge happens to point at a deleted module; "
     "the channel itself survives -- see the note in the figure"),
]

# Surviving edges whose species list Core A' changes. (source, target, discrim):
# list of (old_substring, new_substring) applied in order.
RELABEL = {
    ("Central", "Nucleotide", ""):  [("13dpg, 2dr1p +4", "13dpg, pep")],
    ("Nucleotide", "Central", ""):  [("13dpg, 2dr1p +3", "3pg, pyr")],
    ("Central", "Transport", ""):   [("ac, lac +2", "lac, pep")],
    ("Transport", "CME:tRNAcharging", ""): [("20 aa pools", "aa pool")],
    ("CME:tRNAcharging", "CME:Translation", ""): [("20 charged tRNAs", "1 charged-tRNA pool")],
    ("CME:tRNAcharging", "REBUILD", ""): [("20 charged tRNAs", "charged-tRNA pool")],
    ("CME:Transcription", "CME:Translation", ""): [("455 mRNAs", "17 mRNAs")],
    ("CME:Transcription", "CME:mRNAdecay", ""):   [("455 mRNAs", "17 mRNAs")],
    ("CME:Translation", "HOOK", "fmettrna"):      [("atp, fmettrna, gtp", "atp, gtp")],
    ("CME:Translation", "GROWTH", ""): [("93 membrane protein loci",
                                         "1 membrane protein locus (ptsG)")],
    ("HOOK", "Nucleotide", "dNTP pools"): [("debit NTP / dNTP pools", "debit NTP pools")],
    ("HOOK", "Central", ""): [("mRNAdeg, DNArep, transloc", "mRNAdeg, transloc")],
    ("Nucleotide", "REBUILD", "atp, ctp"): [("atp, ctp, gtp, utp",
                                             "atp, gtp &#160;<i>ctp, utp chemostatted</i>")],
    ("REBUILD", "CME:Translation", ""): [("from live charged-tRNA counts",
                                          "from the live charged-tRNA pool")],
    ("GROWTH", "HOOK", ""): [("dilutes all six ODE modules", "dilutes all three ODE modules")],
}

#!/usr/bin/env python3
"""Regenerate the two Core A' phosphotransferase-transport extracts.

    python3 dev/scripts/extract_pts_transport.py <path-to-Minimal_Cell>

Writes, under src/organisms/coreA/data/:

    pts_transport.tsv           11 rate constants, no uncertainty column
    pts_initial_conditions.tsv   9 initial conditions, with their derivation

Re-running against the same checkout reproduces both files byte for byte;
src/organisms/coreA/data/README.md states the command and the upstream commit.

Nothing here reaches for a third-party package. The one binary input,
proteomics.xlsx, is read with zipfile and ElementTree, because openpyxl is not
in any environment on this machine and compute nodes have no network. That
reader is about twenty lines and is reused by spec phase 10, whose promoter
proxy is the same table's copy number over 180.

Every derived number is cross-checked against a second, independent upstream
statement before it is written. A silently wrong initial condition is the
failure this script exists to prevent.
"""
import csv
import math
import os
import re
import sys
import zipfile
import xml.etree.ElementTree as ET

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.normpath(os.path.join(HERE, "..", ".."))
OUT_DIR = os.path.join(REPO, "src", "organisms", "coreA", "data")

# Kept in step with src/handshake.jl's AVOGADRO and cell_volume_litres. The
# factor is derived rather than typed so that the scoping note's 20,180 is a
# reproduced number; note this is 20180.387, and the archived design's
# seven-decimal table used a rounded 20180, which is why four of its phospho
# values differ from ours in the seventh decimal. See the README.
AVOGADRO = 6.02214076e23
COREA_INITIAL_RADIUS_NM = 200.0


def particles_per_mM(radius_nm=COREA_INITIAL_RADIUS_NM):
    volume_litres = (4 / 3) * math.pi * (radius_nm * 1e-9) ** 3 * 1e3
    return 1e-3 * volume_litres * AVOGADRO


# The eleven constants, in the order the upstream Quantity table lists them.
RATE_IDS = [
    "KF_0_R_GLCpts0", "KR_0_R_GLCpts0",
    "KF_1_R_GLCpts1", "KR_1_R_GLCpts1",
    "KF_2_R_GLCpts2", "KR_2_R_GLCpts2",
    "KF_3_R_GLCpts3", "KR_3_R_GLCpts3",
    "KF_4_R_GLCpts4", "KR_4_R_GLCpts4",
    "P_R_L_LACt2r",
]

# Spot values the archived design records, asserted here so a changed upstream
# table fails loudly rather than importing a different model.
RATE_SPOT_CHECKS = {
    "KF_1_R_GLCpts1": 200000.0,
    "KR_4_R_GLCpts4": 1.00e-05,
    "P_R_L_LACt2r": 5.00e-09,
    "KF_4_R_GLCpts4": 0.88,
}

# The four carriers. `aoe` and `description` identify the proteomics row
# without needing the annotation workbook's MMSYN1 -> JCVISYN2 -> AOE chain;
# `copies` is what that row must round to. `total_mM` is the disabled block at
# setICs_two.py:286-288, an independent plain-text statement of the same
# quantity, used as the cross-check on the whole derivation.
CARRIERS = [
    dict(name="ptsI", gene="JCVISYN3A_0233", aoe="AOE93321.1",
         description="phosphoenolpyruvate--protein phosphotransferase",
         copies=353, unphos="M_ptsi_c", phos="M_ptsi_P_c",
         source="protein_metabolites_frac.csv", key="ptsi", total_mM=0.01750),
    dict(name="ptsH", gene="JCVISYN3A_0694", aoe="AOE93571.1",
         description="phosphocarrier protein HPr",
         copies=290, unphos="M_ptsh_c", phos="M_ptsh_P_c",
         source="protein_metabolites_frac.csv", key="ptsh", total_mM=0.01438),
    dict(name="Crr", gene="JCVISYN3A_0234", aoe="AOE93322.1",
         description="PTS glucose transporter subunit IIA",
         copies=314, unphos="M_crr_c", phos="M_crr_P_c",
         source="protein_metabolites_frac.csv", key="crr", total_mM=0.01557),
    dict(name="ptsG", gene="JCVISYN3A_0779", aoe="AOE93596.1",
         description="PTS sugar transporter",
         copies=831, unphos="M_ptsg_c", phos="M_ptsg_P_c",
         source="membrane_protein_metabolites.csv", key="ptsg", total_mM=0.04121),
]

# Registry order (src/organisms/coreA/registry.jl), so the extract and the
# registry list the eight phospho-states the same way.
REGISTRY_ORDER = ["M_ptsi_c", "M_ptsi_P_c", "M_ptsh_c", "M_ptsh_P_c",
                  "M_crr_c", "M_crr_P_c", "M_ptsg_c", "M_ptsg_P_c"]

TRANSPORT_TSV = os.path.join("CME_ODE", "model_data",
                             "transport_NoH2O_Zane-TB-DB.tsv")
PROTEOMICS_XLSX = os.path.join("CME_ODE", "model_data", "proteomics.xlsx")
SET_ICS = os.path.join("CME_ODE", "program", "setICs_two.py")
FRAC_CSVS = {
    "protein_metabolites_frac.csv":
        os.path.join("CME_ODE", "model_data", "protein_metabolites_frac.csv"),
    "membrane_protein_metabolites.csv":
        os.path.join("CME_ODE", "model_data", "membrane_protein_metabolites.csv"),
}


def fail(message):
    raise SystemExit("extract_pts_transport: " + message)


def read_quantity_table(path):
    """The transport file's Quantity table, as id -> (value, unit, line)."""
    rows = {}
    in_quantity = False
    with open(path) as handle:
        for lineno, line in enumerate(handle, start=1):
            if line.startswith("!!SBtab"):
                in_quantity = "TableType='Quantity'" in line
                continue
            if not in_quantity or line.startswith("!") or not line.strip():
                continue
            cells = line.rstrip("\n").split("\t")
            rows[cells[0].strip()] = (cells[2].strip(), cells[3].strip(), lineno)
    return rows


def read_compound_table(path):
    """The Compound table, as id -> (initial concentration, line)."""
    rows = {}
    in_compound = False
    with open(path) as handle:
        for lineno, line in enumerate(handle, start=1):
            if line.startswith("!!SBtab"):
                in_compound = "TableType='Compound'" in line
                continue
            if not in_compound or line.startswith("!") or not line.strip():
                continue
            cells = line.rstrip("\n").split("\t")
            rows[cells[0].strip()] = (cells[5].strip(), lineno)
    return rows


def read_fractions(root):
    """species -> (gene, proteomics_fraction, file, line), from both CSVs."""
    fractions = {}
    for logical, relative in FRAC_CSVS.items():
        path = os.path.join(root, relative)
        with open(path) as handle:
            for lineno, row in enumerate(csv.DictReader(handle), start=2):
                fractions[row["species"]] = (row["gene"],
                                             float(row["proteomics_fraction"]),
                                             logical, lineno)
    return fractions


def read_proteomics(path):
    """AOE id -> (description, copy number), from proteomics.xlsx.

    One sheet, header on row 2 (the published loaders pass skiprows=[0]).
    Column A is the AOE id, B the NCBI description, and V the absolute
    abundance in copies per cell -- which is the column those loaders reach
    positionally as iloc[:, 21].
    """
    ns = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
    archive = zipfile.ZipFile(path)
    shared = []
    if "xl/sharedStrings.xml" in archive.namelist():
        table = ET.fromstring(archive.read("xl/sharedStrings.xml"))
        shared = ["".join(t.text or "" for t in si.iter(ns + "t"))
                  for si in table.iter(ns + "si")]
    sheet = ET.fromstring(archive.read("xl/worksheets/sheet1.xml"))
    rows = {}
    for row in sheet.iter(ns + "row"):
        cells = {}
        for cell in row.iter(ns + "c"):
            value = cell.find(ns + "v")
            if value is None:
                continue
            column = re.match(r"[A-Z]+", cell.get("r")).group()
            cells[column] = (shared[int(value.text)] if cell.get("t") == "s"
                             else value.text)
        if "A" in cells and "V" in cells:
            rows[cells["A"]] = (cells.get("B", ""), cells["V"])
    return rows


def read_disabled_totals(path):
    """The four total concentrations in setICs_two.py's disabled block.

    Disabled code, so not the citation -- but an independent plain-text
    statement of the same quantity, and therefore the right cross-check on a
    derivation whose inputs are a spreadsheet and two CSVs.
    """
    text = open(path).read()
    totals = {}
    for key in ("ptsi", "ptsh", "crr", "ptsg"):
        match = re.search(r'"%s":([0-9.]+)\*' % key, text)
        if match is None:
            fail("setICs_two.py no longer carries a total for %s; the "
                 "cross-check on the copy numbers cannot be made" % key)
        totals[key] = float(match.group(1))
    return totals


def main():
    if len(sys.argv) != 2:
        raise SystemExit(__doc__)
    root = os.path.abspath(sys.argv[1])
    for relative in [TRANSPORT_TSV, PROTEOMICS_XLSX, SET_ICS] + list(FRAC_CSVS.values()):
        if not os.path.isfile(os.path.join(root, relative)):
            fail("upstream input has moved: %s is not under %s" % (relative, root))

    transport_path = os.path.join(root, TRANSPORT_TSV)
    quantities = read_quantity_table(transport_path)
    compounds = read_compound_table(transport_path)
    fractions = read_fractions(root)
    proteomics = read_proteomics(os.path.join(root, PROTEOMICS_XLSX))
    disabled_totals = read_disabled_totals(os.path.join(root, SET_ICS))
    factor = particles_per_mM()

    # --- the eleven rate constants ------------------------------------------
    missing = [i for i in RATE_IDS if i not in quantities]
    if missing:
        fail("the transport Quantity table no longer holds %s" % ", ".join(missing))
    for identifier, expected in RATE_SPOT_CHECKS.items():
        found = float(quantities[identifier][0])
        if found != expected:
            fail("%s reads %r upstream, not the recorded %r" %
                 (identifier, found, expected))

    placeholder = compounds["ptsi"][0]
    glucose_upstream = compounds["M_glc__D_e"][0]
    lactate_upstream = compounds["M_lac__L_e"][0]

    rate_lines = [
        "!!SBtab SBtabVersion='1.0' TableType='Quantity' "
        "TableName='pts transport constants' Document='coreA'",
        "% Derived extract. See src/organisms/coreA/data/README.md for the",
        "% upstream file, the commit and the regeneration command.",
        "% There is deliberately no GeometricStd column: the upstream file has",
        "% no Parameter table and no uncertainty anywhere, so read_source_table",
        "% derives informedness :asserted for every row from the table's shape.",
        "!ID\t!Mode\t!Unit\t!UpstreamRow",
    ]
    for identifier in RATE_IDS:
        value, unit, lineno = quantities[identifier]
        rate_lines.append("%s\t%s\t%s\t%s:%d" %
                          (identifier, value, unit or "-",
                           os.path.basename(TRANSPORT_TSV), lineno))

    # --- the nine initial conditions ----------------------------------------
    derived = {}
    provenance = {}
    for carrier in CARRIERS:
        description, raw = proteomics.get(carrier["aoe"], (None, None))
        if description is None:
            fail("proteomics.xlsx no longer holds %s (%s)" %
                 (carrier["aoe"], carrier["name"]))
        if not description.startswith(carrier["description"]):
            fail("proteomics row %s describes %r, not %r -- the sheet has "
                 "shifted and the copy numbers would be another protein's" %
                 (carrier["aoe"], description[:60], carrier["description"]))
        copies = round(float(raw))
        if copies != carrier["copies"]:
            fail("%s reads %s copies upstream, not the recorded %d" %
                 (carrier["name"], raw, carrier["copies"]))

        gene, unphos_frac, frac_file, frac_line = fractions[carrier["key"]]
        _, phos_frac, _, phos_line = fractions[carrier["key"] + "_P"]
        if frac_file != carrier["source"]:
            fail("%s's fractions moved to %s" % (carrier["name"], frac_file))
        if abs(unphos_frac + phos_frac - 1.0) > 1e-12:
            fail("%s's two proteomics fractions sum to %r, not 1" %
                 (carrier["name"], unphos_frac + phos_frac))
        if "JCVISYN3A_" + gene.split("_")[1] != carrier["gene"]:
            fail("%s maps to gene %s upstream, not %s" %
                 (carrier["name"], gene, carrier["gene"]))

        # The cross-check: our derivation against setICs_two.py's disabled
        # block, which states the same totals independently. They agree to
        # four decimals -- the block used a rounded 20180 particles per mM.
        total_mM = copies / factor
        published_total = disabled_totals[carrier["key"]]
        if abs(total_mM - published_total) > 1e-4:
            fail("%s derives %.6f mM from %d copies, but setICs_two.py's "
                 "disabled block says %.6f mM -- the two disagree by more "
                 "than the rounding in its conversion factor" %
                 (carrier["name"], total_mM, copies, published_total))
        if round(total_mM * factor) != copies:
            fail("%s does not round-trip: %.10f mM is %.4f copies" %
                 (carrier["name"], total_mM, total_mM * factor))

        derived[carrier["unphos"]] = copies * unphos_frac / factor
        derived[carrier["phos"]] = copies * phos_frac / factor
        provenance[carrier["unphos"]] = (copies, unphos_frac,
                                         "%s:%d" % (frac_file, frac_line))
        provenance[carrier["phos"]] = (copies, phos_frac,
                                       "%s:%d" % (frac_file, phos_line))

    ic_lines = [
        "!!SBtab SBtabVersion='1.0' TableType='Quantity' "
        "TableName='pts initial conditions' Document='coreA'",
        "% Derived extract. Each phospho-state is the published per-cell copy",
        "%% number times the published proteomics_fraction, at %.7f particles"
        % factor,
        "% per mM (a 200 nm sphere). Copies and ProteomicsFraction carry the two",
        "% published inputs, so a derived concentration traces back to both.",
        "%",
        "%% Not the transport reconstruction's uniform %s mM per form, which is"
        % placeholder,
        "%% %d copies behind every carrier -- %s times the proteomics counts for"
        % (round(2 * float(placeholder) * factor),
           ", ".join("%.1f" % (2 * float(placeholder) * factor / c["copies"])
                     for c in CARRIERS[:3])),
        "% ptsI, ptsH and Crr -- so it is a placeholder, not a measurement.",
        "!ID\t!Mode\t!Unit\t!Copies\t!ProteomicsFraction\t!UpstreamRow",
    ]
    for species in REGISTRY_ORDER:
        copies, fraction, row = provenance[species]
        ic_lines.append("conc_%s\t%s\tmM\t%d\t%s\t%s" %
                        (species, repr(derived[species]), copies,
                         repr(fraction), row))
    ic_lines.append(
        "conc_M_lac__L_e\t0.0\tmM\t-\t-\t%s:277" % os.path.basename(SET_ICS))

    os.makedirs(OUT_DIR, exist_ok=True)
    for name, lines in [("pts_transport.tsv", rate_lines),
                        ("pts_initial_conditions.tsv", ic_lines)]:
        with open(os.path.join(OUT_DIR, name), "w") as handle:
            handle.write("\n".join(lines) + "\n")
        print("wrote %s (%d rows)" %
              (os.path.join("src/organisms/coreA/data", name),
               sum(1 for line in lines
                   if not line.startswith(("!", "%")))))

    print("external glucose upstream: %s mM in the Compound table, superseded "
          "by 40 mM at %s:279" % (glucose_upstream, os.path.basename(SET_ICS)))
    print("external lactate upstream: %s mM in the Compound table, superseded "
          "by 0.0 at %s:277" % (lactate_upstream, os.path.basename(SET_ICS)))


if __name__ == "__main__":
    main()

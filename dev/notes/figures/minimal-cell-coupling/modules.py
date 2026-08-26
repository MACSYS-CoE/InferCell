"""Module membership, transcribed from CME_ODE/program/defMetRxns.py @ db048ac."""
import pandas as pd

MD = "mc/CME_ODE/model_data/"

CENTRAL_FILE = MD + "Central_AA_Zane_Balanced_direction_fixed_nounqATP.tsv"
NUCL_FILE    = MD + "Nucleotide_Kinetic_Parameters.tsv"
LIPID_FILE   = MD + "lipid_NoH2O_balanced_model.tsv"
TRANS_FILE   = MD + "transport_NoH2O_Zane-TB-DB.tsv"

# defMetRxns.py:133
CENTRAL = ["PGI","PFK","FBA","TPI","GAPD","PGK","PGM","ENO","PYK",
           "LDH_L","PDH_acald","PDH_E3","PTAr","ACKr","NOX","TALA","TKT1",
           "TKT2","RPE","RPI","PRPPS","PPM","PPM2","DRPA","GAPDP"]
# defMetRxns.py:288  (22 of 23 tRNA-synthetase reactions commented out)
AMINOACID = ["FMETTRS"]
# defMetRxns.py:916
COFACTOR = ["NCTPPRT","NNATr","NADS","NADK","RBFK","FMNAT",
            "5FTHFPGS","FTHFCL","MTHFC","GHMT","MTHFD"]
# defMetRxns.py:744-749
LIPID_ALL = ['GLYK','ACPS','BPNT','FAKr','ACPPAT','APG3PAT','AGPAT','DASYN',
             'PGSA','PGPP','CLPNS','PGMT','PAPA','GALU','UDPG4E','UDPGALM',
             'DAGGALT','PSSYN','DAGPST']
LIP_OFF = ['PAPA','DAGPST']            # lipTurnOff, defMetRxns.py:746
# defMetRxns.py:507-527
NUCL_OFF = ['NTD9','DCDPMP','DCTPMP','DCTPDP','CTPDP','PYK','PGK']

def module_reactions():
    nucl = [x.strip() for x in
            pd.read_csv(MD + "nucleo_rxns_list.txt", header=None)[0].tolist()]
    nucl = [x for x in nucl if x not in NUCL_OFF]
    lip  = [x for x in LIPID_ALL if x not in LIP_OFF]
    return {
        "Central":    ["R_" + r for r in CENTRAL],
        "Nucleotide": ["R_" + r for r in nucl],
        "Lipid":      ["R_" + r for r in lip],
        "Cofactor":   ["R_" + r for r in COFACTOR],
        "AminoAcid":  ["R_" + r for r in AMINOACID],
    }

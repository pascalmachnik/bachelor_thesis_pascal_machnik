import ROOT as r
import yaml
import argparse
from filter import df_filter
from define import BDT_define

# Argparse setup
parser = argparse.ArgumentParser()
parser.add_argument("--input", action="store", nargs="+", dest="input", required=True,type=str)
parser.add_argument("--output", action="store", dest="output", type=str, required=True)
parser.add_argument("--track_type", action="store", dest="track_type", type=str, choices=["ll", "dd"], required=True,)
parser.add_argument("--sample", action="store", dest="sample", type=str, required=True)
args = parser.parse_args()

# load config
with open("analysis_config.yml", "r") as config_file:
    config = yaml.safe_load(config_file)

# Enable multithreading for RDataFrame
r.EnableImplicitMT(8)

# names of branches to read and save
branches = [
    "L1_hasMuon", "L2_hasMuon",
    "L1_hasRich", "L2_hasRich",
    "L1_P", "L2_P",
    "L1_InAccMuon", "L2_InAccMuon",
    "Lambda_M", "Lambda_TAU", "Lambda_DIRA_OWNPV",
    "Lambda_FDCHI2_ORIVX", "Lambda_ENDVERTEX_Z",
    "L1_L0MuonDecision_TOS", "L2_L0MuonDecision_TOS", "Jpsi_L0DiMuonDecision_TOS",
    "Lambdab_Hlt1TrackMVADecision_TOS", "Lambdab_Hlt1TwoTrackMVADecision_TOS", "Lambdab_Hlt1TrackMuonMVADecision_TOS",
    "Lambdab_Hlt2Topo2BodyDecision_TOS", "Lambdab_Hlt2TopoMu2BodyDecision_TOS", "Lambdab_Hlt2TopoMuMu2BodyDecision_TOS",
    "Lambdab_Hlt2Topo3BodyDecision_TOS", "Lambdab_Hlt2TopoMu3BodyDecision_TOS", "Lambdab_Hlt2TopoMuMu3BodyDecision_TOS",
    "Lambdab_Hlt2DiMuonDetachedDecision_TOS", "Lambdab_Hlt2DiMuonDetachedHeavyDecision_TOS",
    "Lambdab_PE", "Lambda_PE", "Lambdab_PX", "Lambda_PX",
    "Lambdab_PY", "Lambda_PY", "Lambdab_PZ", "Lambda_PZ",
    "Lambdab_PT", "Proton_PT", "Pion_PT", "L1_PT", "L2_PT",
    "L1_ProbNNmu", "L2_ProbNNmu", "L1_IPCHI2_OWNPV", "L2_IPCHI2_OWNPV",
    "Pion_IPCHI2_OWNPV", "Proton_IPCHI2_OWNPV", "Lambda_PT", "Lambdab_DIRA_OWNPV",
    "Lambdab_LOKI_MASS_AllConstr", "Pion_TRACK_Type", "Proton_TRACK_Type", "Jpsi_M",
    "Lambdab_DiraCos", "Lambdab_IPCHI2_OWNPV", "Lambda_LOKI_DTF_CHI2NDOF", "Lambdab_ETA", 
    "Lambda_FD_OWNPV", "Pion_ETA", "nTracks", "L1_TRACK_Type", "L2_TRACK_Type"
]

if args.sample == "JpsiMC":
    branches.append("total_weights")
    branches.append("Lambdab_BKGCAT")
if args.sample == "JpKsMC":
    branches.append("Lambda_TRUEID")
    branches.append("Lambda_MC_MOTHER_ID")
    branches.append("Jpsi_MC_MOTHER_ID")
    branches.append("Lambdab_TRUEID")

# create RDataFrame for input files
df = r.RDataFrame("DecayTree", args.input, branches)

# print the number of events in the DataFrame
print(f"Number of events in the DataFrame for {args.sample} {args.track_type}: {df.Count().GetValue()}")

# apply filter and calculate BDT variables
if args.sample == "JpsiMC":
    df = df_filter(df, "truth_matching_Jpsi")
    df = df_filter(df, "trigger")
    df = df_filter(df, "preselection")
    df = df_filter(df, "Jpsi")
    df = df_filter(df, args.track_type)
    df = BDT_define(df)
elif args.sample == "JpKsMC":
    df = df_filter(df, "truth_matching_JpKs")
    df = df_filter(df, "trigger")
    df = df_filter(df, "preselection")
    df = df_filter(df, "Jpsi")
    df = df_filter(df, args.track_type)
    df = BDT_define(df)
elif args.sample == "Jpsi":
    df = df_filter(df, "trigger")
    df = df_filter(df, "preselection")
    df = df_filter(df, "Jpsi")
    df = df_filter(df, args.track_type)
    df = BDT_define(df)

# BDT features
BDT_features = config[f"ll_features"]
BDT_features.extend(config[f"dd_features"])

# names of branches to be saved
for feature in BDT_features:
    if feature not in branches:
        branches.append(feature)
branches.append("q2")

# print the number of events after filtering
print(f"Number of events after filtering for {args.sample} {args.track_type}: {df.Count().GetValue()}")

# save the DataFrame to a ROOT file
df.Snapshot("DecayTree", args.output, branches)
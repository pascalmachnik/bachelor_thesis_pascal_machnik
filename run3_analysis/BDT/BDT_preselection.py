import ROOT as r
import yaml
import argparse
import os
import sys
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__),".."))
sys.path.append(parent_dir)
from utils.filter import df_filter
from utils.define import BDT_define

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

# names of branches to read and save exluding the ones used for BDT training
branches = [
    "L1_PPHASRICH",
    "L2_PPHASRICH",
    "L1_PPHASMUONINFO",
    "L2_PPHASMUONINFO",
    "L1_PT",
    "L2_PT",
    "L1_P",
    "L2_P",
    "L1_INMUON",
    "L2_INMUON",
    "Lambdab_DTF_Lambda_CHI2DOF",
    "Lambda_M",
    "Lambda_BPVLTIME",
    "Lambda_BPVDIRA",
    "Lambda_END_VZ",
    "L1_PE",
    "L2_PE",
    "L1_PX",
    "L2_PX",
    "L1_PY",
    "L2_PY",
    "L1_PZ",
    "L2_PZ",
    "Lambdab_PT",
    "Proton_PT",
    "Pion_PT",
    "Proton_BPVIPCHI2",
    "Pion_BPVIPCHI2",
    "L1_BPVIPCHI2",
    "L2_BPVIPCHI2",
    "Lambda_PT",
    "Lambdab_BPVDIRA",
    "L1_PROBNN_MU",
    "L2_PROBNN_MU",
    "Lambdab_DTF_Lambda_PV_M",
    "Pion_TRACKTYPE",
    "Proton_TRACKTYPE",
    "Lambdab_Hlt1TrackMVADecision_TOS",
    "Lambdab_Hlt1TwoTrackMVADecision_TOS",
    "Lambdab_Hlt1TrackMuonMVADecision_TOS"
]

if args.sample == "MC":
    branches.append("Lambdab_BKGCAT")

# create RDataFrame for input files
df = r.RDataFrame("DecayTree", args.input, branches)

# apply filter and calculate BDT variables
if args.sample == "MC":
    df = df_filter(df, "truth_matching_rd")
    df = df_filter(df, "preselection")
    df = df_filter(df, args.track_type)
    df = BDT_define(df)
    df = df_filter(df, "high_q2")
    df = df_filter(df, "trigger")
else:
    df = df_filter(df, "preselection")
    df = df_filter(df, args.track_type)
    df = BDT_define(df)
    df = df_filter(df, "high_q2")
    df = df_filter(df, "trigger")

# BDT features
BDT_features = config[f"ll_features"]
BDT_features.extend(config[f"dd_features"])

# names of branches to be saved
for feature in BDT_features:
    if feature not in branches:
        branches.append(feature)
branches.append("q2")

# save the DataFrame to a ROOT file
df.Snapshot("DecayTree", args.output, branches)
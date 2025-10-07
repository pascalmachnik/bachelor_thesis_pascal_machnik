import ROOT as r
import yaml
from filter import df_filter
from define import BDT_define
import argparse


# Argparse setup
parser = argparse.ArgumentParser(description="Input and output file")
parser.add_argument("--signal", action="store", nargs="+", dest="signal", type=str, required=True)
parser.add_argument("--background", action="store", nargs="+", dest="background", type=str, required=True)
parser.add_argument("--output_signal_ll", action="store", dest="output_signal_ll", type=str, required=True)
parser.add_argument("--output_signal_dd", action="store", dest="output_signal_dd", type=str, required=True)
parser.add_argument("--output_background_ll", action="store", dest="output_background_ll", type=str,required=True)
parser.add_argument("--output_background_dd", action="store", dest="output_background_dd", type=str, required=True)
args = parser.parse_args()

# Root configurations
r.EnableImplicitMT()

# read the configuration file
with open("analysis_config.yml", "r") as config_file:
    config = yaml.safe_load(config_file)

# create a DataFrame from the input files
df_s = r.RDataFrame("DecayTree", args.signal)
df_b = r.RDataFrame("DecayTree", args.background)

# apply filter to the RDataFrames
df_s = df_filter(df_s, "truth_matching_rd")
df_s = df_filter(df_s, "trigger")
df_s = df_filter(df_s, "preselection_corr")

df_b = df_filter(df_b, "trigger")
df_b = df_filter(df_b, "preselection")
df_b = df_filter(df_b, "background")

# create longstream and downstream datasets
df_s_ll = df_filter(df_s, "ll")
df_s_dd = df_filter(df_s, "dd")
df_b_ll = df_filter(df_b, "ll")
df_b_dd = df_filter(df_b, "dd")

# calculate BDT variables
df_s_ll = BDT_define(df_s_ll, corrected=True)
df_s_dd = BDT_define(df_s_dd, corrected=True)
df_b_ll = BDT_define(df_b_ll)
df_b_dd = BDT_define(df_b_dd)

# Additional filter
df_s_ll = df_filter(df_s_ll, "high_q2")
df_b_ll = df_filter(df_b_ll, "high_q2")
df_s_dd = df_filter(df_s_dd, "high_q2")
df_b_dd = df_filter(df_b_dd, "high_q2")

# Cuts on BDT variables in order to improve the range of the BDT variables for training
df_s_dd = df_s_dd.Filter("Lambdab_LOKI_DTF_CHI2NDOF < 15")
df_b_dd = df_b_dd.Filter("Lambdab_LOKI_DTF_CHI2NDOF < 15")
df_s_dd = df_s_dd.Filter("Lambdab_TAU < 0.01")
df_b_dd = df_b_dd.Filter("Lambdab_TAU < 0.01")
df_s_dd = df_s_dd.Filter("Lambda_FDCHI2_OWNPV < 100000")
df_b_dd = df_b_dd.Filter("Lambda_FDCHI2_OWNPV < 100000")
df_s_dd = df_s_dd.Filter("Lambda_TAUCHI2 < 35")
df_b_dd = df_b_dd.Filter("Lambda_TAUCHI2 < 35")
df_s_dd = df_s_dd.Filter("Lambda_IPCHI2_OWNPV < 1000")
df_b_dd = df_b_dd.Filter("Lambda_IPCHI2_OWNPV < 1000")
df_s_dd = df_s_dd.Filter("Pion_IPCHI2_OWNPV < 9000")
df_b_dd = df_b_dd.Filter("Pion_IPCHI2_OWNPV < 9000")

df_s_ll = df_s_ll.Filter("Lambdab_TAU < 0.01")
df_b_ll = df_b_ll.Filter("Lambdab_TAU < 0.01")
df_s_ll = df_s_ll.Filter("Lambdab_LOKI_DTF_CHI2NDOF < 15")
df_b_ll = df_b_ll.Filter("Lambdab_LOKI_DTF_CHI2NDOF < 15")
df_s_ll = df_s_ll.Filter("Lambda_TAUCHI2 < 35")
df_b_ll = df_b_ll.Filter("Lambda_TAUCHI2 < 35")
df_s_ll = df_s_ll.Filter("Lambda_FDCHI2_OWNPV < 100000")
df_b_ll = df_b_ll.Filter("Lambda_FDCHI2_OWNPV < 100000")
df_s_ll = df_s_ll.Filter("Lambda_IPCHI2_OWNPV < 1000")
df_b_ll = df_b_ll.Filter("Lambda_IPCHI2_OWNPV < 1000")
df_s_ll = df_s_ll.Filter("Pion_IPCHI2_OWNPV < 10000")
df_b_ll = df_b_ll.Filter("Pion_IPCHI2_OWNPV < 10000")

print("Signal longstream events:", df_s_ll.Count().GetValue())
print("Signal downstream events:", df_s_dd.Count().GetValue())
print("Background longstream events:", df_b_ll.Count().GetValue())
print("Background downstream events:", df_b_dd.Count().GetValue())

# list of variables to be saved
branches_ll = config["ll_features"]
branches_dd = config["dd_features"]

# Checking for nans and infs for signal
df_b_ll = df_b_ll.Filter("&&".join([f"TMath::Finite({branch})" for branch in branches_ll]))
df_b_dd = df_b_dd.Filter("&&".join([f"TMath::Finite({branch})" for branch in branches_dd]))

# save the files
print("Saving files...")
df_b_ll.Snapshot("DecayTree", args.output_background_ll, branches_ll)
df_b_dd.Snapshot("DecayTree", args.output_background_dd, branches_dd)

# Add correction weights to the MC signal
branches_ll.append("total_weights")
branches_dd.append("total_weights")

# checking for nans and infs for background
df_s_ll = df_s_ll.Filter("&&".join([f"TMath::Finite({branch})" for branch in branches_ll]))
df_s_dd = df_s_dd.Filter("&&".join([f"TMath::Finite({branch})" for branch in branches_dd]))

df_s_ll.Snapshot("DecayTree", args.output_signal_ll, branches_ll)
df_s_dd.Snapshot("DecayTree", args.output_signal_dd, branches_dd)

print("Entries after NaN/Inf check s_ll:", df_s_ll.Count().GetValue())
print("Entries after NaN/Inf check s_dd:", df_s_dd.Count().GetValue())
print("Entries after NaN/Inf check b_ll:", df_b_ll.Count().GetValue())
print("Entries after NaN/Inf check b_dd:", df_b_dd.Count().GetValue())
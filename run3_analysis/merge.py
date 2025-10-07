from array import array
from utils.filter import df_filter
import ROOT as r
import argparse
import os

# Argparse setup
parser = argparse.ArgumentParser()
parser.add_argument("--input", action="store", nargs="+", dest="input", type=str, required=True)
parser.add_argument("--output", action="store", dest="output", type=str, required=True)
args = parser.parse_args()

# Enable multithreading for RDataFrame
r.EnableImplicitMT()

# dictionary for the block names
blocks = {
    "24c3/Hlt2RD_LbToLMuMu_LL/MagUp": 5,
    "24c3/Hlt2RD_LbToLMuMu_DD/MagUp": 5,
    "24c3/Hlt2RD_LbToLMuMu_LL/MagDown": 6,
    "24c3/Hlt2RD_LbToLMuMu_DD/MagDown": 6,
    "24c4/Hlt2RD_LbToLMuMu_LL/MagUp": 8,
    "24c4/Hlt2RD_LbToLMuMu_DD/MagUp": 8,
    "24c4/Hlt2RD_LbToLMuMu_LL/MagDown": 7,
    "24c4/Hlt2RD_LbToLMuMu_DD/MagDown": 7,
    "W35_37/MagUp": 5,
    "W37_39/MagDown": 6,
    "W40_42/MagDown": 7,
    "W40_42/MagUp": 8,
}

# create dfs from input files
dfs = []
sample = None
for file in args.input:
    for key, value in zip(blocks.keys(), blocks.values()):
        if key in file:
            if "MC" in file:
                sample = "MC"
                tree_dd = "Hlt2RD_LbToLMuMu_DD/DecayTree"
                tree_ll = "Hlt2RD_LbToLMuMu_LL/DecayTree"
                df_dd = r.RDataFrame(tree_dd, file)
                df_ll = r.RDataFrame(tree_ll, file)
                df_dd = df_dd.Define("block", f"int({value})")
                df_ll = df_ll.Define("block", f"int({value})")
                dfs.append(df_dd)
                dfs.append(df_ll)
                break
            else:
                sample = "data"
                df = r.RDataFrame("DecayTree", file)
                df = df.Define("block", f"int({value})")
                dfs.append(df)
                break
    
# Snapshot the dataframes into a temporary files
temp_files = []
for i, df in enumerate(dfs):
    temp_file = f"/ceph/users/pmachnik/temp/{sample}_Lb2Lmumu_temp_{i}.root"
    if os.path.exists(temp_file):
        os.remove(temp_file)
    df.Snapshot("DecayTree", temp_file)
    temp_files.append(temp_file)

# load the temporary file into a new dataframe
merged_df = r.RDataFrame("DecayTree", temp_files)

# save the merged dataframe to the output file
merged_df.Snapshot("DecayTree", args.output)

# remove the temporary files
for temp_file in temp_files:
    if os.path.exists(temp_file):
        os.remove(temp_file)

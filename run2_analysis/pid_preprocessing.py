import ROOT as r
from filter import df_filter
import argparse


# Argparse setup
parser = argparse.ArgumentParser(description="Input and output file")
parser.add_argument("--input", action="store", nargs="+", dest="input", type=str, required=True)
parser.add_argument("--output", action="store", dest="output", type=str, required=True)
args = parser.parse_args()

# Enable multithreading for RDataFrame
r.EnableImplicitMT()

# create a DataFrame from the input file
df = r.RDataFrame("DecayTree", args.input)

# aplly filter to the RDataFrame
df = df_filter(df, "truth_matching_rd")
df = df_filter(df, "fiducial")

# write the DataFrame to the output file
df.Snapshot("DecayTree", args.output)

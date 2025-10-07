import ROOT as r
import yaml
import copy
import argparse
import os
from define import df_define
from filter import df_filter


# Argparse setup
parser = argparse.ArgumentParser(description="Input and output file")
parser.add_argument("--input", action="store", nargs="+", dest="input", type=str, required=True)
parser.add_argument("--output", action="store", dest="output", type=str, required=True)
args = parser.parse_args()

# Enable multithreading for RDataFrame
r.EnableImplicitMT()

# create df from input file
df = r.RDataFrame("DecayTree", args.input)

# Create a new column with the L0 weights
df = df_define(df, "L0_weights")

# Create a new column with the lifetime weights
df = df_define(df, "life_time_weights")

# compute total weights
df = df_define(df, "total_weights")

# update the output file name
output_file = args.output
df.Snapshot("DecayTree", output_file)
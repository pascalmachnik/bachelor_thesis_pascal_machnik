import ROOT as r
import yaml
import argparse
import os
from filter import df_filter
from define import df_define, BDT_define

# Argparse setup
parser = argparse.ArgumentParser(description="Input and output file")
parser.add_argument("--output", action="store", dest="output", type=str, required=True)
parser.add_argument("--workspace", action="store", dest="workspace", type=str, required=True)
parser.add_argument("--track_type", action="store", dest="track_type", type=str, required=True)
args = parser.parse_args()

# Enable multithreading
r.ROOT.EnableImplicitMT()

### Computation of sWeights

# open the workspace
ws = r.TFile(args.workspace)
ws = ws.Get("ws_data_fit")
ws.loadSnapshot("fitResult")

# Extract the model, data and yields from the workspace
model = ws.pdf("model")
data = ws.data("data")
signal_yield = ws.var("signal_yield")
misID_yield = ws.var("misID_yield")
comb_yield = ws.var("comb_yield")

# Set all parameters except yields to constant
for var in model.getParameters(data):
    if not var.GetName().endswith("_yield"):
        var.setConstant(True)

# Use the SPlot tool to create sWeights
sWeights = r.RooStats.SPlot(
    "sWeights",
    "sWeights",
    data,
    model,
    r.RooArgList(signal_yield, misID_yield, comb_yield),
)

# Save the sWeights to the output file
temp_file = r.TFile(f"{args.track_type}_temp_sWeights.root", "RECREATE")
tree = data.GetClonedTree()
tree.SetName("DecayTree")
tree.SetDirectory(temp_file)
tree.Write()
temp_file.Close()


### Calculate the BDT variables and add them to the root file

# create a RDataFrame from the output file
df = r.RDataFrame("DecayTree", f"{args.track_type}_temp_sWeights.root")

# calcilate all BDT variables
df = BDT_define(df)

# save the sWeighted DataFrame to a root file
df = df.Snapshot("DecayTree", args.output)

# remove the temporary file
if os.path.exists(f"{args.track_type}_temp_sWeights.root"):
    os.remove(f"{args.track_type}_temp_sWeights.root")

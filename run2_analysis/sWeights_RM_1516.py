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

# Read configuration files
with open("analysis_config.yml", "r") as analysis_config_file:
    config = yaml.safe_load(analysis_config_file)

# Enable multithreading
r.ROOT.EnableImplicitMT()

### Computation of sWeights

# open the workspace
ws = r.TFile(args.workspace)
ws = ws.Get("ws_rm_data_fit")
ws.loadSnapshot("fitResult")

# Extract the model, data and yields from the workspace
model = ws.pdf("model")
data = ws.data("data")
signal_yield = ws.var("signal_yield")
#misID_yield = ws.var("misID_yield")
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
    r.RooArgList(signal_yield, comb_yield),
)

# Save the sWeights to the output file
temp_file = r.TFile(f"/ceph/users/pmachnik/temp/{args.track_type}_sWeights_1516.root", "RECREATE")
tree = data.GetClonedTree()
tree.SetName("DecayTree")
tree.SetDirectory(temp_file)
tree.Write()
temp_file.Close()

# Create a new RDataFrame from the temporary file
df = r.RDataFrame("DecayTree", temp_file.GetName())

# Calculate the BDT features
df = BDT_define(df)

# Snaoshot the df
df.Snapshot("DecayTree", args.output)

# delete the temporary file
if os.path.exists(temp_file.GetName()):
    os.remove(temp_file.GetName())


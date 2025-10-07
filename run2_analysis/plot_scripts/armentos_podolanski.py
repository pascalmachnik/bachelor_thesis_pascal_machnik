import ROOT as r
import yaml
import numpy as np
import argparse
import sys
import os
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.append(parent_dir)
from filter import df_filter

# Argparse setup
parser = argparse.ArgumentParser(description="Input and output file")
parser.add_argument("--input", action="store", nargs="+", dest="input", type=str, required=True)
parser.add_argument("--output", action="store", dest="output", type=str, required=True)
args = parser.parse_args()

# Enable multithreading for RDataFrame
r.EnableImplicitMT()

# read the configuration files
with open("plot_config.yml", "r") as plot_config_file:
    colors_hex = yaml.safe_load(plot_config_file)

# set global style configurations
r.gROOT.SetBatch(True)
r.gROOT.ProcessLine(".x lhcbStyle.C")
r.gStyle.SetOptStat(0)
r.gStyle.SetLegendFont(132)
r.gStyle.SetLegendTextSize(0.06)
r.gStyle.SetPalette(r.kBird)  # or kRainBow, kViridis, ...

colors = {key: r.TColor.GetColor(val) for key, val in colors_hex.items()}
signal2_fill_color = r.TColor.GetColorTransparent(colors["orange"], 0.7)  # green with alpha
signal1_fill_color = r.TColor.GetColorTransparent(colors["blue"], 0.7)  # blue with alpha

# create a RdataFrame from the input file
df = r.RDataFrame("DecayTree", args.input)

# aplly filters to the DataFrame
df = df_filter(df, "Jpsi")
df = df_filter(df, "preselection")
df = df_filter(df, "trigger")
df = df_filter(df, "fit_range")

# define the norm of the Lambda momentum
df = df.Define(
    "Lambda_P_norm",
    "sqrt(Lambda_PX * Lambda_PX + Lambda_PY * Lambda_PY + Lambda_PZ * Lambda_PZ)",
)
# define the longitudinal momentum of the pion
df = df.Define(
    "Pion_Mother_PL",
    "(Pion_PX * Lambda_PX + Pion_PY * Lambda_PY + Pion_PZ * Lambda_PZ) / Lambda_P_norm",
)

# define the transverse momentum of the pion
df = df.Define(
    "Pion_Mother_PT",
    "sqrt(Pion_PX * Pion_PX + Pion_PY * Pion_PY + Pion_PZ * Pion_PZ - Pion_Mother_PL * Pion_Mother_PL)",
)

# define the longitudinal momentum of the proton
df = df.Define(
    "Proton_Mother_PL",
    "(Proton_PX * Lambda_PX + Proton_PY * Lambda_PY + Proton_PZ * Lambda_PZ) / Lambda_P_norm",
)

# calculate the alpha parameter
df = df.Define(
    "alpha",
    "(Proton_Mother_PL - Pion_Mother_PL) / (Proton_Mother_PL + Pion_Mother_PL)",
)

# Now plot the two dimensional aremntos-podolanski distribution
canv = r.TCanvas("Armentos_Podolanski", "Armentos-Podolanski distribution", 800, 600)
canv.SetRightMargin(0.15)
canv.SetLeftMargin(0.15)

hist = df.Histo2D(
    ("Armentos_Podolanski", "Armentos-Podolanski distribution", 100, 0, 1, 100, 0, 200),
    "alpha",
    "Pion_Mother_PT",
)

hist.GetXaxis().SetTitle("#alpha")
hist.GetYaxis().SetTitle("p_{T}(#pi) [MeV/c]")

hist.Draw("COLZ")

# Save the canvas to the output file
canv.SaveAs(args.output)

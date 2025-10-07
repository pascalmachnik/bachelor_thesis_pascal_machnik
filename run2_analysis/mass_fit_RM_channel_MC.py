import ROOT as r
import yaml
import argparse
import os
from filter import df_filter
from define import df_define

# Argparse setup
parser = argparse.ArgumentParser(description="Input and output file")
parser.add_argument("--input", action="store", nargs="+", dest="input", type=str, required=True)
parser.add_argument("--output", action="store", dest="output", type=str, required=True)
parser.add_argument("--workspace", action="store", dest="workspace", type=str, required=True)
parser.add_argument("--track_type", action="store", dest="track_type", type=str, required=True)
args = parser.parse_args()

# Enable multithreading for RDataFrame
r.EnableImplicitMT()

# read the configuration file
with open("plot_config.yml", "r") as plot_config_file:
    colors_hex = yaml.safe_load(plot_config_file)

with open("analysis_config.yml", "r") as config_file:
    config = yaml.safe_load(config_file)

# set global style configurations
r.gROOT.SetBatch(True)
r.gROOT.ProcessLine(".x lhcbStyle.C")
r.gStyle.SetOptStat(0)
r.gStyle.SetLegendFont(132)
r.gStyle.SetLegendTextSize(0.06)
colors = {key: r.TColor.GetColor(val) for key, val in colors_hex.items()}
signal2_fill_color = r.TColor.GetColorTransparent(colors["orange"], 0.7)  # green with alpha
signal1_fill_color = r.TColor.GetColorTransparent(colors["blue"], 0.7)  # blue with alpha

# create df from input file
df = r.RDataFrame("DecayTree", args.input, ["Lambdab_LOKI_MASS_LConstr", "BDTG"])

# BDT selection
df_ll = df_filter(df, f"BDT_selection_{args.track_type}")


### Do a fit to the Lb mass distribution with RooFit ###

# Define the mass observable
mass = r.RooRealVar("Lambdab_LOKI_MASS_LConstr", "m_{#Lambda_{b}}", 5300, 5940)

# Define the signal model
# The signal is modeled as two double sided Crystal Ball functions
mean = r.RooRealVar("signal_mean", "signal_mean", 5620, 5500, 5740)
sigma1 = r.RooRealVar("siganl_sigma1", "signal_sigma1", 20, 1, 25)
alphal1 = r.RooRealVar("signal_alphaL1", "signal_alphaL1", 4, 2, 30)
alphar1 = r.RooRealVar("signal_alphaR1", "signal_alphaR1", 2, 0.1, 30)
nL1 = r.RooRealVar("signal_nL1", "signal_nL1", 1, 0, 5)
nR1 = r.RooRealVar("signal_nR1", "signal_nR1", 1, 0, 2)
signal1 = r.RooCrystalBall(
    "signal1", "First signal model", mass, mean, sigma1, alphal1, nL1, alphar1, nR1
)

# Define the second signal component
frac_sigma2 = r.RooRealVar("signal_frac_sigma2", "signal_frac_sigma2", 0.5, 0.0, 1.0)
sigma2 = r.RooFormulaVar("signal_sigma2", "@0 * @1", r.RooArgList(frac_sigma2, sigma1))
alphal2 = r.RooRealVar("signal_alphaL2", "signal_alphaL2", 2, 0.1, 30)
alphar2 = r.RooRealVar("signal_alphaR2", "signal_alphaR2", 2, 0.1, 30)
nL2 = r.RooRealVar("signal_nL2", "signal_nL2", 3, 0.1, 30)
nR2 = r.RooRealVar("signal_nR2", "signal_nR2", 0.1, 0, 6)
signal2 = r.RooCrystalBall(
    "signal2", "Second signal model", mass, mean, sigma2, alphal2, nL2, alphar2, nR2
)

# Combine the two signal components
signal_frac = r.RooRealVar("signal_frac", "Signal fraction", 0.5, 0.0, 1.0)
signal = r.RooAddPdf(
    "signal", "Signal model", r.RooArgList(signal1, signal2), r.RooArgList(signal_frac)
)

# Create a RooDataSet
df.Snapshot("DecayTree", f"/ceph/users/pmachnik/temp/{args.track_type}_temp_rm_mc.root", ["Lambdab_LOKI_MASS_LConstr"])
f = r.TFile(f"/ceph/users/pmachnik/temp/{args.track_type}_temp_rm_mc.root")
tree = f.Get("DecayTree")
data = r.RooDataSet("data", "Data", tree, r.RooArgSet(mass))


# Fit the model to the data
fit_result = signal.fitTo(data, r.RooFit.Save(), r.RooFit.PrintLevel(-1))


### Plot the fit result ###
# Create a canvas with two pads for the main plot and the pull plot
canvas = r.TCanvas("c", "", 1400, 1400)
pad1 = r.TPad("pad1", "pad1", 0, 0.3, 1, 1)
pad1.SetBottomMargin(0.02)
pad1.Draw()
pad2 = r.TPad("pad2", "pad2", 0, 0, 1, 0.3)
pad2.SetBottomMargin(0.4)
pad2.Draw()
padratio = 0.7 / 0.3

# Draw the fit on the  main pad
pad1.cd()
frame = mass.frame(r.RooFit.Bins(100))

data.plotOn(frame, r.RooFit.Invisible())  # Plot data invisibly first to set up the frame
# signal.plotOn(
# frame,
# r.RooFit.Components("signal1"),
# r.RooFit.FillColor(signal1_fill_color),
# r.RooFit.DrawOption("F"),
# r.RooFit.Name("Signal1"),
# )
# signal.plotOn(
# frame,
# r.RooFit.Components("signal2"),
# r.RooFit.FillColor(signal2_fill_color),
# r.RooFit.DrawOption("F"), 
# r.RooFit.Name("Signal2"),
# )
signal.plotOn(frame, r.RooFit.LineColor(colors["red"]), r.RooFit.Name("Signal"))
data.plotOn(frame, r.RooFit.Name("Data"))
r.gStyle.SetStatFontSize(0.06)  # Set the font size for the statistics box
signal.paramOn(frame, r.RooFit.Layout(0.18, 0.25, 0.9))

frame.Draw()

frame.GetXaxis().SetTitle("m(#Lambda^{0} #mu^{+}#mu^{-}) (MeV/c^{2})")
frame.GetXaxis().SetLabelSize(0)
frame.GetXaxis().SetTitleSize(0)

frame.GetYaxis().SetTitle("Events")
frame.GetYaxis().SetTitleOffset(1.2)
frame.GetYaxis().SetLabelSize(0.06)
frame.GetYaxis().SetTitleSize(0.06)

# frame.findObject("Signal1").SetLineColor(0)
# frame.findObject("Signal2").SetLineColor(0)

legend = r.TLegend(0.625, 0.7, 0.925, 0.9)
legend.AddEntry(frame.findObject("Data"), "Data", "lep")
legend.AddEntry(frame.findObject("Signal"), "Total Signal", "l")
# legend.AddEntry(frame.findObject("Signal1"), "dsCB1", "f")
# legend.AddEntry(frame.findObject("Signal2"), "dsCB2", "f")
legend.Draw()


# pull plot on the lower pad
pad2.cd()
pull_frame = mass.frame()

pull_hist = frame.pullHist("Data", "Signal")
pull_frame.addPlotable(pull_hist, "P")

pull_frame.GetXaxis().SetTitle("m_{#Lambda_{b}} [MeV/c^{2}]")
pull_frame.GetXaxis().SetLabelSize(padratio * 0.06)
pull_frame.GetXaxis().SetTitleSize(padratio * 0.06)

pull_frame.GetYaxis().SetTitle("Pull")
pull_frame.GetYaxis().SetTitleOffset(1.2 / padratio)
pull_frame.GetYaxis().SetLabelSize(padratio * 0.06)
pull_frame.GetYaxis().SetTitleSize(padratio * 0.06)
pull_frame.GetYaxis().SetNdivisions(505)
pull_frame.SetMinimum(-5.5)
pull_frame.SetMaximum(5.5)

pull_frame.Draw()

# Save the canvas to a file 
canvas.SaveAs(args.output)

# Create a RooWorkspace and import the model and data
ws = r.RooWorkspace("ws_mc_rm")
getattr(ws, "import")(mass)
getattr(ws, "import")(signal)
getattr(ws, "import")(data)
ws.saveSnapshot("fitResult", ws.allVars())
ws.writeToFile(args.workspace)

# remove the temporary file
os.remove(f"/ceph/users/pmachnik/temp/{args.track_type}_temp_rm_mc.root")

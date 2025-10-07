import ROOT as r
import yaml
import copy
import argparse
import os
from filter import df_filter


# Argparse setup
parser = argparse.ArgumentParser(description="Input and output file")
parser.add_argument("--input", action="store", nargs="+", dest="input", type=str, required=True)
parser.add_argument("--output", action="store", dest="output", type=str, required=True)
parser.add_argument("--Jpsi", action="store", dest="Jpsi", type=str, required=True)
parser.add_argument("--JpKs", action="store", dest="JpKs", type=str, required=True)
parser.add_argument("--workspace", action="store", dest="workspace", type=str, required=True)
parser.add_argument("--track_type", action="store", dest="track_type", type=str, required=True)
parser.add_argument("--BDTG", action="store_true", help="Apply BDTG cut")
args = parser.parse_args()

# Enable multithreading for RDataFrame
r.EnableImplicitMT()

# read the configuration files
with open("plot_config.yml", "r") as plot_config_file:
    colors_hex = yaml.safe_load(plot_config_file)

# read the worlspace parameters
ws_Jpsi_MC = r.TFile(args.Jpsi)
ws_Jpsi_MC = ws_Jpsi_MC.Get("ws_mc_jpsi")
ws_Jpsi_MC.loadSnapshot("fitResult")

ws_JpKs_MC = r.TFile(args.JpKs)
ws_JpKs_MC = ws_JpKs_MC.Get("ws_mc_jpks")
ws_JpKs_MC.loadSnapshot("fitResult")

# set global style configurations
r.gROOT.SetBatch(True)
r.gROOT.ProcessLine(".x lhcbStyle.C")
r.gStyle.SetOptStat(0)
r.gStyle.SetLegendFont(132)
r.gStyle.SetLegendTextSize(0.06)

colors = {key: r.TColor.GetColor(val) for key, val in colors_hex.items()}
signal_fill_color = r.TColor.GetColorTransparent(colors["orange"], 1.0)
comb_fill_color = r.TColor.GetColorTransparent(colors["blue"], 1.0)
misID_fill_color = r.TColor.GetColorTransparent(colors["pink"], 1.0)


# create df from input files
df = r.RDataFrame("DecayTree", args.input)

# Apply filter to the DataFrame
df = df_filter(df, "Jpsi")
df = df_filter(df, "trigger")
df = df_filter(df, "preselection")
df = df_filter(df, "fit_range")
df = df_filter(df, f"{args.track_type}")
if args.BDTG:
    df = df_filter(df, f"BDT_selection_{args.track_type}")


### Do a fit to the Lb mass distribution with RooFit ### 

# Define the mass observable
mass = r.RooRealVar("Lambdab_LOKI_MASS_AllConstr", "m_{#Lambda_{b}}", 5300, 5940)

# Define the signal model
# The signal is modeled as two double sided Crystal Ball functions
signal_mean = r.RooRealVar("signal_mean", "signal_mean", 5620, 5500, 5740)
signal_sigma = r.RooRealVar("signal_sigma", "signal_sigma", 5.7, 1, 20)
# Load parameters from the Jpsi MC workspace
signal_alphaL1 = ws_Jpsi_MC.var("signal_alphaL1")
signal_alphaR1 = ws_Jpsi_MC.var("signal_alphaR1")
signal_nL1 = ws_Jpsi_MC.var("signal_nL1")
signal_nR1 = ws_Jpsi_MC.var("signal_nR1")
# Set the external parameters to constant
signal_alphaL1.setConstant(True)
signal_alphaR1.setConstant(True)
signal_nL1.setConstant(True)
signal_nR1.setConstant(True)
signal1 = r.RooCrystalBall(
    "signal1", "First signal model", mass, signal_mean, signal_sigma, signal_alphaL1, signal_nL1, signal_alphaR1, signal_nR1
)

# Read the paramters for the second signal component
signal_frac_sigma2 = ws_Jpsi_MC.var("signal_frac_sigma2")
signal_sigma2 = r.RooFormulaVar(
    "signal_sigma2",
    "@0 * @1",
    r.RooArgList(signal_frac_sigma2, signal_sigma),
)
signal_alphaL2 = ws_Jpsi_MC.var("signal_alphaL2")
signal_alphaR2 = ws_Jpsi_MC.var("signal_alphaR2")
signal_nL2 = ws_Jpsi_MC.var("signal_nL2")
signal_nR2 = ws_Jpsi_MC.var("signal_nR2")
# Set the external parameters to constant
signal_frac_sigma2.setConstant(True)
signal_alphaL2.setConstant(True)
signal_alphaR2.setConstant(True)
signal_nL2.setConstant(True)
signal_nR2.setConstant(True)
signal2 = r.RooCrystalBall(
    "signal2", "Second signal model", mass, signal_mean, signal_sigma2, signal_alphaL2, signal_nL2, signal_alphaR2, signal_nR2
)

# Combine the two signal components
signal_frac = ws_Jpsi_MC.var("signal_frac")
signal_frac.setConstant(True)
signal = r.RooAddPdf(
    "signal", "Signal model", r.RooArgList(signal1, signal2), r.RooArgList(signal_frac)
)

# Define the combinatorial background model
# The combinatorial background is modeled as an exponential function
tau = r.RooRealVar("comb_tau", "comb_tau", -0.001, -1, 0)
comb_bkg = r.RooExponential("comb_bkg", "Combinatorial background", mass, tau)

# Define the model for the misID background
# The misID background is modeled as one double sided Crystal Ball function
misID_mean = r.RooRealVar("misID_mean", "misID_mean", 5450, 5300, 5740)
misID_sigma = r.RooRealVar("misID_sigma", "misID_sigma", 10, 0.5, 30)
misID_alphaL = ws_JpKs_MC.var("misID_alphaL")
misID_alphaR = ws_JpKs_MC.var("misID_alphaR")
misID_nL = ws_JpKs_MC.var("misID_nL")
misID_nR = ws_JpKs_MC.var("misID_nR")
misID_alphaL.setConstant(True)
misID_alphaR.setConstant(True)
misID_nL.setConstant(True)
misID_nR.setConstant(True)
misID_bkg = r.RooCrystalBall(
    "misID_bkg",
    "MisID background model",
    mass,
    misID_mean,
    misID_sigma,
    misID_alphaL,
    misID_nL,
    misID_alphaR,
    misID_nR,
)

# Define yields for signal and background
signal_yield = r.RooRealVar("signal_yield", "Signal yield", 1000, 0, 1e6)
combinatorial_yield = r.RooRealVar("comb_yield", "Combinatorial yield", 10000, 0, 1e6)
misID_yield = r.RooRealVar("misID_yield", "MisID yield", 5000, 0, 1e6)
yields = r.RooArgList(signal_yield, combinatorial_yield, misID_yield)

# combine the signal and background models
model = r.RooAddPdf(
    "model",
    "Total model",
    r.RooArgList(signal, comb_bkg, misID_bkg),
    yields,
)

# Store the data in a RooDataSet
# For the fit only the mass is needed, but for further analysis,
# more variables are needed

branches = [
    "Lambdab_TAU",
    "Lambdab_TAUCHI2",
    "Lambdab_ENDVERTEX_CHI2",
    "Lambdab_LOKI_DTF_CHI2NDOF",
    "Lambdab_DIRA_OWNPV",
    "Lambdab_DiraCos",
    "Lambdab_IPCHI2_OWNPV",
    "Lambda_TAU",
    "Lambda_TAUCHI2",
    "Lambda_FDCHI2_OWNPV",
    "Lambda_ENDVERTEX_Z",
    "Lambda_LOKI_DTF_CHI2NDOF",
    "Lambda_IPCHI2_OWNPV",
    "Lambda_FD_TOPPV",
    "Lambda_PT",
    "Lambdab_ETA",
    "Lambda_FD_OWNPV",
    "Pion_PT",
    "Pion_IPCHI2_OWNPV",
    "Pion_ETA",
    "Proton_PT",
    "Proton_IPCHI2_OWNPV",
    "L1_PT",
    "L2_PT",
    "nTracks",
    "Lambdab_PE",
    "Lambdab_PX",
    "Lambdab_PY",
    "Lambdab_PZ",
    "Lambda_PE",
    "Lambda_PX",
    "Lambda_PY",
    "Lambda_PZ",
    "Lambdab_PT",
    "L1_ProbNNmu",
    "L2_ProbNNmu",
    "L1_IPCHI2_OWNPV",
    "L2_IPCHI2_OWNPV",
    "Pion_TRACK_Type",
    "Proton_TRACK_Type",
    "L1_TRACK_Type",
    "L2_TRACK_Type",
    "Jpsi_L0DiMuonDecision_TOS",
    "L1_L0MuonDecision_TOS"
]

# create Argset for the branches
vars = [r.RooRealVar(branch, branch, -1e16,1e16) for branch in branches]
vars.append(mass)  # Add the mass variable to the list
# Add mass to the branches list
branches.append(mass.GetName()) 

# Create a RooDataSet
df.Snapshot("DecayTree", f"/ceph/users/pmachnik/temp/{args.BDTG}_{args.track_type}_temp_data.root", branches).GetValue()
f = r.TFile(f"/ceph/users/pmachnik/temp/{args.BDTG}_{args.track_type}_temp_data.root")
tree = f.Get("DecayTree")
data = r.RooDataSet("data", "Data", tree, r.RooArgSet(*vars))

# Fit the model to the data
fit_result = model.fitTo(data, r.RooFit.Save(), r.RooFit.PrintLevel(-1))


### Plot the fit result ###

canvas = r.TCanvas("c", "", 1500, 1400)
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
model.plotOn(
    frame,
    r.RooFit.Components("comb_bkg"),
    r.RooFit.FillColor(comb_fill_color),
    r.RooFit.DrawOption("F"),
    r.RooFit.Name("Combinatorial"),
)
model.plotOn(
    frame,
    r.RooFit.Components("signal"),
    r.RooFit.FillColor(signal_fill_color),
    r.RooFit.DrawOption("F"),
    r.RooFit.Name("Signal"),
)
model.plotOn(
    frame,
    r.RooFit.Components("misID_bkg"),
    r.RooFit.FillColor(misID_fill_color),
    r.RooFit.DrawOption("F"),
    r.RooFit.Name("MisID"),
)
model.plotOn(frame, r.RooFit.LineColor(colors["red"]), r.RooFit.Name("Total Model"))
r.gStyle.SetStatFontSize(0.06)  # Set the font size for the statistics box
model.paramOn(frame, r.RooFit.Layout(0.18, 0.25, 0.9))
data.plotOn(frame, r.RooFit.Name("Data"))

frame.Draw()

frame.GetYaxis().SetTitle("Events")
frame.GetXaxis().SetLabelSize(0)
frame.GetXaxis().SetTitleSize(0)

frame.GetXaxis().SetTitle("m_{#Lambda_{b}} (MeV/c^{2})")
frame.GetYaxis().SetTitleOffset(1.2)
frame.GetYaxis().SetLabelSize(0.06)
frame.GetYaxis().SetTitleSize(0.06)

frame.findObject("Combinatorial").SetLineColor(0)
frame.findObject("MisID").SetLineColor(0)
frame.findObject("Signal").SetLineColor(0)

legend = r.TLegend(0.6, 0.6, 0.9, 0.9)
legend.AddEntry(frame.findObject("Signal"), "Signal", "f")
legend.AddEntry(frame.findObject("Combinatorial"), "Combinatorial", "f")
legend.AddEntry(frame.findObject("MisID"), "MisID", "f")
legend.AddEntry(frame.findObject("Data"), "Data", "lep")
legend.AddEntry(frame.findObject("Total Model"), "Total Model", "l")
legend.Draw()


# pull plot on the lower pad  
pad2.cd()
pull_frame = mass.frame()

pull_hist = frame.pullHist("Data", "Total Model")
pull_frame.addPlotable(pull_hist, "P")

pull_frame.GetXaxis().SetTitle("m(#Lambda^{0} J/#psi) [MeV/c^{2}]")
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

# Write the fit in a workspace
ws = r.RooWorkspace("ws_data_fit")
getattr(ws, "import")(model)
getattr(ws, "import")(data)
ws.saveSnapshot("fitResult", ws.allVars())
ws.writeToFile(args.workspace)

# Clean up the temporary file
if os.path.exists(f"/ceph/users/pmachnik/temp/{args.BDTG}_{args.track_type}_temp_data.root"):
    os.remove(f"/ceph/users/pmachnik/temp/{args.BDTG}_{args.track_type}_temp_data.root")

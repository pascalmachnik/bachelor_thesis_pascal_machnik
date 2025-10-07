import ROOT as r
import yaml
import os
import sys
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__),".."))
sys.path.append(parent_dir) 
from filter import df_filter

############### Global Configurations ################

# Enable multithreading for RDataFrame
r.EnableImplicitMT(16)

# read the configuration files
with open("plot_config.yml", "r") as plot_config_file:
    colors_hex = yaml.safe_load(plot_config_file)


# set global style configurations
r.gROOT.SetBatch(True)
r.gROOT.ProcessLine(".x lhcbStyle.C")
r.gStyle.SetOptStat(0)
r.gStyle.SetLegendFont(132)
r.gStyle.SetLegendTextSize(0.06)
colors_list = {key: r.TColor.GetColor(val) for key, val in colors_hex.items()}


############### Plotting Functions ################

def plot_feature(dfs, features, xlabel, ylabel, xmin, xmax, nbins, output_file, styles, markers, colors, names = None, marker_colors = None,
                 weights = None, xwidth = 800, yheight = 600, legend_position = [0.6, 0.7, 0.8, 0.9], log_scale = False, normalize = False):
    """
    Plots histograms for given features from multiple DataFrames.
    
    Parameters:
    - dfs: List of DataFrames to plot
    - features: List of features to plot
    - xlabel: Label for x-axis
    - ylabel: Label for y-axis
    - xmin: Minimum x value for the histogram
    - xmax: Maximum x value for the histogram
    - nbins: Number of bins for the histogram
    - output_file: File path to save the plot
    - styles: List of styles for the histograms (e.g., "HIST", "E")
    - marker_styles: List of marker styles for the histograms
    - colors: List of colors for the histograms
    - names: List of names for the histograms. If None, no legend will be created.
    - marker_colors: Optional list of colors for the markers (default None)
    - weighrs: Optional list of Weights for the histograms (default None)
    - xwidth: Width of the canvas (default 800)
    - yheight: Height of the canvas (default 600)
    - legend_position: Coordinates for the legend (default 0.7, 0.7, 0.9, 0.9)
    - log_scale: Whether to set the y-axis to log scale (default False)
    - normalize: Whether to normalize the histograms (default False)
    """
    
    # if sweights are provided, use them for the histograms
    if weights is not None:
        histograms = []
        for df, feature, weight in zip(dfs, features, weights):
            if weight is not None:
                hist = df.Histo1D((feature, "", nbins, xmin, xmax), feature, weight)
            else:
                hist = df.Histo1D((feature, "", nbins, xmin, xmax), feature)
            histograms.append(hist)
    else:
        histograms = [df.Histo1D((feature, "", nbins, xmin, xmax), feature) for df, feature in zip(dfs, features)]
    
    # set histogram styles and colors
    for hist, marker, color in zip(histograms, markers, colors):
        hist.SetLineColor(colors_list[color])
        hist.SetMarkerStyle(marker)
        if marker_colors is not None:
            hist.SetMarkerColor(colors_list[marker_colors[histograms.index(hist)]])
    
    # Create canvas and draw histograms
    c = r.TCanvas("c", "Canvas", xwidth, yheight)
    c.SetTopMargin(0.05) 
    # Normalize histograms if specified
    if normalize:
        for hist in histograms:
            hist.Scale(1.0 / hist.Integral())

    # Set log scale if specified
    if log_scale:
        c.SetLogy()

    # Draw the histograms
    max_values = [hist.GetMaximum() for hist in histograms]

    sorted_indices = sorted(range(len(max_values)), key=lambda i: max_values[i], reverse=True)

    histograms[sorted_indices[0]].Draw(styles[sorted_indices[0]])
    histograms[sorted_indices[0]].GetXaxis().SetTitle(xlabel)
    histograms[sorted_indices[0]].GetXaxis().SetTitleSize(0.06)
    histograms[sorted_indices[0]].GetXaxis().SetLabelSize(0.06)
    histograms[sorted_indices[0]].GetXaxis().SetTitleOffset(1.2)
    histograms[sorted_indices[0]].GetYaxis().SetTitle(ylabel)
    histograms[sorted_indices[0]].GetYaxis().SetTitleOffset(1.2)
    histograms[sorted_indices[0]].GetYaxis().SetLabelSize(0.06)
    histograms[sorted_indices[0]].GetYaxis().SetTitleSize(0.06)
    
    for i in sorted_indices:
        if i == sorted_indices[0]:
            continue
        else:
            histograms[i].Draw(styles[i] + " SAME")

    # dict for legend markers
    legend_markers = {
        "HIST": "l",
        "E": "lep"
    }
    
    # Add legend
    if names is not None:
        legend = r.TLegend(*legend_position)

        for hist, name, style in zip(histograms, names, styles):
            legend.AddEntry(hist.GetValue(), name, legend_markers[style])
        
        legend.Draw()
    
    # Save the plot
    c.SaveAs(output_file)


def plot_feature_and_pull(dfs, features, xlabel, ylabel, xmin, xmax, nbins, output_file, styles, markers, colors, pull_index, names=None,
                 weights = None, xwidth = 800, yheight = 600, legend_position = [0.55, 0.7, 0.8, 0.9], log_scale = False, normalize = False,
                 pad1_height = 0.7, marker_colors = None):
    """
    Plots histograms for given features from multiple DataFrames.
    
    Parameters:
    - dfs: List of DataFrames to plot
    - features: List of features to plot
    - xlabel: Label for x-axis
    - ylabel: Label for y-axis
    - xmin: Minimum x value for the histogram
    - xmax: Maximum x value for the histogram
    - nbins: Number of bins for the histogram
    - output_file: File path to save the plot
    - styles: List of styles for the histograms (e.g., "HIST", "E")
    - marker_styles: List of marker styles for the histograms
    - colors: List of colors for the histograms
    - pull_index: List of indices for the histograms to be used for the pull plot ([index_data, index_model])
    - names: List of names for the histograms. If None, no legend will be created.
    - weights: Optional list of sWeights for the histograms (default None)
    - xwidth: Width of the canvas (default 800)
    - yheight: Height of the canvas (default 600)
    - legend_position: Coordinates for the legend (default 0.7, 0.7, 0.9, 0.9)
    - log_scale: Whether to set the y-axis to log scale (default False)
    - normalize: Whether to normalize the histograms (default False)
    - pad1_height: Height of the upper pad for the main histograms (default 0.7)
    - marker_colors: Optional list of colors for the markers (default None)
    """
    
    # if sweights are provided, use them for the histograms
    if weights is not None:
        histograms = []
        for df, feature, weight in zip(dfs, features, weights):
            if weight is not None:
                hist = df.Histo1D((feature, "", nbins, xmin, xmax), feature, weight)
            else:
                hist = df.Histo1D((feature, "", nbins, xmin, xmax), feature)
            histograms.append(hist)
    else:
        histograms = [df.Histo1D((feature, "", nbins, xmin, xmax), feature) for df, feature in zip(dfs, features)]
    
    # set histogram styles and colors
    for hist, marker, color in zip(histograms, markers, colors):
        hist.SetLineColor(colors_list[color])
        hist.SetMarkerStyle(marker)
        if marker_colors is not None:
            hist.SetMarkerColor(colors_list[marker_colors[histograms.index(hist)]])
    
    # Create canvas and pads
    c = r.TCanvas("c", "Canvas", xwidth, yheight)
    pad1 = r.TPad("pad1", "pad1", 0, 1 - pad1_height, 1, 1)
    pad1.SetBottomMargin(0.02)  # Set bottom margin to 0.02
    pad1.Draw()
    pad1.SetTopMargin(0.05) 
    pad2 = r.TPad("pad2", "pad2", 0, 0, 1, 1 - pad1_height)
    pad2.SetBottomMargin(0.3)  # Set bottom margin to 0.3
    pad2.Draw()
    padratio = pad1_height / (1-pad1_height)

    # Normalize histograms if specified
    if normalize:
        for hist in histograms:
            hist.Scale(1.0 / hist.Integral())

    # Draw the histograms
    pad1.cd()

    # Set log scale if specified
    if log_scale:
        r.gPad.SetLogy()

    max_values = [hist.GetMaximum() for hist in histograms]

    sorted_indices = sorted(range(len(max_values)), key=lambda i: max_values[i], reverse=True)

    histograms[sorted_indices[0]].Draw(styles[sorted_indices[0]])

    histograms[sorted_indices[0]].GetXaxis().SetTitle(xlabel)
    histograms[sorted_indices[0]].GetXaxis().SetTitleSize(0)
    histograms[sorted_indices[0]].GetXaxis().SetLabelSize(0)

    histograms[sorted_indices[0]].GetYaxis().SetTitle(ylabel)
    histograms[sorted_indices[0]].GetYaxis().SetTitleOffset(1.2)
    histograms[sorted_indices[0]].GetYaxis().SetLabelSize(0.06)
    histograms[sorted_indices[0]].GetYaxis().SetTitleSize(0.06)

    for i in sorted_indices:
        if i == sorted_indices[0]:
            continue
        else:
            histograms[i].Draw(styles[i] + " SAME")

    # dict for legend markers
    legend_markers = {
        "HIST": "l",
        "E": "lep"
    }
    
    # Add legend if names is provided
    if names is not None:
        legend = r.TLegend(*legend_position)

        for hist, name, style in zip(histograms, names, styles):
            legend.AddEntry(hist.GetValue(), name, legend_markers[style])

        legend.Draw()

    # Pull plot
    pad2.cd()

    hist_data = histograms[pull_index[0]].Clone("hist_data")
    hist_model = histograms[pull_index[1]].Clone("hist_model")
    hist_pull = hist_data.Clone("hist_pull")
    hist_pull.Reset()

    for i in range(1, hist_data.GetNbinsX() + 1):
        d = hist_data.GetBinContent(i)
        m = hist_model.GetBinContent(i)
        ed = hist_data.GetBinError(i)
        em = hist_model.GetBinError(i)
        denom = (ed**2 + em**2) ** 0.5
        if denom > 0:
            pull = (d - m) / denom
        else:
            pull = 0
        hist_pull.SetBinContent(i, pull)
        hist_pull.SetBinError(i, 0)

    hist_pull.SetLineColor(r.kBlack)
    #hist_pull.SetFillColor(r.kBlack)
    hist_pull.SetMarkerStyle(20)

    hist_pull.GetXaxis().SetTitle(xlabel)
    hist_pull.GetXaxis().SetLabelSize(padratio * 0.06)
    hist_pull.GetXaxis().SetTitleSize(padratio * 0.06)

    hist_pull.GetYaxis().SetTitle("Pull")
    hist_pull.GetYaxis().SetTitleOffset(1.2 / padratio)
    hist_pull.GetYaxis().SetLabelSize(padratio * 0.06)
    hist_pull.GetYaxis().SetTitleSize(padratio * 0.06)
    hist_pull.GetYaxis().SetNdivisions(505)

    hist_pull.SetMinimum(-5.5)
    hist_pull.SetMaximum(5.5)
    
    hist_pull.Draw("P")
    
    # Save the plot
    c.SaveAs(output_file)

############### Main Script for Plotting Features ################

#f_data = [
#    f"/ceph/users/jnicolini/scripts/Lambdab/Samples/angles/Lb2Lmm_Data_{year}_{magnet}.root"
#    for year in ["2015", "2016", "2017", "2018"]
#    for magnet in ["MU", "MD"]
#]
#
#df = r.RDataFrame("DecayTree", f_data)
#df = df_filter(df, "preselection")
#df = df_filter(df, "trigger")
##df = df_filter(df, "background")
#df = BDT_define(df)
#df = df_filter(df, "high_q2")
#
#plot_feature(
#    dfs = [df],
#    features=["Lambdab_LOKI_MASS_AllConstr"],
#    xlabel= "Lambdab LOKI Mass AllConstr",
#    ylabel= "Events",
#    xmin=4000,
#    xmax=7000,
#    nbins=100,
#    output_file="../plots/Lambdab_LOKI_MASS_AllConstr.pdf",
#    colors=["red"],
#    styles=["HIST"],
#    markers=[20],
#)

# L0 trigger weights
#f_MC = [
#    f"/ceph/users/jnicolini/scripts/Lambdab/Samples/angles/Lb2Lmm_PHMC_{year}_{magnet}.root" 
#    for year in config["years"] 
#    for magnet in config["magnet"]
#    ]
#df_MC = r.RDataFrame("DecayTree", f_MC)
#
#f_MC_corr = [
#    f"/ceph/users/pmachnik/Lambdab/Lb2Lmm_PHMC_{year}_{mag}_corr_all.root"
#    for year in config["years"]
#    for mag in config["magnet"]
#]
#df_MC_corr = r.RDataFrame("DecayTree", f_MC_corr)
#
#f_sw = "/ceph/users/pmachnik/Lambdab/Jpsi_sWeights.root"
#df_sw = r.RDataFrame("DecayTree", f_sw)
#
#df_MC = df_filter(df_MC, "truth_matching_rd")
#df_MC = df_filter(df_MC, "trigger")
#df_MC = df_filter(df_MC, "preselection")
#df_MC = df_filter(df_MC, "Jpsi")
#df_MC = df_filter(df_MC, "fit_range")
#
#df_MC_corr = df_filter(df_MC_corr, "truth_matching_rd")
#df_MC_corr = df_filter(df_MC_corr, "trigger")
#df_MC_corr = df_filter(df_MC_corr, "preselection_corr")
#df_MC_corr = df_filter(df_MC_corr, "Jpsi")
#df_MC_corr = df_filter(df_MC_corr, "fit_range")

#plot_feature(
#    dfs = [df_MC_corr, df_MC_corr, df_sw],
#    names = ["L0MuonWeight", "L0DiMuonWeight"],
#    features=["L0MuonWeight_Run2_v1", "L0DiMuonWeight_Run2_v1"],
#    xlabel= "Weight",
#    ylabel= "Events",
#    xmin=0.85,
#    xmax=1.2,
#    nbins=50,
#    output_file="../plots/weights/L0_trigger_weights.pdf",
#    styles=["HIST", "HIST"],
#    markers=[21, 22],
#    colors=["blue", "red"],
#    legend_position=[0.5, 0.725, 0.7, 0.925],
#    xwidth=1100,
#    yheight=800,
#)

#plot_feature(
#    dfs = [df_MC_corr, df_MC, df_sw],
#    names = ["Corrected MC", "Uncorrected MC", "sWeighted Data"],
#    features=["L1_L0MuonDecision_TOS", "L1_L0MuonDecision_TOS", "L1_L0MuonDecision_TOS"],
#    xlabel= "L1 L0MuonDecision TOS",
#    ylabel= "Normalized Events",
#    xmin=-0.5,
#    xmax=1.5,
#    nbins=2,
#    output_file="../plots/weights/L1_L0MuonDecision_TOS.pdf",
#    styles=["E", "E", "E"],
#    markers=[20, 21, 22],
#    colors=["blue", "red", "black"],
#    weights=["L0MuonWeight_Run2_v1", None, "signal_yield_sw"],
#    normalize=True,
#    legend_position=[0.2, 0.6, 0.4, 0.8],
#)

#plot_feature(
#    dfs = [df_MC_corr, df_MC, df_sw],
#    names = ["Corrected MC", "Uncorrected MC", "sWeighted Data"],
#    features=["Jpsi_L0DiMuonDecision_TOS", "Jpsi_L0DiMuonDecision_TOS", "Jpsi_L0DiMuonDecision_TOS"],
#    xlabel= "Jpsi L0DiMuonDecision TOS",
#    ylabel= "Normalized Events",
#    xmin=-0.5,
#    xmax=1.5,
#    nbins=2,
#    output_file="../plots/weights/L1_L0DiMuonDecision_TOS.pdf",
#    styles=["E", "E", "E"],
#    markers=[20, 21, 22],
#    colors=["blue", "red", "black"],
#    weights=["L0DiMuonWeight_Run2_v1", None, "signal_yield_sw"],
#    normalize=True,
#    legend_position=[0.2, 0.65, 0.4, 0.85],
#)

# Tracktype comparison

#plot_feature(
#    dfs=[df_sw, df_MC],
#    names=["sWeighted Data", "MC"],
#    features=["Pion_TRACK_Type", "Pion_TRACK_Type"],
#    xlabel="Pion Track Type",
#    ylabel="Normalized Events",
#    xmin=0,
#    xmax=8,
#    nbins=8,
#    output_file="../plots/track_types/Pion_TRACK_TYPE.pdf",
#    styles=["E", "E"],
#    markers=[20, 21],
#    colors=["black", "red"],
#    weights=["signal_yield_sw", None],
#    normalize=True,
#    legend_position=[0.2, 0.7, 0.4, 0.9],
#    marker_colors=["black", "red"]
#)

#plot_feature(
#    dfs=[df_sw, df_MC],
#    names=["sWeighted Data", "MC"],
#    features=["Proton_TRACK_Type", "Proton_TRACK_Type"],
#    xlabel="Proton Track Type",
#    ylabel="Normalized Events",
#    xmin=0,
#    xmax=8,
#    nbins=8,
#    output_file="../plots/track_types/Proton_TRACK_TYPE.pdf",
#    styles=["E", "E"],
#    markers=[20, 21],
#    colors=["black", "red"],
#    weights=["signal_yield_sw", None],
#    normalize=True,
#    legend_position=[0.2, 0.7, 0.4, 0.9],
#    marker_colors=["black", "red"]
#)

# pid correction comparison

#plot_feature_and_pull(
#    dfs = [df_sw, df_MC, df_MC_corr],
#    names = ["sWeighted Data", "MC", "MC with PID correction"],
#    features = ["L1_ProbNNmu", "L1_ProbNNmu", "L1_ProbNNmu_pidcorr_default"],
#    xlabel = "ProbNNmu",
#    ylabel = "Normalized Events",
#    xmin = 0,
#    xmax = 1,
#    nbins = 50,
#    output_file = "../plots/pid_correction/ProbNNmu_comparison.pdf",
#    styles = ["E", "HIST", "HIST"],
#    markers = [20, 21, 22],
#    colors = ["black", "blue", "red"],
#    weights = ["signal_yield_sw", None, None],
#    normalize = True,
#    log_scale= True,
#    legend_position= [0.2, 0.7, 0.4, 0.9],
#    pull_index=[0, 2],
#    xwidth=800,
#    yheight=1000
#)

# Lifetime distibution comparison

#plot_feature_and_pull(
#    dfs = [df_MC, df_MC_corr, df_sw],
#    names = ["MC", "MC with lifetime correction", "sWeighted Data"],
#    features= ["Lambdab_TRUETAU", "Lambdab_TRUETAU", "Lambdab_TAU"],
#    xlabel = "Lambdab Lifetime [ns]",
#    ylabel = "Normalized Events",
#    xmin = 0,
#    xmax = 0.01,
#    nbins = 50,
#    output_file= "../plots/weights/Lambdab_TRUETAU_comparison.pdf",
#    styles = ["HIST", "HIST", "E"],
#    markers = [20, 21, 20],
#    colors = ["blue", "red", "black"],
#    weights = [None, "total_weights", "signal_yield_sw"],
#    normalize = True,
#    pull_index=[2, 1],
#    legend_position=[0.5, 0.7, 0.7, 0.9],
#    log_scale=True
#    
#)
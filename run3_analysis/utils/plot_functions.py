import ROOT as r
import yaml

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
r.gStyle.SetPadRightMargin(0.1)
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
    c.SetTopMargin(0.06) 
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
    histograms[sorted_indices[0]].GetYaxis().SetTitle(ylabel)
    histograms[sorted_indices[0]].GetYaxis().SetTitleSize(0.06)
    histograms[sorted_indices[0]].GetYaxis().SetLabelSize(0.06)
    histograms[sorted_indices[0]].GetYaxis().SetTitleOffset(1.2)

    
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
                 weights = None, xwidth = 800, yheight = 600, legend_position = [0.6, 0.7, 0.8, 0.9], log_scale = False, normalize = False,
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
    pad1.SetTopMargin(0.06)  # Set top margin to 0.02
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
    hist_pull.SetFillColor(r.kBlack)
    hist_pull.SetMarkerStyle(20)

    hist_pull.Draw("HIST")

    hist_pull.GetXaxis().SetTitle(xlabel)
    hist_pull.GetXaxis().SetLabelSize(padratio * 0.06)
    hist_pull.GetXaxis().SetTitleSize(padratio * 0.06)

    hist_pull.GetYaxis().SetTitle("Pull")
    hist_pull.GetYaxis().SetTitleOffset(1.2 / padratio)
    hist_pull.GetYaxis().SetLabelSize(padratio * 0.06)
    hist_pull.GetYaxis().SetTitleSize(padratio * 0.06)
    hist_pull.GetYaxis().SetNdivisions(505)

    
    # Save the plot
    c.SaveAs(output_file)

import plot_functions as pf
import ROOT as r
import yaml
import os
import sys
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__),".."))
sys.path.append(parent_dir) 
from filter import df_filter


# read the configuration files
with open("analysis_config.yml", "r") as analysis_config_file:
    config = yaml.safe_load(analysis_config_file)

# Create RDataFrames from the ROOT files
f_MC = [
    f"/ceph/users/jnicolini/scripts/Lambdab/Samples/angles/Lb2Lmm_PHMC_{year}_{magnet}.root" 
    for year in config["years"] 
    for magnet in config["magnet"]
    ]
df_MC = r.RDataFrame("DecayTree", f_MC)

f_MC_corr = [
    f"/ceph/users/pmachnik/Lambdab/MC/Lb2Lmm_PHMC_{year}_{mag}_corr_all.root"
    for year in config["years"]
    for mag in config["magnet"]
]
df_MC_corr = r.RDataFrame("DecayTree", f_MC_corr)

f_sw = "/ceph/users/pmachnik/Lambdab/Jpsi/Jpsi_sWeights.root"
df_sw = r.RDataFrame("DecayTree", f_sw)

df_MC = df_filter(df_MC, "truth_matching_rd")
df_MC = df_filter(df_MC, "trigger")
df_MC = df_filter(df_MC, "preselection")
df_MC = df_filter(df_MC, "Jpsi")
df_MC = df_filter(df_MC, "fit_range")

df_MC_corr = df_filter(df_MC_corr, "truth_matching_rd")
df_MC_corr = df_filter(df_MC_corr, "trigger")
df_MC_corr = df_filter(df_MC_corr, "preselection_corr")
df_MC_corr = df_filter(df_MC_corr, "Jpsi")
df_MC_corr = df_filter(df_MC_corr, "fit_range")

# plot
pf.plot_feature_and_pull(
    dfs = [df_sw, df_MC, df_MC_corr],
    names = ["sWeighted Data", "MC", "MC with PID correction"],
    features = ["L1_ProbNNmu", "L1_ProbNNmu", "L1_ProbNNmu_pidcorr_default"],
    xlabel = "ProbNNmu",
    ylabel = "Normalized Events",
    xmin = 0.1,
    xmax = 1,
    nbins = 50,
    output_file = "../plots/pid_correction/ProbNNmu_comparison.pdf",
    styles = ["E", "E", "E"],
    markers = [20, 20, 20],
    marker_colors= ["black", "blue", "red"],
    colors = ["black", "blue", "red"],
    weights = ["signal_yield_sw", None, None],
    normalize = True,
    log_scale= True,
    legend_position= [0.2, 0.7, 0.4, 0.9],
    pull_index=[0, 2]
)
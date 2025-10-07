import plot_functions as pf
import ROOT as r
import yaml
import os
import sys
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__),".."))
sys.path.append(parent_dir) 
from filter import df_filter
from define import BDT_define

# read the configuration files
with open("analysis_config.yml", "r") as analysis_config_file:
    config = yaml.safe_load(analysis_config_file)

# Create RDataFrames from the ROOT files
f_MC_corr = [
    f"/ceph/users/pmachnik/Lambdab/MC/Lb2Lmm_PHMC_{year}_{mag}_corr_all.root"
    for year in config["years"]
    for mag in config["magnet"]
]

f_data = [
     f"/ceph/users/jnicolini/scripts/Lambdab/Samples/angles/Lb2Lmm_Data_{year}_{magnet}.root"
     for year in config["years"]
     for magnet in config["magnet"]
]

df_MC_corr = r.RDataFrame("DecayTree", f_MC_corr)
df_data = r.RDataFrame("DecayTree", f_data)

# Apply filters
df_MC_corr = df_filter(df_MC_corr, "truth_matching_rd")
df_MC_corr = df_filter(df_MC_corr, "trigger")
df_MC_corr = df_filter(df_MC_corr, "preselection_corr")

df_data = df_filter(df_data, "trigger")
df_data = df_filter(df_data, "preselection")
df_data = df_filter(df_data, "background")

# Create BDT variables
df_MC_corr = BDT_define(df_MC_corr, corrected=True)
df_data = BDT_define(df_data)

# Apply high q2 filter
df_MC_corr = df_filter(df_MC_corr, "high_q2")
df_data = df_filter(df_data, "high_q2")

# plot the features
# Lambdab_TAU
pf.plot_feature(
    dfs=[df_MC_corr, df_data],
    normalize=True,
    styles=["HIST", "E"],
    colors=["red", "black"],
    weights=["total_weights", None],
    names=["Corrected MC", "Background"],
    markers=[20, 22],
    features=["Lambdab_TAU", "Lambdab_TAU"],
    xlabel="#tau_{#Lambda_{b}} [ns]",
    ylabel="Normalized Events",
    xmin=0,
    xmax=7e-3,
    nbins=50,
    output_file="../plots/BDT/Lambdab_TAU.pdf",
    legend_position=[0.6, 0.7, 0.9, 0.9],
)

# Lambdab_TAUCHI2
pf.plot_feature(
    dfs=[df_MC_corr, df_data],
    normalize=True,
    styles=["HIST", "E"],
    colors=["red", "black"],
    weights=["total_weights", None],
    names=["Corrected MC", "Background"],
    markers=[20, 22],
    features=["Lambdab_TAUCHI2", "Lambdab_TAUCHI2"],
    xlabel="#chi^{2}_{#tau, #Lambda_{b}}",
    ylabel="Normalized Events",
    xmin=0,
    xmax=35,
    nbins=50,
    output_file="../plots/BDT/Lambdab_TAUCHI2.pdf",
)

# Lambdab_ENDVERTEX_CHI2
pf.plot_feature(
    normalize=True,
    dfs=[df_MC_corr, df_data],
    styles=["HIST", "E"],
    colors=["red", "black"],
    weights=["total_weights", None],
    names=["Corrected MC", "Background"],
    markers=[20, 22],
    features=["Lambdab_ENDVERTEX_CHI2", "Lambdab_ENDVERTEX_CHI2"],
    xlabel="#chi^{2}_{Endvertex, #Lambda_{b}}",
    ylabel="Normalized Events",
    xmin=0,
    xmax=28,
    nbins=50,
    output_file="../plots/BDT/Lambdab_ENDVERTEX_CHI2.pdf",
)

pf.plot_feature(
    dfs=[df_MC_corr, df_data],
    normalize=True,
    styles=["HIST", "E"],
    colors=["red", "black"],
    weights=["total_weights", None],
    names=["Corrected MC", "Background"],
    markers=[20, 22],
    features=["log10_Lambdab_pT", "log10_Lambdab_pT"],
    xlabel="log10_Lambdab_pT",
    ylabel="Normalized Events",
    xmin=1,
    xmax=5,
    nbins=100,
    output_file="../plots/BDT/Lambdab_log10_pT.pdf",
    legend_position=[0.2, 0.7, 0.5, 0.9],
)

pf.plot_feature(
    dfs=[df_MC_corr, df_data],
    normalize=True,
    styles=["HIST", "E"],
    colors=["red", "black"],
    weights=["total_weights", None],
    names=["Corrected MC", "Background"],
    markers=[20, 22],
    features=["Lambdab_LOKI_DTF_CHI2NDOF", "Lambdab_LOKI_DTF_CHI2NDOF"],
    xlabel="Lambdab_LOKI_DTF_CHI2NDOF",
    ylabel="Normalized Events",
    xmin=0,
    xmax=15,
    nbins=100,
    output_file="../plots/BDT/Lambdab_LOKI_DTF_CHI2NDOF.pdf",
)

pf.plot_feature(
    dfs=[df_MC_corr, df_data],
    normalize=True,
    styles=["HIST", "E"],
    colors=["red", "black"],
    weights=["total_weights", None],
    names=["Corrected MC", "Background"],
    markers=[20, 22],
    features=["Lambdab_DIRA_OWNPV", "Lambdab_DIRA_OWNPV"],
    xlabel="Lambdab_DIRA_OWNPV",
    ylabel="Normalized Events",
    xmin=0.9995,
    xmax=1,
    nbins=100,
    output_file="../plots/BDT/Lambdab_DIRA_OWNPV.pdf",
    log_scale=True,
)

pf.plot_feature(
    dfs=[df_MC_corr, df_data],
    normalize=True,
    styles=["HIST", "E"],
    colors=["red", "black"],
    weights=["total_weights", None],
    names=["Corrected MC", "Background"],
    markers=[20, 22],
    features=["Lambdab_DiraCos", "Lambdab_DiraCos"],
    xlabel="Lambdab_DiraCos",
    ylabel="Normalized Events",
    xmin=0.9999,
    xmax=1.0000001,
    nbins=100,
    output_file="../plots/BDT/Lambdab_DiraCos.pdf",
    log_scale=True,
)

pf.plot_feature(
    dfs=[df_MC_corr, df_data],
    normalize=True,
    styles=["HIST", "E"],
    colors=["red", "black"],
    weights=["total_weights", None],
    names=["Corrected MC", "Background"],
    markers=[20, 22],
    features=["Lambdab_IPCHI2_OWNPV", "Lambdab_IPCHI2_OWNPV"],
    xlabel="Lambdab_IPCHI2_OWNPV",
    ylabel="Normalized Events",
    xmin=0,
    xmax=26,
    nbins=100,
    output_file="../plots/BDT/Lambdab_IPCHI2_OWNPV.pdf",
)

# Lambda features
pf.plot_feature(
    dfs=[df_MC_corr, df_data],
    normalize=True,
    styles=["HIST", "E"],
    colors=["red", "black"],
    weights=["total_weights", None],
    names=["Corrected MC", "Background"],
    markers=[20, 22],
    features=["Lambda_TAU", "Lambda_TAU"],
    xlabel="Lambda_TAU",
    ylabel="Normalized Events",
    xmin=0,
    xmax=0.5,
    nbins=100,
    output_file="../plots/BDT/Lambda_TAU.pdf",
)

pf.plot_feature(
    dfs=[df_MC_corr, df_data],
    normalize=True,
    styles=["HIST", "E"],
    colors=["red", "black"],
    weights=["total_weights", None],
    names=["Corrected MC", "Background"],
    markers=[20, 22],
    features=["Lambda_TAUCHI2", "Lambda_TAUCHI2"],
    xlabel="Lambda_TAUCHI2",
    ylabel="Normalized Events",
    xmin=0,
    xmax=30,
    nbins=100,
    output_file="../plots/BDT/Lambda_TAUCHI2.pdf",
)

pf.plot_feature(
    dfs=[df_MC_corr, df_data],
    normalize=True,
    styles=["HIST", "E"],
    colors=["red", "black"],
    weights=["total_weights", None],
    names=["Corrected MC", "Background"],
    markers=[20, 22],
    features=["Lambda_FDCHI2_OWNPV", "Lambda_FDCHI2_OWNPV"],
    xlabel="Lambda_FDCHI2_OWNPV",
    ylabel="Normalized Events",
    xmin=0,
    xmax=1700,
    nbins=100,
    output_file="../plots/BDT/Lambda_FDCHI2_OWNPV.pdf",
)

pf.plot_feature(
    dfs=[df_MC_corr, df_data],
    normalize=True,
    styles=["HIST", "E"],
    colors=["red", "black"],
    weights=["total_weights", None],
    names=["Corrected MC", "Background"],
    markers=[20, 22],
    features=["Lambda_ENDVERTEX_Z", "Lambda_ENDVERTEX_Z"],
    xlabel="Lambda_ENDVERTEX_Z",
    ylabel="Normalized Events",
    xmin=0,
    xmax=2400,
    nbins=100,
    output_file="../plots/BDT/Lambda_ENDVERTEX_Z.pdf",
)

pf.plot_feature(
    dfs=[df_MC_corr, df_data],
    normalize=True,
    styles=["HIST", "E"],
    colors=["red", "black"],
    weights=["total_weights", None],
    names=["Corrected MC", "Background"],
    markers=[20, 22],
    features=["Lambda_LOKI_DTF_CHI2NDOF", "Lambda_LOKI_DTF_CHI2NDOF"],
    xlabel="Lambda_LOKI_DTF_CHI2NDOF",
    ylabel="Normalized Events",
    xmin=0,
    xmax=20,
    nbins=100,
    output_file="../plots/BDT/Lambda_LOKI_DTF_CHI2NDOF.pdf",
)

pf.plot_feature(
    dfs=[df_MC_corr, df_data],
    normalize=True,
    styles=["HIST", "E"],
    colors=["red", "black"],
    weights=["total_weights", None],
    names=["Corrected MC", "Background"],
    markers=[20, 22],
    features=["log10_Lambda_pT", "log10_Lambda_pT"],
    xlabel="log10_Lambda_pT",
    ylabel="Normalized Events",
    xmin=1,
    xmax=5,
    nbins=100,
    output_file="../plots/BDT/Lambda_log10_pT.pdf",
)

pf.plot_feature(
    dfs=[df_MC_corr, df_data],
    normalize=True,
    styles=["HIST", "E"],
    colors=["red", "black"],
    weights=["total_weights", None],
    names=["Corrected MC", "Background"],
    markers=[20, 22],
    features=["Lambda_IPCHI2_ORIVX", "Lambda_IPCHI2_ORIVX"],
    xlabel="Lambda_IPCHI2_ORIVX",
    ylabel="Normalized Events",
    xmin=0,
    xmax=30,
    nbins=100,
    output_file="../plots/BDT/Lambda_IPCHI2_ORIVX.pdf",
)

# pion features
pf.plot_feature(
    dfs=[df_MC_corr, df_data],
    normalize=True,
    styles=["HIST", "E"],
    colors=["red", "black"],
    weights=["total_weights", None],
    names=["Corrected MC", "Background"],
    markers=[20, 22],
    features=["Pion_PT", "Pion_PT"],
    xlabel="Pion_PT",
    ylabel="Normalized Events",
    xmin=0,
    xmax=1200,
    nbins=100,
    output_file="../plots/BDT/Pion_PT.pdf",
)

pf.plot_feature(
    dfs=[df_MC_corr, df_data],
    normalize=True,
    styles=["HIST", "E"],
    colors=["red", "black"],
    weights=["total_weights", None],
    names=["Corrected MC", "Background"],
    markers=[20, 22],
    features=["Pion_IPCHI2_ORIVX", "Pion_IPCHI2_ORIVX"],
    xlabel="Pion_IPCHI2_ORIVX",
    ylabel="Normalized Events",
    xmin=0,
    xmax=3,
    nbins=100,
    output_file="../plots/BDT/Pion_IPCHI2_ORIVX.pdf",
)

pf.plot_feature(
    dfs=[df_MC_corr, df_data],
    normalize=True,
    styles=["HIST", "E"],
    colors=["red", "black"],
    weights=["total_weights", None],
    names=["Corrected MC", "Background"],
    markers=[20, 22],
    features=["Pion_ETA", "Pion_ETA"],
    xlabel="Pion_ETA",
    ylabel="Normalized Events",
    xmin=1.5,
    xmax=6,
    nbins=100,
    output_file="../plots/BDT/Pion_ETA.pdf",
)

# proton features
pf.plot_feature(
    dfs=[df_MC_corr, df_data],
    normalize=True,
    styles=["HIST", "E"],
    colors=["red", "black"],
    weights=["total_weights", None],
    names=["Corrected MC", "Background"],
    markers=[20, 22],
    features=["Proton_PT", "Proton_PT"],
    xlabel="Proton_PT",
    ylabel="Normalized Events",
    xmin=0,
    xmax=6000,
    nbins=100,
    output_file="../plots/BDT/Proton_PT.pdf",
)

pf.plot_feature(
    dfs=[df_MC_corr, df_data],
    normalize=True,
    styles=["HIST", "E"],
    colors=["red", "black"],
    weights=["total_weights", None],
    names=["Corrected MC", "Background"],
    markers=[20, 22],
    features=["log10_Proton_pT", "log10_Proton_pT"],
    xlabel="log10_Proton_pT",
    ylabel="Normalized Events",
    xmin=2,
    xmax=4.5,
    nbins=100,
    output_file="../plots/BDT/Proton_log10_pT.pdf",
)

pf.plot_feature(
    dfs=[df_MC_corr, df_data],
    normalize=True,
    styles=["HIST", "E"],
    colors=["red", "black"],
    weights=["total_weights", None],
    names=["Corrected MC", "Background"],
    markers=[20, 22],
    features=["Proton_IPCHI2_ORIVX", "Proton_IPCHI2_ORIVX"],
    xlabel="Proton_IPCHI2_ORIVX",
    ylabel="Normalized Events",
    xmin=5,
    xmax=40,
    nbins=100,
    output_file="../plots/BDT/Proton_IPCHI2_ORIVX.pdf",
)

# muon features
pf.plot_feature(
    dfs=[df_MC_corr, df_data],
    normalize=True,
    styles=["HIST", "E"],
    colors=["red", "black"],
    weights=["total_weights", None],
    names=["Corrected MC", "Background"],
    markers=[20, 22],
    features=["L1_PT", "L1_PT"],
    xlabel="L1_PT",
    ylabel="Normalized Events",
    xmin=800,
    xmax=7000,
    nbins=100,
    output_file="../plots/BDT/L1_PT.pdf",
)

pf.plot_feature(
    dfs=[df_MC_corr, df_data],
    normalize=True,
    styles=["HIST", "E"],
    colors=["red", "black"],
    weights=["total_weights", None],
    names=["Corrected MC", "Background"],
    markers=[20, 22],
    features=["L2_PT", "L2_PT"],
    xlabel="L2_PT",
    ylabel="Normalized Events",
    xmin=800,
    xmax=7000,
    nbins=100,
    output_file="../plots/BDT/L2_PT.pdf",
)

pf.plot_feature(
    dfs=[df_MC_corr, df_data],
    normalize=True,
    styles=["HIST", "E"],
    colors=["red", "black"],
    weights=["total_weights", None],
    names=["Corrected MC", "Background"],
    markers=[20, 22],
    features=["ProbNNmu_max", "ProbNNmu_max"],
    xlabel="ProbNNmu_max",
    ylabel="Normalized Events",
    xmin=0.1,
    xmax=1,
    nbins=100,
    output_file="../plots/BDT/ProbNNmu_max.pdf",
    log_scale=True,
)

pf.plot_feature(
    dfs=[df_MC_corr, df_data],
    normalize=True,
    styles=["HIST", "E"],
    colors=["red", "black"],
    weights=["total_weights", None],
    names=["Corrected MC", "Background"],
    markers=[20, 22],
    features=["ProbNNmu_min", "ProbNNmu_min"],
    xlabel="ProbNNmu_min",
    ylabel="Normalized Events",
    xmin=0.1,
    xmax=1,
    nbins=100,
    output_file="../plots/BDT/ProbNNmu_min.pdf",
    log_scale=True,
)

pf.plot_feature(
    dfs=[df_MC_corr, df_data],
    normalize=True,
    styles=["HIST", "E"],
    colors=["red", "black"],
    weights=["total_weights", None],
    names=["Corrected MC", "Background"],
    markers=[20, 22],
    features=["chi2ip_max", "chi2ip_max"],
    xlabel="chi2ip_max",
    ylabel="Normalized Events",
    xmin=0,
    xmax=8,
    nbins=100,
    output_file="../plots/BDT/chi2ip_max.pdf",
)

pf.plot_feature(
    dfs=[df_MC_corr, df_data],
    normalize=True,
    styles=["HIST", "E"],
    colors=["red", "black"],
    weights=["total_weights", None],
    names=["Corrected MC", "Background"],
    markers=[20, 22],
    features=["chi2ip_min", "chi2ip_min"],
    xlabel="chi2ip_min",
    ylabel="Normalized Events",
    xmin=0,
    xmax=3,
    nbins=100,
    output_file="../plots/BDT/chi2ip_min.pdf",
)

# mixed features
pf.plot_feature(
    dfs=[df_MC_corr, df_data],
    normalize=True,
    styles=["HIST", "E"],
    colors=["red", "black"],
    weights=["total_weights", None],
    names=["Corrected MC", "Background"],
    markers=[20, 22],
    features=["chi2ip_min_pi_p", "chi2ip_min_pi_p"],
    xlabel="chi2ip_min_pi_p",
    ylabel="Normalized Events",
    xmin=9,
    xmax=50,
    nbins=100,
    output_file="../plots/BDT/chi2ip_min_pi_p.pdf",
)

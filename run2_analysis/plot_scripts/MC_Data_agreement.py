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
    f"/ceph/users/pmachnik/Lambdab/MC/Lb2Lmm_JpsiMC_{year}_{mag}_corr_all.root"
    for year in config["years"]
    for mag in config["magnet"]
]

f_MC = [
    f"/ceph/users/jnicolini/scripts/Lambdab/Samples/angles/Lb2Lmm_JpsiMC_{year}_{mag}.root"
    for year in config["years"]
    for mag in config["magnet"]
]

f_sw_ll = f"/ceph/users/pmachnik/Lambdab/Jpsi/Jpsi_sWeights_ll.root"
f_sw_dd = f"/ceph/users/pmachnik/Lambdab/Jpsi/Jpsi_sWeights_dd.root"

df_sw_ll = r.RDataFrame("DecayTree", f_sw_ll)
df_sw_dd = r.RDataFrame("DecayTree", f_sw_dd)
df_MC_corr = r.RDataFrame("DecayTree", f_MC_corr)
df_MC = r.RDataFrame("DecayTree", f_MC)

# Apply filters
df_MC_corr = df_filter(df_MC_corr, "truth_matching_Jpsi")
df_MC_corr = df_filter(df_MC_corr, "trigger")
df_MC_corr = df_filter(df_MC_corr, "preselection")
df_MC_corr = df_filter(df_MC_corr, "fit_range")
df_MC_corr = df_filter(df_MC_corr, "Jpsi")

df_MC = df_filter(df_MC, "truth_matching_Jpsi")
df_MC = df_filter(df_MC, "trigger")
df_MC = df_filter(df_MC, "preselection")
df_MC = df_filter(df_MC, "fit_range")
df_MC = df_filter(df_MC, "Jpsi")

# Apply BDT definitions
df_MC_corr = BDT_define(df_MC_corr)
df_MC = BDT_define(df_MC)

# create ll and dd MC dataframes
df_MC_ll = df_filter(df_MC, "ll")
df_MC_dd = df_filter(df_MC, "dd")
df_MC_corr_ll = df_filter(df_MC_corr, "ll")
df_MC_corr_dd = df_filter(df_MC_corr, "dd")

dfs_MC = [df_MC_ll, df_MC_dd]
dfs_MC_corr = [df_MC_corr_ll, df_MC_corr_dd]
dfs_sw = [df_sw_ll, df_sw_dd]
tt = ["ll", "dd"]


xmin_dict = {
    "Lambdab_TAU": {"ll": 1e-4, "dd": 1e-4},
    "log10_Lambdab_pT": {"ll": 1.5, "dd": 1.5},
    "Lambdab_LOKI_DTF_CHI2NDOF": {"ll": 0, "dd": 0},
    "Lambdab_ENDVERTEX_CHI2": {"ll": 0, "dd": 0},
    "Lambdab_TAUCHI2": {"ll": 0, "dd": 0},
    "log10_Lambdab_DIRA_OWNPV": {"ll": -2.2e-4, "dd": -2.2e-4},
    "Lambda_FDCHI2_OWNPV": {"ll": 0, "dd": 0},
    "Lambda_FD_TOPPV": {"ll": 0, "dd": 200},
    "Lambda_IPCHI2_OWNPV": {"ll": 0, "dd": 0},
    "Lambda_TAU": {"ll": 0, "dd": 0},
    "Lambda_TAUCHI2": {"ll": 0, "dd": 0},
    "Pion_IPCHI2_OWNPV": {"ll": 0, "dd": 0},
    "log10_Proton_pT": {"ll": 1.8, "dd": 1.8},
    "log10_L1_pT": {"ll": 2.8, "dd": 2.8},
    "ProbNNmu_max": {"ll": 0.1, "dd": 0.1},
    "ProbNNmu_min": {"ll": 0.1, "dd": 0.1},
    "log10_chi2ip_max": {"ll": 0.8, "dd": 0.8},
    "log10_chi2ip_min": {"ll": 0.8, "dd": 0.8},
    "log10_chi2ip_min_pi_p": {"ll": 0.8, "dd": 0.8},
    "L1_P": {"ll": 3, "dd": 3},
    "Lambda_M": {"ll": 1108, "dd": 1108},
    "Lambda_ENDVERTEX_Z": {"ll": 0, "dd": 200},
    "ProbNNmu": {"ll": 0, "dd": 0}, 
}
xmax_dict = {
    "Lambdab_TAU": {"ll": 8e-3, "dd": 8e-3},
    "log10_Lambdab_pT": {"ll": 5, "dd": 5},
    "Lambdab_LOKI_DTF_CHI2NDOF": {"ll": 10, "dd": 10},
    "Lambdab_ENDVERTEX_CHI2": {"ll": 20, "dd": 20},
    "Lambdab_TAUCHI2": {"ll": 20, "dd": 20},
    "log10_Lambdab_DIRA_OWNPV": {"ll": 0, "dd": 0},
    "Lambda_FDCHI2_OWNPV": {"ll": 3000, "dd": 3000},
    "Lambda_FD_TOPPV": {"ll": 750, "dd": 2650},
    "Lambda_IPCHI2_OWNPV": {"ll": 200, "dd": 25},
    "Lambda_TAU": {"ll": 2.1e-1, "dd": 6e-1},
    "Lambda_TAUCHI2": {"ll": 10, "dd": 15},
    "Pion_IPCHI2_OWNPV": {"ll": 3000, "dd": 1500},
    "log10_Proton_pT": {"ll": 4.2, "dd": 4.5},
    "log10_L1_pT": {"ll": 4.8, "dd": 4.8},
    "ProbNNmu_max": {"ll": 1, "dd": 1},
    "ProbNNmu_min": {"ll": 1, "dd": 1},
    "log10_chi2ip_max": {"ll": 5, "dd": 5},
    "log10_chi2ip_min": {"ll": 5, "dd": 4.8},
    "log10_chi2ip_min_pi_p": {"ll": 4, "dd": 3},
    "L1_P": {"ll": 120, "dd": 120},
    "Lambda_M": {"ll": 1124, "dd": 1124},
    "Lambda_ENDVERTEX_Z": {"ll": 700, "dd": 2400},
    "ProbNNmu": {"ll": 1, "dd": 1},
}
# plot the features
# Lambdab_TAU
for df_MC, df_MC_corr, df_sw, tracktype in zip(dfs_MC, dfs_MC_corr, dfs_sw, tt):
    feature = "Lambdab_TAU"
    pf.plot_feature(
        dfs=[df_MC, df_MC_corr, df_sw],
        normalize=True,
        styles=["E", "E", "E"],
        colors=["blue", "red", "black"],
        weights=[None, "total_weights", "signal_yield_sw"],
        names=["MC", "Corrected MC", "sWeighted Data"],
        markers=[20, 20, 20],
        marker_colors=["blue", "red", "black"],
        features=[feature, feature, feature],
        xlabel="#tau(#Lambda^{0}_{b}) [ns]",
        ylabel="Normalized Events",
        xmin=xmin_dict[feature][tracktype],
        xmax=xmax_dict[feature][tracktype],
        nbins=50,
        output_file=f"../plots/MC_Data_agreement/R2/{tracktype}/Lambdab_TAU.pdf",
        legend_position=[0.6, 0.7, 0.9, 0.9],
    )

    # Lambdab_TAUCHI2
    feature = "Lambda_TAUCHI2"
    pf.plot_feature(
        dfs=[df_MC, df_MC_corr, df_sw],
        normalize=True,
        styles=["E", "E", "E"],
        colors=["blue", "red", "black"],
        weights=[None, "total_weights", "signal_yield_sw"],
        names=["MC", "Corrected MC", "sWeighted Data"],
        markers=[20, 20, 20],
        marker_colors= ["blue", "red", "black"],
        features=["Lambdab_TAUCHI2", "Lambdab_TAUCHI2", "Lambdab_TAUCHI2"],
        xlabel="#chi^{2}_{#tau}(#Lambda^{0}_{b})",
        ylabel="Normalized Events",
        xmin= xmin_dict[feature][tracktype],
        xmax= xmax_dict[feature][tracktype],
        nbins=50,
        output_file=f"../plots/MC_Data_agreement/R2/{tracktype}/Lambdab_TAUCHI2.pdf",
    )
    
    # Lambdab_ENDVERTEX_CHI2
    feature = "Lambdab_ENDVERTEX_CHI2"
    pf.plot_feature(
        normalize=True,
        dfs=[df_MC, df_MC_corr, df_sw],
        styles=["E", "E", "E"],
        colors=["blue", "red", "black"],
        weights=[None, "total_weights", "signal_yield_sw"],
        names=["MC", "Corrected MC", "sWeighted Data"],
        markers=[20, 20, 20],
        marker_colors= ["blue", "red", "black"],
        features=["Lambdab_ENDVERTEX_CHI2", "Lambdab_ENDVERTEX_CHI2", "Lambdab_ENDVERTEX_CHI2"],
        xlabel="#chi^{2}_{endvertex}(#Lambda^{0}_{b})",
        ylabel="Normalized Events",
        xmin=xmin_dict[feature][tracktype],
        xmax=xmax_dict[feature][tracktype],
        nbins=50,
        output_file=f"../plots/MC_Data_agreement/R2/{tracktype}/Lambdab_ENDVERTEX_CHI2.pdf",
    )

    feature = "log10_Lambdab_pT"
    pf.plot_feature(
        dfs=[df_MC, df_MC_corr, df_sw],
        normalize=True,
        styles=["E", "E", "E"],
        colors=["blue", "red", "black"],
        weights=[None, "total_weights", "signal_yield_sw"],
        names=["MC", "Corrected MC", "sWeighted Data"],
        markers=[20, 20, 20],
        marker_colors=["blue", "red", "black"],
        features=[feature, feature, feature],
        xlabel="log_{10}(p_{T}(#Lambda^{0}_{b}))",
        ylabel="Normalized Events",
        xmin=xmin_dict[feature][tracktype],
        xmax=xmax_dict[feature][tracktype],
        nbins=50,
        output_file=f"../plots/MC_Data_agreement/R2/{tracktype}/Lambdab_log10_pT.pdf",
        legend_position=[0.2, 0.7, 0.5, 0.9],
    )

    feature = "Lambdab_LOKI_DTF_CHI2NDOF"
    pf.plot_feature(
        dfs=[df_MC, df_MC_corr, df_sw],
        normalize=True,
        styles=["E", "E", "E"],
        colors=["blue", "red", "black"],
        weights=[None, "total_weights", "signal_yield_sw"],
        names=["MC", "Corrected MC", "sWeighted Data"],
        markers=[20, 20, 20],
        marker_colors=["blue", "red", "black"],
        features=[feature, feature, feature],
        xlabel="#chi^{2}_{m}(#Lambda^{0}_{b})",
        ylabel="Normalized Events",
        xmin=xmin_dict[feature][tracktype],
        xmax=xmax_dict[feature][tracktype],
        nbins=50,
        output_file=f"../plots/MC_Data_agreement/R2/{tracktype}/Lambdab_LOKI_DTF_CHI2NDOF.pdf",
    )

    feature = "log10_Lambdab_DIRA_OWNPV"
    pf.plot_feature(
        dfs=[df_MC, df_MC_corr, df_sw],
        normalize=True,
        styles=["E", "E", "E"],
        colors=["blue", "red", "black"],
        weights=[None, "total_weights", "signal_yield_sw"],
        names=["MC", "Corrected MC", "sWeighted Data"],
        markers=[20, 20, 20],
        marker_colors=["blue", "red", "black"],
        features=[feature, feature, feature],
        xlabel="log_{10}(#Theta_{DIRA}(#Lambda^{0}_{b}))",
        ylabel="Normalized Events",
        xmin=xmin_dict[feature][tracktype],
        xmax=xmax_dict[feature][tracktype],
        nbins=50,
        output_file=f"../plots/MC_Data_agreement/R2/{tracktype}/log10_Lambdab_DIRA_OWNPV.pdf",
        log_scale=True,
        legend_position=[0.2, 0.7, 0.5, 0.9],
    )

    #pf.plot_feature(
    #    dfs=[df_MC, df_MC_corr, df_sw],
    #    normalize=True,
    #    styles=["E", "E", "E"],
    #    colors=["blue", "red", "black"],
    #    weights=[None, "total_weights", "signal_yield_sw"],
    #    names=["MC", "Corrected MC", "sWeighted Data"],
    #    markers=[20, 01, 02],
    #    features=["Lambdab_DiraCos", "Lambdab_DiraCos", "Lambdab_DiraCos"],
    #    xlabel="Lambdab_DiraCos",
    #    ylabel="Normalized Events",
    #    xmin=0.9999,
    #    xmax=1.0000001,
    #    nbins=50,
    #    output_file=f"../plots/MC_Data_agreement/R2/{tracktype}/Lambdab_DiraCos.pdf",
    #    log_scale=True,
    #)

    #pf.plot_feature(
    #    dfs=[df_MC, df_MC_corr, df_sw],
    #    normalize=True,
    #    styles=["E", "E", "E"],
    #    colors=["blue", "red", "black"],
    #    weights=[None, "total_weights", "signal_yield_sw"],
    #    names=["MC", "Corrected MC", "sWeighted Data"],
    #    markers=[20, 01, 02],
    #    features=["Lambdab_IPCHI2_OWNPV", "Lambdab_IPCHI2_OWNPV", "Lambdab_IPCHI2_OWNPV"],
    #    xlabel="Lambdab_IPCHI2_OWNPV",
    #    ylabel="Normalized Events",
    #    xmin=0,
    #    xmax=36,
    #    nbins=50,
    #    output_file=f"../plots/MC_Data_agreement/R2/{tracktype}/Lambdab_IPCHI2_OWNPV.pdf",
    #)

    # Lambda features
    feature = "Lambda_TAU"
    pf.plot_feature(
        dfs=[df_MC, df_MC_corr, df_sw],
        normalize=True,
        styles=["E", "E", "E"],
        colors=["blue", "red", "black"],
        weights=[None, "total_weights", "signal_yield_sw"],
        names=["MC", "Corrected MC", "sWeighted Data"],
        markers=[20, 20, 20],
        marker_colors=["blue", "red", "black"],
        features=[feature, feature, feature],
        xlabel="#tau(#Lambda^{0})",
        ylabel="Normalized Events",
        xmin=xmin_dict[feature][tracktype],
        xmax=xmax_dict[feature][tracktype],
        nbins=50,
        output_file=f"../plots/MC_Data_agreement/R2/{tracktype}/Lambda_TAU.pdf",
    )

    # Lambda_FD_TOPPV
    feature = "Lambda_FD_TOPPV"
    legend_position = [0.3, 0.3, 0.6, 0.5] if tracktype == "dd" else [0.6, 0.7, 0.9, 0.9]
    pf.plot_feature(
        dfs=[df_MC, df_MC_corr, df_sw],
        normalize=True,
        styles=["E", "E", "E"],
        colors=["blue", "red", "black"],
        weights=[None, "total_weights", "signal_yield_sw"],
        names=["MC", "Corrected MC", "sWeighted Data"],
        markers=[20, 20, 20],
        marker_colors=["blue", "red", "black"],
        features=[feature, feature, feature],
        xlabel="FD(#Lambda^{0}) [mm]",
        ylabel="Normalized Events",
        xmin=xmin_dict[feature][tracktype],
        xmax=xmax_dict[feature][tracktype],
        nbins=50,
        output_file=f"../plots/MC_Data_agreement/R2/{tracktype}/Lambda_FD_TOPPV.pdf",
        legend_position=legend_position,
    )

    # Lambda_TAUCHI2
    feature = "Lambda_TAUCHI2"
    pf.plot_feature(
        dfs=[df_MC, df_MC_corr, df_sw],
        normalize=True,
        styles=["E", "E", "E"],
        colors=["blue", "red", "black"],
        weights=[None, "total_weights", "signal_yield_sw"],
        names=["MC", "Corrected MC", "sWeighted Data"],
        markers=[20, 20, 20],
        marker_colors= ["blue", "red", "black"],
        features=["Lambda_TAUCHI2", "Lambda_TAUCHI2", "Lambda_TAUCHI2"],
        xlabel="#chi^{2}_{#tau}(#Lambda^{0})",
        ylabel="Normalized Events",
        xmin= xmin_dict[feature][tracktype],
        xmax= xmax_dict[feature][tracktype],
        nbins=50,
        output_file=f"../plots/MC_Data_agreement/R2/{tracktype}/Lambda_TAUCHI2.pdf",
    )

    #pf.plot_feature(
    #    dfs=[df_MC, df_MC_corr, df_sw],
    #    normalize=True,
    #    styles=["E", "E", "E"],
    #    colors=["blue", "red", "black"],
    #    weights=[None, "total_weights", "signal_yield_sw"],
    #    names=["MC", "Corrected MC", "sWeighted Data"],
    #    markers=[20, 01, 02],
    #    features=["log10_Lambda_pT", "log10_Lambda_pT", "log10_Lambda_pT"],
    #    xlabel="log10_Lambda_pT",
    #    ylabel="Normalized Events",
    #    xmin=1,
    #    xmax=5,
    #    nbins=50,
    #    output_file=f"../plots/MC_Data_agreement/R2/{tracktype}/Lambda_log10_pT.pdf",
    #)

    #pf.plot_feature(
    #    dfs=[df_MC, df_MC_corr, df_sw],
    #    normalize=True,
    #    styles=["E", "E", "E"],
    #    colors=["blue", "red", "black"],
    #    weights=[None, "total_weights", "signal_yield_sw"],
    #    names=["MC", "Corrected MC", "sWeighted Data"],
    #    markers=[20, 01, 02],
    #    features=["Lambda_FDCHI2_OWNPV", "Lambda_FDCHI2_OWNPV", "Lambda_FDCHI2_OWNPV"],
    #    xlabel="Lambda_FDCHI2_OWNPV",
    #    ylabel="Normalized Events",
    #    xmin=0,
    #    xmax=1300,
    #    nbins=50,
    #    output_file=f"../plots/MC_Data_agreement/R2/{tracktype}/Lambda_FDCHI2_OWNPV.pdf",
    #)

    #pf.plot_feature(
    #    dfs=[df_MC, df_MC_corr, df_sw],
    #    normalize=True,
    #    styles=["E", "E", "E"],
    #    colors=["blue", "red", "black"],
    #    weights=[None, "total_weights", "signal_yield_sw"],
    #    names=["MC", "Corrected MC", "sWeighted Data"],
    #    markers=[20, 01, 02],
    #    features=["Lambda_ENDVERTEX_Z", "Lambda_ENDVERTEX_Z", "Lambda_ENDVERTEX_Z"],
    #    xlabel="Lambda_ENDVERTEX_Z",
    #    ylabel="Normalized Events",
    #    xmin=0,
    #    xmax=2500,
    #    nbins=50,
    #    output_file=f"../plots/MC_Data_agreement/R2/{tracktype}/Lambda_ENDVERTEX_Z.pdf",
    #)

    #pf.plot_feature(
    #    dfs=[df_MC, df_MC_corr, df_sw],
    #    normalize=True,
    #    styles=["E", "E", "E"],
    #    colors=["blue", "red", "black"],
    #    weights=[None, "total_weights", "signal_yield_sw"],
    #    names=["MC", "Corrected MC", "sWeighted Data"],
    #    markers=[20, 01, 02],
    #    features=["Lambda_LOKI_DTF_CHI2NDOF", "Lambda_LOKI_DTF_CHI2NDOF", "Lambda_LOKI_DTF_CHI2NDOF"],
    #    xlabel="Lambda_LOKI_DTF_CHI2NDOF",
    #    ylabel="Normalized Events",
    #    xmin=0,
    #    xmax=20,
    #    nbins=50,
    #    output_file=f"../plots/MC_Data_agreement/R2/{tracktype}/Lambda_LOKI_DTF_CHI2NDOF.pdf",
    #)

    feature = "Lambda_IPCHI2_OWNPV"
    pf.plot_feature(
        dfs=[df_MC, df_MC_corr, df_sw],
        normalize=True,
        styles=["E", "E", "E"],
        colors=["blue", "red", "black"],
        weights=[None, "total_weights", "signal_yield_sw"],
        names=["MC", "Corrected MC", "sWeighted Data"],
        markers=[20, 20, 20],
        marker_colors=["blue", "red", "black"],
        features=[feature, feature, feature],
        xlabel="#chi^{2}_{IP}(#Lambda^{0})",
        ylabel="Normalized Events",
        xmin=xmin_dict[feature][tracktype],
        xmax=xmax_dict[feature][tracktype],
        nbins=50,
        output_file=f"../plots/MC_Data_agreement/R2/{tracktype}/Lambda_IPCHI2_OWNPV.pdf",
    )

    # pion features
    #pf.plot_feature(
    #    dfs=[df_MC, df_MC_corr, df_sw],
    #    normalize=True,
    #    styles=["E", "E", "E"],
    #    colors=["blue", "red", "black"],
    #    weights=[None, "total_weights", "signal_yield_sw"],
    #    names=["MC", "Corrected MC", "sWeighted Data"],
    #    markers=[20, 01, 02],
    #    features=["Pion_PT", "Pion_PT", "Pion_PT"],
    #    xlabel="Pion_PT",
    #    ylabel="Normalized Events",
    #    xmin=0,
    #    xmax=1200,
    #    nbins=50,
    #    output_file=f"../plots/MC_Data_agreement/R2/{tracktype}/Pion_PT.pdf",
    #)

    feature = "Pion_IPCHI2_OWNPV"
    pf.plot_feature(
        dfs=[df_MC, df_MC_corr, df_sw],
        normalize=True,
        styles=["E", "E", "E"],
        colors=["blue", "red", "black"],
        weights=[None, "total_weights", "signal_yield_sw"],
        names=["MC", "Corrected MC", "sWeighted Data"],
        markers=[20, 20, 20],
        marker_colors=["blue", "red", "black"],
        features=[feature, feature, feature],
        xlabel="#chi^{2}_{IP}(#pi)",
        ylabel="Normalized Events",
        xmin=xmin_dict[feature][tracktype],
        xmax=xmax_dict[feature][tracktype],
        nbins=50,
        output_file=f"../plots/MC_Data_agreement/R2/{tracktype}/Pion_IPCHI2_OWNPV.pdf",
    )

    #pf.plot_feature(
    #    dfs=[df_MC, df_MC_corr, df_sw],
    #    normalize=True,
    #    styles=["E", "E", "E"],
    #    colors=["blue", "red", "black"],
    #    weights=[None, "total_weights", "signal_yield_sw"],
    #    names=["MC", "Corrected MC", "sWeighted Data"],
    #    markers=[20, 01, 02],
    #    features=["Pion_ETA", "Pion_ETA", "Pion_ETA"],
    #    xlabel="Pion_ETA",
    #    ylabel="Normalized Events",
    #    xmin=1.5,
    #    xmax=6,
    #    nbins=50,
    #    output_file=f"../plots/MC_Data_agreement/R2/{tracktype}/Pion_ETA.pdf",
    #)

    # proton features
    #pf.plot_feature(
    #    dfs=[df_MC, df_MC_corr, df_sw],
    #    normalize=True,
    #    styles=["E", "E", "E"],
    #    colors=["blue", "red", "black"],
    #    weights=[None, "total_weights", "signal_yield_sw"],
    #    names=["MC", "Corrected MC", "sWeighted Data"],
    #    markers=[20, 01, 02],
    #    features=["log10_Proton_PT", "log10_Proton_PT", "log10_Proton_PT"],
    #    xlabel="Proton_PT",
    #    ylabel="Normalized Events",
    #    xmin=0,
    #    xmax=4,
    #    nbins=50,
    #    output_file=f"../plots/MC_Data_agreement/R2/{tracktype}/log10_Proton_PT.pdf",
    #)

    # Proton ETA
    #pf.plot_feature(
    #    dfs=[df_MC, df_MC_corr, df_sw],
    #    normalize=True,
    #    styles=["E", "E", "E"],
    #    colors=["blue", "red", "black"],
    #    weights=[None, "total_weights", "signal_yield_sw"],
    #    names=["MC", "Corrected MC", "sWeighted Data"],
    #    markers=[20, 01, 02],
    #    features=["Proton_ETA", "Proton_ETA", "Proton_ETA"],
    #    xlabel="Proton_ETA",
    #    ylabel="Normalized Events",
    #    xmin=1.9,
    #    xmax=5.5,
    #    nbins=50,
    #    output_file=f"../plots/MC_Data_agreement/R2/{tracktype}/Proton_ETA.pdf",
    #)


    feature = "log10_Proton_pT"
    pf.plot_feature(
        dfs=[df_MC, df_MC_corr, df_sw],
        normalize=True,
        styles=["E", "E", "E"],
        colors=["blue", "red", "black"],
        weights=[None, "total_weights", "signal_yield_sw"],
        names=["MC", "Corrected MC", "sWeighted Data"],
        markers=[20, 20, 20],
        marker_colors=["blue", "red", "black"],
        features=[feature, feature, feature],
        xlabel="log_{10}(p_{T}(p))",
        ylabel="Normalized Events",
        xmin=xmin_dict[feature][tracktype],
        xmax=xmax_dict[feature][tracktype],
        nbins=50,
        output_file=f"../plots/MC_Data_agreement/R2/{tracktype}/log10_Proton_pT.pdf",
        legend_position=[0.2, 0.7, 0.5, 0.9],
    )

    #pf.plot_feature(
    #    dfs=[df_MC, df_MC_corr, df_sw],
    #    normalize=True,
    #    styles=["E", "E", "E"],
    #    colors=["blue", "red", "black"],
    #    weights=[None, "total_weights", "signal_yield_sw"],
    #    names=["MC", "Corrected MC", "sWeighted Data"],
    #    markers=[20, 01, 02],
    #    features=["Proton_IPCHI2_OWNPV", "Proton_IPCHI2_OWNPV", "Proton_IPCHI2_OWNPV"],
    #    xlabel="Proton_IPCHI2_OWNPV",
    #    ylabel="Normalized Events",
    #    xmin=8,
    #    xmax=300,
    #    nbins=50,
    #    output_file=f"../plots/MC_Data_agreement/R2/{tracktype}/Proton_IPCHI2_OWNPV.pdf",
    #)

    # muon features
    feature = "log10_L1_pT"
    pf.plot_feature(
        dfs=[df_MC, df_MC_corr, df_sw],
        normalize=True,
        styles=["E", "E", "E"],
        colors=["blue", "red", "black"],
        weights=[None, "total_weights", "signal_yield_sw"],
        names=["MC", "Corrected MC", "sWeighted Data"],
        markers=[20, 20, 20],
        marker_colors=["blue", "red", "black"],
        features=[feature, feature, feature],
        xlabel="log_{10}(p_{T}(#mu))",
        ylabel="Normalized Events",
        xmin=xmin_dict[feature][tracktype],
        xmax=xmax_dict[feature][tracktype],
        nbins=50,
        output_file=f"../plots/MC_Data_agreement/R2/{tracktype}/log10_L1_pT.pdf",
    )

    #pf.plot_feature(
    #    dfs=[df_MC, df_MC_corr, df_sw],
    #    normalize=True,
    #    styles=["E", "E", "E"],
    #    colors=["blue", "red", "black"],
    #    weights=[None, "total_weights", "signal_yield_sw"],
    #    names=["MC", "Corrected MC", "sWeighted Data"],
    #    markers=[20, 01, 02],
    #    features=["L2_PT", "L2_PT", "L2_PT"],
    #    xlabel="L2_PT",
    #    ylabel="Normalized Events",
    #    xmin=800,
    #    xmax=6000,
    #    nbins=50,
    #    output_file=f"../plots/MC_Data_agreement/R2/{tracktype}/L2_PT.pdf",
    #)

    feature = "ProbNNmu_max"
    pf.plot_feature(
        dfs=[df_MC, df_MC_corr, df_sw],
        normalize=True,
        styles=["E", "E", "E"],
        colors=["blue", "red", "black"],
        weights=[None, "total_weights", "signal_yield_sw"],
        names=["MC", "Corrected MC", "sWeighted Data"],
        markers=[20, 20, 20],
        marker_colors=["blue", "red", "black"],
        features=[feature, feature, feature],
        xlabel="ProbNNmu_{max}",
        ylabel="Normalized Events",
        xmin=xmin_dict[feature][tracktype],
        xmax=xmax_dict[feature][tracktype],
        nbins=50,
        output_file=f"../plots/MC_Data_agreement/R2/{tracktype}/ProbNNmu_max.pdf",
        log_scale=True,
        legend_position=[0.2, 0.7, 0.5, 0.9],
    )

    feature = "ProbNNmu_min"
    pf.plot_feature(
        dfs=[df_MC, df_MC_corr, df_sw],
        normalize=True,
        styles=["E", "E", "E"],
        colors=["blue", "red", "black"],
        weights=[None, "total_weights", "signal_yield_sw"],
        names=["MC", "Corrected MC", "sWeighted Data"],
        markers=[20, 20, 20],
        marker_colors=["blue", "red", "black"],
        features=[feature, feature, feature],
        xlabel="ProbNNmu_{min}",
        ylabel="Normalized Events",
        xmin=xmin_dict[feature][tracktype],
        xmax=xmax_dict[feature][tracktype],
        nbins=50,
        output_file=f"../plots/MC_Data_agreement/R2/{tracktype}/ProbNNmu_min.pdf",
        log_scale=True,
        legend_position=[0.2, 0.7, 0.5, 0.9],
    )

    feature = "log10_chi2ip_max"
    pf.plot_feature(
        dfs=[df_MC, df_MC_corr, df_sw],
        normalize=True,
        styles=["E", "E", "E"],
        colors=["blue", "red", "black"],
        weights=[None, "total_weights", "signal_yield_sw"],
        names=["MC", "Corrected MC", "sWeighted Data"],
        markers=[20, 20, 20],
        marker_colors=["blue", "red", "black"],
        features=[feature, feature, feature],
        xlabel="log_{10}(#chi^{2}_{IP, max})",
        ylabel="Normalized Events",
        xmin=xmin_dict[feature][tracktype],
        xmax=xmax_dict[feature][tracktype],
        nbins=50,
        output_file=f"../plots/MC_Data_agreement/R2/{tracktype}/log10_chi2ip_max.pdf",
        legend_position=[0.4, 0.2, 0.6, 0.5],
    )

    feature = "log10_chi2ip_min"
    pf.plot_feature(
        dfs=[df_MC, df_MC_corr, df_sw],
        normalize=True,
        styles=["E", "E", "E"],
        colors=["blue", "red", "black"],
        weights=[None, "total_weights", "signal_yield_sw"],
        names=["MC", "Corrected MC", "sWeighted Data"],
        markers=[20, 20, 20],
        marker_colors=["blue", "red", "black"],
        features=[feature, feature, feature],
        xlabel="log_{10}(#chi^{2}_{IP, min})",
        ylabel="Normalized Events",
        xmin=xmin_dict[feature][tracktype],
        xmax=xmax_dict[feature][tracktype],
        nbins=50,
        output_file=f"../plots/MC_Data_agreement/R2/{tracktype}/log10_chi2ip_min.pdf",
    )

    # mixed features
    feature = "log10_chi2ip_min_pi_p"
    pf.plot_feature(
        dfs=[df_MC, df_MC_corr, df_sw],
        normalize=True,
        styles=["E", "E", "E"],
        colors=["blue", "red", "black"],
        weights=[None, "total_weights", "signal_yield_sw"],
        names=["MC", "Corrected MC", "sWeighted Data"],
        markers=[20, 20, 20],
        marker_colors=["blue", "red", "black"],
        features=[feature, feature, feature],
        xlabel="log_{10}(#chi^{2}_{IP, min}(#pi, p))",
        ylabel="Normalized Events",
        xmin=xmin_dict[feature][tracktype],
        xmax=xmax_dict[feature][tracktype],
        nbins=50,
        output_file=f"../plots/MC_Data_agreement/R2/{tracktype}/log10_chi2ip_min_pi_p.pdf",
        legend_position=[0.4, 0.2, 0.6, 0.5],
    )

    # ProbNNmu
    pf.plot_feature(
        dfs=[df_MC, df_MC_corr, df_sw],
        normalize=True,
        styles=["E", "E", "E"],
        colors=["blue", "red", "black"],
        weights=[None, "total_weights", "signal_yield_sw"],
        names=["MC", "Corrected MC", "sWeighted Data"],
        markers=[20, 20, 20],
        marker_colors= ["blue", "red", "black"],
        features=["L1_ProbNNmu", "L1_ProbNNmu", "L1_ProbNNmu"],
        xlabel="ProbNNmu",
        ylabel="Normalized Events",
        xmin=0.1,
        xmax=1,
        nbins=50,
        output_file=f"../plots/MC_Data_agreement/R2/{tracktype}/ProbNNmu.pdf",
        legend_position=[0.2, 0.7, 0.5, 0.9],
        log_scale=True,
    )
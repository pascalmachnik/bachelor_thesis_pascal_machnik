import ROOT as r
import sys
import yaml
import os
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.append(parent_dir)
from utils.filter import df_filter
from utils.define import BDT_define, df_define
from utils.plot_functions import plot_feature, plot_feature_and_pull

# filenames
f_MC_R3 = "/ceph/users/pmachnik/Lambdab_Run3/merged/Lb2Lmumu_MC.root"
f_MC_R2 = [f"/ceph/users/pmachnik/Lambdab/MC/Lb2Lmm_PHMC_{year}_{magnet}_corr_all.root" 
           for year in ["2015", "2016", "2017", "2018"] 
           for magnet in ["MU", "MD"]]
f_MC_R2_unc = [f"/ceph/users/jnicolini/scripts/Lambdab/Samples/angles/Lb2Lmm_PHMC_{year}_{magnet}.root"
               for year in ["2015", "2016", "2017", "2018"]
               for magnet in ["MU", "MD"]]

# create dataframes
df_R2 = r.RDataFrame("DecayTree", f_MC_R2)
df_R2_unc = r.RDataFrame("DecayTree", f_MC_R2_unc)
df_R3 = r.RDataFrame("DecayTree", f_MC_R3)

# Apply filter on Run2 MC
df_R2 = df_filter(df_R2, "truth_matching_rd")
df_R2 = df_filter(df_R2, "preselection_corr_R2")
df_R2 = df_filter(df_R2, "trigger_R2")

df_R2_unc = df_filter(df_R2_unc, "truth_matching_rd")
df_R2_unc = df_filter(df_R2_unc, "preselection_R2")
df_R2_unc = df_filter(df_R2_unc, "trigger_R2")

# Apply filter on Run3 MC (Trigger are already applied in the Run3 data)
df_R3 = df_filter(df_R3, "truth_matching_rd")
df_R3 = df_filter(df_R3, "preselection")

# Define BDT variables for Run2 MC
df_R2 = BDT_define(df_R2, corrected=True, run3=False)
df_R2_unc = BDT_define(df_R2_unc, corrected=False, run3=False)

# Define BDT variables for Run3 MC
df_R3 = BDT_define(df_R3, corrected=False)

# High q2 filter
df_R2 = df_filter(df_R2, "high_q2")
df_R2_unc = df_filter(df_R2_unc, "high_q2")
df_R3 = df_filter(df_R3, "high_q2")

# print number of events in each dataframe
print(f"Run2 cor MC events: {df_R2.Count().GetValue()}")
print(f"Run2 unc MC events: {df_R2_unc.Count().GetValue()}")
print(f"Run3 unc MC events: {df_R3.Count().GetValue()}")

# Create ll and dd dfs and save them in a list
df_R2_ll = df_filter(df_R2, "ll_R2")
df_R2_dd = df_filter(df_R2, "dd_R2")
df_R2_unc_ll = df_filter(df_R2_unc, "ll_R2")
df_R2_unc_dd = df_filter(df_R2_unc, "dd_R2")
df_R3_ll = df_filter(df_R3, "ll")
df_R3_dd = df_filter(df_R3, "dd")
dfs_R2 = [df_R2_ll, df_R2_dd]
dfs_R3 = [df_R3_ll, df_R3_dd]
dfs_R2_unc = [df_R2_unc_ll, df_R2_unc_dd]
tracktype = ["ll", "dd"]

# print combined ll and dd events
print(f"Run2 unc MC ll events: {df_R2_unc_ll.Count().GetValue()}")
print(f"Run2 unc MC dd events: {df_R2_unc_dd.Count().GetValue()}")
print(f"Run3 unc MC ll events: {df_R3_ll.Count().GetValue()}")
print(f"Run3 unc MC dd events: {df_R3_dd.Count().GetValue()}")

# Start of plotting

xmin_dict = {
    "Lambdab_TAU": {"ll": 1e-4, "dd": 1e-4},
    "log10_Lambdab_pT": {"ll": 1.5, "dd": 1.5},
    "Lambdab_LOKI_DTF_CHI2NDOF": {"ll": 0, "dd": 0},
    "log10_Lambdab_DIRA_OWNPV": {"ll": -2.2e-4, "dd": -2.2e-4},
    "Lambda_FDCHI2_OWNPV": {"ll": 0, "dd": 0},
    "Lambda_FD_TOPPV": {"ll": 0, "dd": 200},
    "Lambda_IPCHI2_OWNPV": {"ll": 0, "dd": 0},
    "Pion_IPCHI2_OWNPV": {"ll": 0, "dd": 0},
    "log10_Proton_pT": {"ll": 2.4, "dd": 2.4},
    "log10_L1_pT": {"ll": 2.8, "dd": 2.8},
    "ProbNNmu_max": {"ll": 0, "dd": 0},
    "ProbNNmu_min": {"ll": 0, "dd": 0},
    "log10_chi2ip_max": {"ll": 0.8, "dd": 0.8},
    "log10_chi2ip_min": {"ll": 0.8, "dd": 0.8},
    "log10_chi2ip_min_pi_p": {"ll": 0.8, "dd": 0.8},
    "Lambda_TAU": {"ll": 0, "dd": 0},
    "L1_P": {"ll": 3, "dd": 3},
    "Lambda_M": {"ll": 1108, "dd": 1108},
    "Lambda_ENDVERTEX_Z": {"ll": 0, "dd": 200},
    "ProbNNmu": {"ll": 0, "dd": 0}, 
}
xmax_dict = {
    "Lambdab_TAU": {"ll": 8e-3, "dd": 8e-3},
    "log10_Lambdab_pT": {"ll": 5, "dd": 5},
    "Lambdab_LOKI_DTF_CHI2NDOF": {"ll": 20, "dd": 10},
    "log10_Lambdab_DIRA_OWNPV": {"ll": 0, "dd": 0},
    "Lambda_FDCHI2_OWNPV": {"ll": 3000, "dd": 3000},
    "Lambda_FD_TOPPV": {"ll": 750, "dd": 2500},
    "Lambda_IPCHI2_OWNPV": {"ll": 150, "dd": 50},
    "Pion_IPCHI2_OWNPV": {"ll": 3000, "dd": 2000},
    "log10_Proton_pT": {"ll": 4.2, "dd": 4.2},
    "log10_L1_pT": {"ll": 4.8, "dd": 4.8},
    "ProbNNmu_max": {"ll": 1, "dd": 1},
    "ProbNNmu_min": {"ll": 1, "dd": 1},
    "log10_chi2ip_max": {"ll": 5, "dd": 5},
    "log10_chi2ip_min": {"ll": 5, "dd": 5},
    "log10_chi2ip_min_pi_p": {"ll": 4, "dd": 4},
    "Lambda_TAU": {"ll": 2.1e-1, "dd": 5.1e-1},
    "L1_P": {"ll": 120, "dd": 120},
    "Lambda_M": {"ll": 1124, "dd": 1124},
    "Lambda_ENDVERTEX_Z": {"ll": 700, "dd": 2400},
    "ProbNNmu": {"ll": 1, "dd": 1},
}
for df_2, df_3, df_2_unc, tt in zip(dfs_R2, dfs_R3, dfs_R2_unc, tracktype):
    # Lambdab Tau
    feature = "Lambdab_TAU"
    xmin = xmin_dict[feature][tt]
    xmax = xmax_dict[feature][tt]
    plot_feature(
        dfs=[df_2, df_3, df_2_unc],
        normalize=True,
        styles=["E", "E", "E"],
        colors=["red", "black", "blue"],
        marker_colors=["red", "black", "blue"],
        weights=["total_weights", None, None],
        names=["Run2 cor MC", "Run3 unc MC", "Run2 unc MC"],
        markers=[20, 20, 20],
        features=["Lambdab_TAU", "Lambdab_BPVLTIME", "Lambdab_TAU"],
        xlabel="#tau(#Lambda^{0}_{b}) [ns]",
        ylabel="Normalized Events",
        xmin=xmin,
        xmax=xmax,
        nbins=50,
        legend_position=[0.6, 0.7, 0.8, 0.9],
        output_file=f"../plots/Run2_Run3_comparison/{tt}/Lambdab_lifetime.pdf",
    )

    # Lambdab Tau Chi2 (not used in Run3 but try Lambdab_DTF_Lamabda_CTAU)

    # Lambdab Endvertex Chi2
    #df_3 = df_3.Define("Lambdab_ENDVERTEX_CHI2", "Lambdab_ENDVERTEX_CHI2DOF * 4")
    #plot_feature(
    #    dfs=[df_2, df_3, df_2_unc],
    #    normalize=True,
    #    styles=["E", "E", "E"],
    #    colors=["red", "black", "blue"],
    #    marker_colors=["red", "black", "blue"],
    #    weights=["total_weights", None, None],
    #    names=["Run2 cor MC", "Run3 unc MC", "Run2 unc MC"],
    #    markers=[20, 20, 20],
    #    features=["Lambdab_ENDVERTEX_NDOF", "Lambdab_ENDVERTEX_CHI2DOF", "Lambdab_ENDVERTEX_NDOF"],
    #    xlabel="Lambdab endvertex Chi2",
    #    ylabel="Normalized Events",
    #    xmin=0,
    #    xmax=30,
    #    nbins=50,
    #    legend_position=[0.5, 0.7, 0.8, 0.9],
    #    output_file=f"../plots/Run2_Run3_comparison/{tt}/Lambdab_endvertex_chi2.pdf",
    #)

    # log10_Lambdab_pT
    feature = "log10_Lambdab_pT"
    xmin = xmin_dict[feature][tt]
    xmax = xmax_dict[feature][tt]
    plot_feature(
        dfs=[df_2, df_3, df_2_unc],
        normalize=True,
        styles=["E", "E", "E"],
        colors=["red", "black", "blue"],
        marker_colors=["red", "black", "blue"],
        weights=["total_weights", None, None],
        names=["Run2 cor MC", "Run3 unc MC", "Run2 unc MC"],
        markers=[20, 20, 20],
        features=["log10_Lambdab_pT", "log10_Lambdab_pT", "log10_Lambdab_pT"],
        xlabel="log10(p_{T}(#Lambda^{0}_{b}))",
        ylabel="Normalized Events",
        xmin=xmin,
        xmax=xmax,
        nbins=50,
        legend_position=[0.2, 0.7, 0.5, 0.9],
        output_file=f"../plots/Run2_Run3_comparison/{tt}/log10_Lambdab_pT.pdf",
    )

    # Lambadb DTF CHI2
    feature = "Lambdab_LOKI_DTF_CHI2NDOF"
    xmin = xmin_dict[feature][tt]
    xmax = xmax_dict[feature][tt]
    plot_feature(
        dfs=[df_2, df_3, df_2_unc],
        normalize=True,
        styles=["E", "E", "E"],
        colors=["red", "black", "blue"],
        marker_colors=["red", "black", "blue"],
        weights=["total_weights", None, None],
        names=["Run2 cor MC", "Run3 unc MC", "Run2 unc MC"],
        markers=[20, 20, 20],
        features=["Lambdab_LOKI_DTF_CHI2NDOF", "Lambdab_DTF_Lambda_CHI2DOF", "Lambdab_LOKI_DTF_CHI2NDOF"],
        xlabel="#chi^{2}_{m}(#Lambda^{0}_{b})",
        ylabel="Normalized Events",
        xmin=xmin,
        xmax=xmax,
        nbins=50,
        legend_position=[0.6, 0.7, 0.8, 0.9],
        output_file=f"../plots/Run2_Run3_comparison/{tt}/Lambdab_DTF_chi2.pdf",
    )

    # log10_Lambdab_BPVDIRA XXXX
    feature = "log10_Lambdab_DIRA_OWNPV"
    xmin = xmin_dict[feature][tt]
    xmax = xmax_dict[feature][tt]
    plot_feature(
        dfs=[df_2, df_3, df_2_unc],
        normalize=True,
        styles=["E", "E", "E"],
        colors=["red", "black", "blue"],
        marker_colors=["red", "black", "blue"],
        weights=["total_weights", None, None],
        names=["Run2 cor MC", "Run3 unc MC", "Run2 unc MC"],
        markers=[20, 20, 20],
        features=["log10_Lambdab_DIRA_OWNPV", "log10_Lambdab_BPVDIRA", "log10_Lambdab_DIRA_OWNPV"],
        xlabel="log10(#Theta_{DIRA}(#Lambda^{0}_{b}))",
        ylabel="Normalized Events",
        xmin=xmin,
        xmax=xmax,
        nbins=50,
        legend_position=[0.2, 0.7, 0.5, 0.9],
        output_file=f"../plots/Run2_Run3_comparison/{tt}/log10_Lambdab_BPVDIRA.pdf",
        log_scale=True,
    )

    # Lambda_FDCHI2_OWNPV
    feature = "Lambda_FDCHI2_OWNPV"
    xmin = xmin_dict[feature][tt]
    xmax = xmax_dict[feature][tt]
    plot_feature(
        dfs=[df_2, df_3, df_2_unc],
        normalize=True,
        styles=["E", "E", "E"],
        colors=["red", "black", "blue"],
        marker_colors=["red", "black", "blue"],
        weights=["total_weights", None, None],
        names=["Run2 cor MC", "Run3 unc MC", "Run2 unc MC"],
        markers=[20, 20, 20],
        features=["Lambda_FDCHI2_OWNPV", "Lambda_BPVFDCHI2", "Lambda_FDCHI2_OWNPV"],
        xlabel="#chi^{2}_{FD}(#Lambda^{0})",
        ylabel="Normalized Events",
        xmin=xmin,
        xmax=xmax,
        nbins=50,
        legend_position=[0.6, 0.7, 0.8, 0.9],
        output_file=f"../plots/Run2_Run3_comparison/{tt}/Lambda_FD_chi2.pdf",
    )

    # Lambda FD to PV
    feature = "Lambda_FD_TOPPV"
    xmin = xmin_dict[feature][tt]
    xmax = xmax_dict[feature][tt]
    plot_feature(
        dfs=[df_2, df_3, df_2_unc],
        normalize=True,
        styles=["E", "E", "E"],
        colors=["red", "black", "blue"],
        marker_colors=["red", "black", "blue"],
        weights=["total_weights", None, None],
        names=["Run2 cor MC", "Run3 unc MC", "Run2 unc MC"],
        markers=[20, 20, 20],
        features=["Lambda_FD_TOPPV", "Lambda_BPVFD", "Lambda_FD_TOPPV"],
        xlabel="FD(#Lambda^{0}) [mm]",
        ylabel="Normalized Events",
        xmin=xmin,
        xmax=xmax,
        nbins=50,
        legend_position=[0.6, 0.7, 0.8, 0.9],
        output_file=f"../plots/Run2_Run3_comparison/{tt}/Lambda_FD_PV.pdf",
    )

    # Lambda_IPCHI2_OWNPV
    feature = "Lambda_IPCHI2_OWNPV"
    xmin = xmin_dict[feature][tt]
    xmax = xmax_dict[feature][tt]
    plot_feature(
        dfs=[df_2, df_3, df_2_unc],
        normalize=True,
        styles=["E", "E", "E"],
        colors=["red", "black", "blue"],
        marker_colors=["red", "black", "blue"],
        weights=["total_weights", None, None],
        names=["Run2 cor MC", "Run3 unc MC", "Run2 unc MC"],
        markers=[20, 20, 20],
        features=["Lambda_IPCHI2_OWNPV", "Lambda_BPVIPCHI2", "Lambda_IPCHI2_OWNPV"],
        xlabel="#chi^{2}_{IP}(#Lambda^{0})",
        ylabel="Normalized Events",
        xmin=xmin,
        xmax=xmax,
        nbins=50,
        legend_position=[0.6, 0.7, 0.8, 0.9],
        output_file=f"../plots/Run2_Run3_comparison/{tt}/Lambda_IP_chi2.pdf",
        log_scale=True,
    )

    # Pion_IPCHI2_OWNPV
    feature = "Pion_IPCHI2_OWNPV"
    xmin = xmin_dict[feature][tt]
    xmax = xmax_dict[feature][tt]
    plot_feature(
        dfs=[df_2, df_3, df_2_unc],
        normalize=True,
        styles=["E", "E", "E"],
        colors=["red", "black", "blue"],
        marker_colors=["red", "black", "blue"],
        weights=["total_weights", None, None],
        names=["Run2 cor MC", "Run3 unc MC", "Run2 unc MC"],
        markers=[20, 20, 20],
        features=["Pion_IPCHI2_OWNPV", "Pion_BPVIPCHI2", "Pion_IPCHI2_OWNPV"],
        xlabel="#chi^{2}_{IP}(#pi)",
        ylabel="Normalized Events",
        xmin=xmin,
        xmax=xmax,
        nbins=50,
        legend_position=[0.6, 0.7, 0.8, 0.9],
        output_file=f"../plots/Run2_Run3_comparison/{tt}/Pion_IP_chi2.pdf",
    )

    # log10_Proton_pT
    feature = "log10_Proton_pT"
    xmin = xmin_dict[feature][tt]
    xmax = xmax_dict[feature][tt]
    plot_feature(
        dfs=[df_2, df_3, df_2_unc],
        normalize=True,
        styles=["E", "E", "E"],
        colors=["red", "black", "blue"],
        marker_colors=["red", "black", "blue"],
        weights=["total_weights", None, None],
        names=["Run2 cor MC", "Run3 unc MC", "Run2 unc MC"],
        markers=[20, 20, 20],
        features=["log10_Proton_pT", "log10_Proton_pT", "log10_Proton_pT"],
        xlabel="log10(p_{T}(p)) ",
        ylabel="Normalized Events",
        xmin=xmin,
        xmax=xmax,
        nbins=50,
        legend_position=[0.6, 0.7, 0.8, 0.9],
        output_file=f"../plots/Run2_Run3_comparison/{tt}/log10_Proton_pT.pdf",
    )

    # log10_L1_pT
    feature = "log10_L1_pT"
    xmin = xmin_dict[feature][tt]
    xmax = xmax_dict[feature][tt]
    plot_feature(
        dfs=[df_2, df_3, df_2_unc],
        normalize=True,
        styles=["E", "E", "E"],
        colors=["red", "black", "blue"],
        marker_colors=["red", "black", "blue"],
        weights=["total_weights", None, None],
        names=["Run2 cor MC", "Run3 unc MC", "Run2 unc MC"],
        markers=[20, 20, 20],
        features=["log10_L1_pT", "log10_L1_pT", "log10_L1_pT"],
        xlabel="log10(p_{T}(#mu))",
        ylabel="Normalized Events",
        xmin=xmin,
        xmax=xmax,
        nbins=50,
        legend_position=[0.6, 0.7, 0.8, 0.9],
        output_file=f"../plots/Run2_Run3_comparison/{tt}/log10_L1_pT.pdf",
    )

    # ProbNNmu_max
    feature = "ProbNNmu_max"
    xmin = xmin_dict[feature][tt]
    xmax = xmax_dict[feature][tt]
    plot_feature(
        dfs=[df_2, df_3, df_2_unc],
        normalize=True,
        styles=["E", "E", "E"],
        colors=["red", "black", "blue"],
        marker_colors=["red", "black", "blue"],
        weights=["total_weights", None, None],
        names=["Run2 cor MC", "Run3 unc MC", "Run2 unc MC"],
        markers=[20, 20, 20],
        features=["ProbNNmu_max", "ProbNNmu_max", "ProbNNmu_max"],
        xlabel="ProbNNmu_{max}",
        ylabel="Normalized Events",
        xmin=xmin,
        xmax=xmax,
        nbins=50,
        legend_position=[0.2, 0.7, 0.5, 0.9],
        output_file=f"../plots/Run2_Run3_comparison/{tt}/ProbNNmu_max.pdf",
        log_scale=True,
    )

    # ProbNNmu_min
    feature = "ProbNNmu_min"
    xmin = xmin_dict[feature][tt]
    xmax = xmax_dict[feature][tt]
    plot_feature(
        dfs=[df_2, df_3, df_2_unc],
        normalize=True,
        styles=["E", "E", "E"],
        colors=["red", "black", "blue"],
        marker_colors=["red", "black", "blue"],
        weights=["total_weights", None, None],
        names=["Run2 cor MC", "Run3 unc MC", "Run2 unc MC"],
        markers=[20, 20, 20],
        features=["ProbNNmu_min", "ProbNNmu_min", "ProbNNmu_min"],
        xlabel="ProbNNmu_{min}",
        ylabel="Normalized Events",
        xmin=xmin,
        xmax=xmax,
        nbins=50,
        legend_position=[0.2, 0.7, 0.5, 0.9],
        output_file=f"../plots/Run2_Run3_comparison/{tt}/ProbNNmu_min.pdf",
        log_scale=True,
    )

    # log10_chi2ip_max
    feature = "log10_chi2ip_max"
    xmin = xmin_dict[feature][tt]
    xmax = xmax_dict[feature][tt]
    plot_feature(
        dfs=[df_2, df_3, df_2_unc],
        normalize=True,
        styles=["E", "E", "E"],
        colors=["red", "black", "blue"],
        marker_colors=["red", "black", "blue"],
        weights=["total_weights", None, None],
        names=["Run2 cor MC", "Run3 unc MC", "Run2 unc MC"],
        markers=[20, 20, 20],
        features=["log10_chi2ip_max", "log10_chi2ip_max", "log10_chi2ip_max"],
        xlabel="log10(#chi^{2}_{IP,max}(#mu))",
        ylabel="Normalized Events",
        xmin=xmin,
        xmax=xmax,
        nbins=50,
        legend_position=[0.15, 0.7, 0.4, 0.9],
        output_file=f"../plots/Run2_Run3_comparison/{tt}/log10_chi2ip_max.pdf",
    )

    # log10_chi2ip_min
    feature = "log10_chi2ip_min"
    xmin = xmin_dict[feature][tt]
    xmax = xmax_dict[feature][tt]
    plot_feature(
        dfs=[df_2, df_3, df_2_unc],
        normalize=True,
        styles=["E", "E", "E"],
        colors=["red", "black", "blue"],
        marker_colors=["red", "black", "blue"],
        weights=["total_weights", None, None],
        names=["Run2 cor MC", "Run3 unc MC", "Run2 unc MC"],
        markers=[20, 20, 20],
        features=["log10_chi2ip_min", "log10_chi2ip_min", "log10_chi2ip_min"],
        xlabel="log10(#chi^{2}_{IP,min}(#mu))",
        ylabel="Normalized Events",
        xmin=xmin,
        xmax=xmax,
        nbins=50,
        legend_position=[0.6, 0.7, 0.8, 0.9],
        output_file=f"../plots/Run2_Run3_comparison/{tt}/log10_chi2ip_min.pdf",
    )

    # log10_chi2ip_min_pi_p
    feature = "log10_chi2ip_min_pi_p"
    xmin = xmin_dict[feature][tt]
    xmax = xmax_dict[feature][tt]
    legend_position = [0.15, 0.7, 0.4, 0.9] if tt == "ll" else [0.6, 0.7, 0.9, 0.9]
    plot_feature(
        dfs=[df_2, df_3, df_2_unc],
        normalize=True,
        styles=["E", "E", "E"],
        colors=["red", "black", "blue"],
        marker_colors=["red", "black", "blue"],
        weights=["total_weights", None, None],
        names=["Run2 cor MC", "Run3 unc MC", "Run2 unc MC"],
        markers=[20, 20, 20],
        features=["log10_chi2ip_min_pi_p", "log10_chi2ip_min_pi_p", "log10_chi2ip_min_pi_p"],
        xlabel="log10(#chi^{2}_{IP,min}(#pi, p))",
        ylabel="Normalized Events",
        xmin=xmin,
        xmax=xmax,
        nbins=50,
        legend_position= legend_position,
        output_file=f"../plots/Run2_Run3_comparison/{tt}/log10_chi2ip_min_pi_p.pdf",
    )

    # Lambda Tau
    feature = "Lambda_TAU"
    xmin = xmin_dict[feature][tt]
    xmax = xmax_dict[feature][tt]
    plot_feature(
        dfs=[df_2, df_3, df_2_unc],
        normalize=True,
        styles=["E", "E", "E"],
        colors=["red", "black", "blue"],
        marker_colors=["red", "black", "blue"],
        weights=["total_weights", None, None],
        names=["Run2 cor MC", "Run3 unc MC", "Run2 unc MC"],
        markers=[20, 20, 20],
        features=["Lambda_TAU", "Lambda_BPVLTIME", "Lambda_TAU"],
        xlabel="#tau(#Lambda^{0}) [ns]",
        ylabel="Normalized Events",
        xmin=xmin,
        xmax=xmax,
        nbins=50,
        legend_position=[0.6, 0.7, 0.8, 0.9],
        output_file=f"../plots/Run2_Run3_comparison/{tt}/Lambda_lifetime.pdf",
    )

    # L1 P
    feature = "L1_P"
    xmin = xmin_dict[feature][tt]
    xmax = xmax_dict[feature][tt]
    df_2 = df_2.Define("L1_P_GeV", "L1_P/1000")
    df_2_unc = df_2_unc.Define("L1_P_GeV", "L1_P/1000")
    df_3 = df_3.Define("L1_P_GeV", "L1_P/1000")
    plot_feature(
        dfs=[df_2, df_3, df_2_unc],
        normalize=True,
        styles=["E", "E", "E"],
        colors=["red", "black", "blue"],
        marker_colors=["red", "black", "blue"],
        weights=["total_weights", None, None],
        names=["Run2 cor MC", "Run3 unc MC", "Run2 unc MC"],
        markers=[20, 20, 20],
        features=["L1_P_GeV", "L1_P_GeV", "L1_P_GeV"],
        xlabel="p(#mu) [GeV/c]",
        ylabel="Normalized Events",
        xmin=xmin,
        xmax=xmax,
        nbins=50,
        legend_position=[0.6, 0.7, 0.8, 0.9],
        output_file=f"../plots/Run2_Run3_comparison/{tt}/L1_momentum.pdf",
    )

    # Lambda M
    feature = "Lambda_M"
    xmin = xmin_dict[feature][tt]
    xmax = xmax_dict[feature][tt]
    plot_feature(
        dfs=[df_2, df_3, df_2_unc],
        normalize=True,
        styles=["E", "E", "E"],
        colors=["red", "black", "blue"],
        marker_colors=["red", "black", "blue"],
        weights=["total_weights", None, None],
        names=["Run2 cor MC", "Run3 unc MC", "Run2 unc MC"],
        markers=[20, 20, 20],
        features=["Lambda_M", "Lambda_M", "Lambda_M"],
        xlabel="m_{#Lambda^{0}} [MeV/c^2]",
        ylabel="Normalized Events",
        xmin=xmin,
        xmax=xmax,
        nbins=50,
        legend_position=[0.6, 0.7, 0.8, 0.9],
        output_file=f"../plots/Run2_Run3_comparison/{tt}/Lambda_mass.pdf",
    )

    # Lambda Endvertex Z
    feature = "Lambda_ENDVERTEX_Z"
    xmin = xmin_dict[feature][tt]
    xmax = xmax_dict[feature][tt]
    plot_feature(
        dfs=[df_2, df_3, df_2_unc],
        normalize=True,
        styles=["E", "E", "E"],
        colors=["red", "black", "blue"],
        marker_colors=["red", "black", "blue"],
        weights=["total_weights", None, None],
        names=["Run2 cor MC", "Run3 unc MC", "Run2 unc MC"],
        markers=[20, 20, 20],
        features=["Lambda_ENDVERTEX_Z", "Lambda_END_VZ", "Lambda_ENDVERTEX_Z"],
        xlabel="z_{endvertex} [mm]",
        ylabel="Normalized Events",
        xmin=xmin,
        xmax=xmax,
        nbins=50,
        legend_position=[0.6, 0.7, 0.8, 0.9],
        output_file=f"../plots/Run2_Run3_comparison/{tt}/Lambda_endvertex_z.pdf",
    )

    # ProbNNmu
    feature = "ProbNNmu"
    xmin = xmin_dict[feature][tt]
    xmax = xmax_dict[feature][tt]
    plot_feature(
        dfs=[df_2, df_3, df_2_unc],
        normalize=True,
        styles=["E", "E", "E"],
        colors=["red", "black", "blue"],
        marker_colors=["red", "black", "blue"],
        weights=["total_weights", None, None],
        names=["Run2 cor MC", "Run3 unc MC", "Run2 unc MC"],
        markers=[20, 20, 20],
        features=["L1_ProbNNmu_pidcorr_default", "L1_PROBNN_MU", "L1_ProbNNmu"],
        xlabel="ProbNNmu",
        ylabel="Normalized Events",
        xmin=xmin,
        xmax=xmax,
        nbins=50,
        legend_position=[0.2, 0.7, 0.4, 0.9],
        output_file=f"../plots/Run2_Run3_comparison/{tt}/ProbNNmu.pdf",
        log_scale=True,
    )
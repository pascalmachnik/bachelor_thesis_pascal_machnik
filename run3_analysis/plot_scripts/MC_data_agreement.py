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
f_MC = "/ceph/users/pmachnik/Lambdab_Run3/merged/Lb2LJpsi_MC.root"
f_sWeights_ll = "/ceph/users/pmachnik/Lambdab_Run3/Jpsi/Jpsi_sWeights_ll.root"
f_sWeights_dd = "/ceph/users/pmachnik/Lambdab_Run3/Jpsi/Jpsi_sWeights_dd.root"

# create dataframes
df_MC = r.RDataFrame("DecayTree", f_MC)
df_sWeights_R3_ll = r.RDataFrame("DecayTree", f_sWeights_ll)
df_sWeights_R3_dd = r.RDataFrame("DecayTree", f_sWeights_dd)

# Apply filter on Run3 MC (Trigger are already applied in the Run3 data)
df_MC = df_filter(df_MC, "truth_matching_Jpsi")
df_MC = df_filter(df_MC, "preselection")
df_MC = df_filter(df_MC, "fit_range")
df_MC = df_filter(df_MC, "Jpsi")

# Define BDT variables for Run3 MC
df_MC = BDT_define(df_MC, corrected=False)

# Create ll and dd dataframes for MC
df_MC_ll = df_filter(df_MC, "ll")
df_MC_dd = df_filter(df_MC, "dd")

# Create lists of dataframes
tracktype = ["ll", "dd"]
dfs_MC = [df_MC_ll, df_MC_dd]
dfs_sWeights = [df_sWeights_R3_ll, df_sWeights_R3_dd]

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
    "log10_Lambdab_pT": {"ll": 4.8, "dd": 4.8},
    "Lambdab_LOKI_DTF_CHI2NDOF": {"ll": 20, "dd": 10},
    "log10_Lambdab_DIRA_OWNPV": {"ll": 0, "dd": 0},
    "Lambda_FDCHI2_OWNPV": {"ll": 3000, "dd": 3000},
    "Lambda_FD_TOPPV": {"ll": 750, "dd": 2500},
    "Lambda_IPCHI2_OWNPV": {"ll": 150, "dd": 20},
    "Pion_IPCHI2_OWNPV": {"ll": 3000, "dd": 2000},
    "log10_Proton_pT": {"ll": 4, "dd": 4},
    "log10_L1_pT": {"ll": 4.6, "dd": 4.6},
    "ProbNNmu_max": {"ll": 1, "dd": 1},
    "ProbNNmu_min": {"ll": 1, "dd": 1},
    "log10_chi2ip_max": {"ll": 5, "dd": 5},
    "log10_chi2ip_min": {"ll": 3.8, "dd": 3.8},
    "log10_chi2ip_min_pi_p": {"ll": 3.5, "dd": 3.5},
    "Lambda_TAU": {"ll": 2.1e-1, "dd": 5.1e-1},
    "L1_P": {"ll": 120, "dd": 120},
    "Lambda_M": {"ll": 1124.5, "dd": 1124.5},
    "Lambda_ENDVERTEX_Z": {"ll": 700, "dd": 2450},
    "ProbNNmu": {"ll": 1, "dd": 1},
}

for MC, sWeights, tt in zip(dfs_MC, dfs_sWeights, tracktype):
    # Lambdab Tau
    feature = "Lambdab_TAU"
    xmin = xmin_dict[feature][tt]
    xmax = xmax_dict[feature][tt]
    plot_feature(
        dfs=[MC, sWeights],
        normalize=True,
        styles=["E", "E"],
        colors=["red", "black"],
        marker_colors=["red", "black"],
        weights=[None, "signal_yield_sw"],
        names=["MC", "sWeighted data"],
        markers=[20, 20],
        features=["Lambdab_BPVLTIME","Lambdab_BPVLTIME"],
        xlabel="#tau(#Lambda^{0}_{b}) [ns]",
        ylabel="Normalized Events",
        xmin=xmin,
        xmax=xmax,
        nbins=50,
        legend_position=[0.6, 0.7, 0.8, 0.9],
        output_file=f"../plots/MC_Data_agreement/R3/{tt}/Lambdab_lifetime.pdf",
    )

    # Lambdab Tau Chi2 (not used in Run3 but try Lambdab_DTF_Lamabda_CTAU)

    # Lambdab Endvertex Chi2
    #sWeights = sWeights.Define("Lambdab_ENDVERTEX_CHI2", "Lambdab_ENDVERTEX_CHI2DOF * 4")
    #plot_feature(
    #    dfs=[MC, sWeights],
    #    normalize=True,
    #    styles=["HIST", "E"],
    #    colors=["red", "black"],
    #    marker_colors=["red", "black"],
    #    weights=[None, "signal_yield_sw"],
    #    names=["MC", "sWeighted data"],
    #    markers=[20, 24],
    #    features=["Lambdab_ENDVERTEX_NDOF", "Lambdab_ENDVERTEX_CHI2DOF", "Lambdab_ENDVERTEX_NDOF"],
    #    xlabel="Lambdab endvertex Chi2",
    #    ylabel="Normalized Events",
    #    xmin=0,
    #    xmax=30,
    #    nbins=50,
    #    legend_position=[0.5, 0.7, 0.8, 0.9],
    #    output_file=f"../plots/MC_Data_agreement/R3/{tt}/Lambdab_endvertex_chi2.pdf",
    #)

    # log10_Lambdab_pT
    feature = "log10_Lambdab_pT"
    xmin = xmin_dict[feature][tt]
    xmax = xmax_dict[feature][tt]
    plot_feature(
        dfs=[MC, sWeights],
        normalize=True,
        styles=["E", "E"],
        colors=["red", "black"],
        marker_colors=["red", "black"],
        weights=[None, "signal_yield_sw"],
        names=["MC", "sWeighted data"],
        markers=[20, 20],
        features=["log10_Lambdab_pT", "log10_Lambdab_pT"],
        xlabel="log_{10}(p_{T}(#Lambda^{0}_{b}))",
        ylabel="Normalized Events",
        xmin=xmin,
        xmax=xmax,
        nbins=50,
        legend_position=[0.2, 0.7, 0.5, 0.9],
        output_file=f"../plots/MC_Data_agreement/R3/{tt}/log10_Lambdab_pT.pdf",
    )

    # Lambadb DTF CHI2
    feature = "Lambdab_LOKI_DTF_CHI2NDOF"
    xmin = xmin_dict[feature][tt]
    xmax = xmax_dict[feature][tt]
    plot_feature(
        dfs=[MC, sWeights],
        normalize=True,
        styles=["E", "E"],
        colors=["red", "black"],
        marker_colors=["red", "black"],
        weights=[None, "signal_yield_sw"],
        names=["MC", "sWeighted data"],
        markers=[20, 20],
        features=["Lambdab_DTF_Lambda_CHI2DOF", "Lambdab_DTF_Lambda_CHI2DOF"],
        xlabel="#chi^{2}_{m}(#Lambda^{0}_{b})",
        ylabel="Normalized Events",
        xmin=xmin,
        xmax=xmax,
        nbins=50,
        legend_position=[0.6, 0.7, 0.8, 0.9],
        output_file=f"../plots/MC_Data_agreement/R3/{tt}/Lambdab_DTF_chi2.pdf",
    )

    # log10_Lambdab_BPVDIRA XXXX
    feature = "log10_Lambdab_DIRA_OWNPV"
    xmin = xmin_dict[feature][tt]
    xmax = xmax_dict[feature][tt]
    plot_feature(
        dfs=[MC, sWeights],
        normalize=True,
        styles=["E", "E"],
        colors=["red", "black"],
        marker_colors=["red", "black"],
        weights=[None, "signal_yield_sw"],
        names=["MC", "sWeighted data"],
        markers=[20, 20],
        features=["log10_Lambdab_BPVDIRA", "log10_Lambdab_BPVDIRA"],
        xlabel="log_{10}(#Theta_{DIRA}(#Lambda^{0}_{b}))",
        ylabel="Normalized Events",
        xmin=xmin,
        xmax=xmax,
        nbins=50,
        legend_position=[0.2, 0.7, 0.5, 0.9],
        output_file=f"../plots/MC_Data_agreement/R3/{tt}/log10_Lambdab_BPVDIRA.pdf",
        log_scale=True,
    )

    # Lambda_FDCHI2_OWNPV
    feature = "Lambda_FDCHI2_OWNPV"
    xmin = xmin_dict[feature][tt]
    xmax = xmax_dict[feature][tt]
    plot_feature(
        dfs=[MC, sWeights],
        normalize=True,
        styles=["E", "E"],
        colors=["red", "black"],
        marker_colors=["red", "black"],
        weights=[None, "signal_yield_sw"],
        names=["MC", "sWeighted data"],
        markers=[20, 20],
        features=["Lambda_BPVFDCHI2", "Lambda_BPVFDCHI2"],
        xlabel="#chi^{2}_{FD}(#Lambda^{0})",
        ylabel="Normalized Events",
        xmin=xmin,
        xmax=xmax,
        nbins=50,
        legend_position=[0.6, 0.7, 0.8, 0.9],
        output_file=f"../plots/MC_Data_agreement/R3/{tt}/Lambda_FD_chi2.pdf",
    )

    # Lambda FD to PV
    feature = "Lambda_FD_TOPPV"
    xmin = xmin_dict[feature][tt]
    xmax = xmax_dict[feature][tt]
    legend_position = [0.6, 0.7, 0.8, 0.9] if tt == "ll" else [0.45, 0.3, 0.5, 0.6]
    plot_feature(
        dfs=[MC, sWeights],
        normalize=True,
        styles=["E", "E"],
        colors=["red", "black"],
        marker_colors=["red", "black"],
        weights=[None, "signal_yield_sw"],
        names=["MC", "sWeighted data"],
        markers=[20, 20],
        features=["Lambda_BPVFD", "Lambda_BPVFD"],
        xlabel="FD(#Lambda^{0}) [mm]",
        ylabel="Normalized Events",
        xmin=xmin,
        xmax=xmax,
        nbins=50,
        legend_position=legend_position,
        output_file=f"../plots/MC_Data_agreement/R3/{tt}/Lambda_FD_PV.pdf",
    )

    # Lambda_IPCHI2_OWNPV
    feature = "Lambda_IPCHI2_OWNPV"
    xmin = xmin_dict[feature][tt]
    xmax = xmax_dict[feature][tt]
    plot_feature(
        dfs=[MC, sWeights],
        normalize=True,
        styles=["E", "E"],
        colors=["red", "black"],
        marker_colors=["red", "black"],
        weights=[None, "signal_yield_sw"],
        names=["MC", "sWeighted data"],
        markers=[20, 20],
        features=["Lambda_BPVIPCHI2", "Lambda_BPVIPCHI2"],
        xlabel="#chi^{2}_{IP}(#Lambda^{0})",
        ylabel="Normalized Events",
        xmin=xmin,
        xmax=xmax,
        nbins=50,
        legend_position=[0.6, 0.7, 0.8, 0.9],
        output_file=f"../plots/MC_Data_agreement/R3/{tt}/Lambda_IP_chi2.pdf",
    )

    # Pion_IPCHI2_OWNPV
    feature = "Pion_IPCHI2_OWNPV"
    xmin = xmin_dict[feature][tt]
    xmax = xmax_dict[feature][tt]
    plot_feature(
        dfs=[MC, sWeights],
        normalize=True,
        styles=["E", "E"],
        colors=["red", "black"],
        marker_colors=["red", "black"],
        weights=[None, "signal_yield_sw"],
        names=["MC", "sWeighted data"],
        markers=[20, 20],
        features=["Pion_BPVIPCHI2", "Pion_BPVIPCHI2"],
        xlabel="#chi^{2}_{IP}(#pi)",
        ylabel="Normalized Events",
        xmin=xmin,
        xmax=xmax,
        nbins=50,
        legend_position=[0.6, 0.7, 0.8, 0.9],
        output_file=f"../plots/MC_Data_agreement/R3/{tt}/Pion_IP_chi2.pdf",
    )

    # log10_Proton_pT
    feature = "log10_Proton_pT"
    xmin = xmin_dict[feature][tt]
    xmax = xmax_dict[feature][tt]
    plot_feature(
        dfs=[MC, sWeights],
        normalize=True,
        styles=["E", "E"],
        colors=["red", "black"],
        marker_colors=["red", "black"],
        weights=[None, "signal_yield_sw"],
        names=["MC", "sWeighted data"],
        markers=[20, 20],
        features=["log10_Proton_pT", "log10_Proton_pT"],
        xlabel="log_{10}(p_{T}(p)) ",
        ylabel="Normalized Events",
        xmin=xmin,
        xmax=xmax,
        nbins=50,
        legend_position=[0.2, 0.7, 0.6, 0.9],
        output_file=f"../plots/MC_Data_agreement/R3/{tt}/log10_Proton_pT.pdf",
    )

    # log10_L1_pT
    feature = "log10_L1_pT"
    xmin = xmin_dict[feature][tt]
    xmax = xmax_dict[feature][tt]
    plot_feature(
        dfs=[MC, sWeights],
        normalize=True,
        styles=["E", "E"],
        colors=["red", "black"],
        marker_colors=["red", "black"],
        weights=[None, "signal_yield_sw"],
        names=["MC", "sWeighted data"],
        markers=[20, 20],
        features=["log10_L1_pT", "log10_L1_pT"],
        xlabel="log_{10}(p_{T}(#mu))",
        ylabel="Normalized Events",
        xmin=xmin,
        xmax=xmax,
        nbins=50,
        legend_position=[0.6, 0.7, 0.8, 0.9],
        output_file=f"../plots/MC_Data_agreement/R3/{tt}/log10_L1_pT.pdf",
    )

    # ProbNNmu_max
    feature = "ProbNNmu_max"
    xmin = xmin_dict[feature][tt]
    xmax = xmax_dict[feature][tt]
    plot_feature(
        dfs=[MC, sWeights],
        normalize=True,
        styles=["E", "E"],
        colors=["red", "black"],
        marker_colors=["red", "black"],
        weights=[None, "signal_yield_sw"],
        names=["MC", "sWeighted data"],
        markers=[20, 20],
        features=["ProbNNmu_max", "ProbNNmu_max"],
        xlabel="ProbNNmu_{max}",
        ylabel="Normalized Events",
        xmin=xmin,
        xmax=xmax,
        nbins=50,
        legend_position=[0.2, 0.7, 0.5, 0.9],
        output_file=f"../plots/MC_Data_agreement/R3/{tt}/ProbNNmu_max.pdf",
        log_scale=True,
    )

    # ProbNNmu_min
    feature = "ProbNNmu_min"
    xmin = xmin_dict[feature][tt]
    xmax = xmax_dict[feature][tt]
    plot_feature(
        dfs=[MC, sWeights],
        normalize=True,
        styles=["E", "E"],
        colors=["red", "black"],
        marker_colors=["red", "black"],
        weights=[None, "signal_yield_sw"],
        names=["MC", "sWeighted data"],
        markers=[20, 20],
        features=["ProbNNmu_min", "ProbNNmu_min"],
        xlabel="ProbNNmu_{min}",
        ylabel="Normalized Events",
        xmin=xmin,
        xmax=xmax,
        nbins=50,
        legend_position=[0.2, 0.7, 0.5, 0.9],
        output_file=f"../plots/MC_Data_agreement/R3/{tt}/ProbNNmu_min.pdf",
        log_scale=True,
    )

    # log10_chi2ip_max
    feature = "log10_chi2ip_max"
    xmin = xmin_dict[feature][tt]
    xmax = xmax_dict[feature][tt]
    plot_feature(
        dfs=[MC, sWeights],
        normalize=True,
        styles=["E", "E"],
        colors=["red", "black"],
        marker_colors=["red", "black"],
        weights=[None, "signal_yield_sw"],
        names=["MC", "sWeighted data"],
        markers=[20, 20],
        features=["log10_chi2ip_max", "log10_chi2ip_max"],
        xlabel="log_{10}(#chi^{2}_{IP, max}(#mu))",
        ylabel="Normalized Events",
        xmin=xmin,
        xmax=xmax,
        nbins=50,
        legend_position=[0.4, 0.3, 0.6, 0.6],
        output_file=f"../plots/MC_Data_agreement/R3/{tt}/log10_chi2ip_max.pdf",
    )

    # log10_chi2ip_min
    feature = "log10_chi2ip_min"
    xmin = xmin_dict[feature][tt]
    xmax = xmax_dict[feature][tt]
    plot_feature(
        dfs=[MC, sWeights],
        normalize=True,
        styles=["E", "E"],
        colors=["red", "black"],
        marker_colors=["red", "black"],
        weights=[None, "signal_yield_sw"],
        names=["MC", "sWeighted data"],
        markers=[20, 20],
        features=["log10_chi2ip_min", "log10_chi2ip_min"],
        xlabel="log_{10}(#chi^{2}_{IP, min}(#mu))",
        ylabel="Normalized Events",
        xmin=xmin,
        xmax=xmax,
        nbins=50,
        legend_position=[0.4, 0.3, 0.6, 0.6],
        output_file=f"../plots/MC_Data_agreement/R3/{tt}/log10_chi2ip_min.pdf",
    )

    # log10_chi2ip_min_pi_p
    feature = "log10_chi2ip_min_pi_p"
    xmin = xmin_dict[feature][tt]
    xmax = xmax_dict[feature][tt]
    legend_position = [0.45, 0.3, 0.5, 0.6]
    plot_feature(
        dfs=[MC, sWeights],
        normalize=True,
        styles=["E", "E"],
        colors=["red", "black"],
        marker_colors=["red", "black"],
        weights=[None, "signal_yield_sw"],
        names=["MC", "sWeighted data"],
        markers=[20, 20],
        features=["log10_chi2ip_min_pi_p", "log10_chi2ip_min_pi_p"],
        xlabel="log_{10}(#chi^{2}_{IP, min}(#pi, p))",
        ylabel="Normalized Events",
        xmin=xmin,
        xmax=xmax,
        nbins=50,
        legend_position= legend_position,
        output_file=f"../plots/MC_Data_agreement/R3/{tt}/log10_chi2ip_min_pi_p.pdf",
    )

    # Lambda Tau
    feature = "Lambda_TAU"
    xmin = xmin_dict[feature][tt]
    xmax = xmax_dict[feature][tt]
    plot_feature(
        dfs=[MC, sWeights],
        normalize=True,
        styles=["E", "E"],
        colors=["red", "black"],
        marker_colors=["red", "black"],
        weights=[None, "signal_yield_sw"],
        names=["MC", "sWeighted data"],
        markers=[20, 20],
        features=["Lambda_BPVLTIME", "Lambda_BPVLTIME"],
        xlabel="#tau(#Lambda^{0}) [ns]",
        ylabel="Normalized Events",
        xmin=xmin,
        xmax=xmax,
        nbins=50,
        legend_position=[0.6, 0.7, 0.8, 0.9],
        output_file=f"../plots/MC_Data_agreement/R3/{tt}/Lambda_lifetime.pdf",
    )

    # L1 P
    feature = "L1_P"
    xmin = xmin_dict[feature][tt]
    xmax = xmax_dict[feature][tt]
    MC = MC.Define("L1_P_GeV", "L1_P/1000")
    sWeights = sWeights.Define("L1_P_GeV", "L1_P/1000")
    plot_feature(
        dfs=[MC, sWeights],
        normalize=True,
        styles=["E", "E"],
        colors=["red", "black"],
        marker_colors=["red", "black"],
        weights=[None, "signal_yield_sw"],
        names=["MC", "sWeighted data"],
        markers=[20, 20],
        features=["L1_P_GeV", "L1_P_GeV"],
        xlabel="p(#mu) [GeV/c]",
        ylabel="Normalized Events",
        xmin=xmin,
        xmax=xmax,
        nbins=50,
        legend_position=[0.6, 0.7, 0.8, 0.9],
        output_file=f"../plots/MC_Data_agreement/R3/{tt}/L1_momentum.pdf",
    )

    # Lambda M
    feature = "Lambda_M"
    xmin = xmin_dict[feature][tt]
    xmax = xmax_dict[feature][tt]
    plot_feature(
        dfs=[MC, sWeights],
        normalize=True,
        styles=["E", "E"],
        colors=["red", "black"],
        marker_colors=["red", "black"],
        weights=[None, "signal_yield_sw"],
        names=["MC", "sWeighted data"],
        markers=[20, 20],
        features=["Lambda_M", "Lambda_M"],
        xlabel= "m_{#Lambda^{0}} [MeV/c^2]",
        ylabel="Normalized Events",
        xmin=xmin,
        xmax=xmax,
        nbins=50,
        legend_position=[0.6, 0.7, 0.8, 0.9],
        output_file=f"../plots/MC_Data_agreement/R3/{tt}/Lambda_mass.pdf",
    )

    # Lambda Endvertex Z
    feature = "Lambda_ENDVERTEX_Z"
    xmin = xmin_dict[feature][tt]
    xmax = xmax_dict[feature][tt]
    plot_feature(
        dfs=[MC, sWeights],
        normalize=True,
        styles=["E", "E"],
        colors=["red", "black"],
        marker_colors=["red", "black"],
        weights=[None, "signal_yield_sw"],
        names=["MC", "sWeighted data"],
        markers=[20, 20],
        features=["Lambda_END_VZ", "Lambda_END_VZ"],
        xlabel="z_{endvertex}(#Lambda^{0}) [mm]",
        ylabel="Normalized Events",
        xmin=xmin,
        xmax=xmax,
        nbins=50,
        legend_position=[0.45, 0.3, 0.6, 0.6],
        output_file=f"../plots/MC_Data_agreement/R3/{tt}/Lambda_endvertex_z.pdf",
    )

    # ProbNNmu
    feature = "ProbNNmu"
    xmin = xmin_dict[feature][tt]
    xmax = xmax_dict[feature][tt]
    plot_feature(
        dfs=[MC, sWeights],
        normalize=True,
        styles=["E", "E"],
        colors=["red", "black"],
        marker_colors=["red", "black"],
        weights=[None, "signal_yield_sw"],
        names=["MC", "sWeighted data"],
        markers=[20, 20],
        features=["L1_PROBNN_MU", "L1_PROBNN_MU"],
        xlabel="ProbNNmu",
        ylabel="Normalized Events",
        xmin=xmin,
        xmax=xmax,
        nbins=50,
        legend_position=[0.2, 0.7, 0.4, 0.9],
        output_file=f"../plots/MC_Data_agreement/R3/{tt}/ProbNNmu.pdf",
        log_scale=True,
    )
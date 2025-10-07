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

# file names
#file_ll = "/ceph/users/pmachnik/Lambdab/BDT/RM_sWeights_ll_1516_BDTG.root"
#file_dd = "/ceph/users/pmachnik/Lambdab/BDT/RM_sWeights_dd_1516_BDTG.root"
#file_MC_ll = "/ceph/users/pmachnik/Lambdab/BDT/RM_MC_1516_BDTG_ll.root"
#file_MC_dd = "/ceph/users/pmachnik/Lambdab/BDT/RM_MC_1516_BDTG_dd.root"

file_MC_ll = "/ceph/users/pmachnik/Lambdab/BDT/MC_ll_BDTG.root"
file_MC_dd = "/ceph/users/pmachnik/Lambdab/BDT/MC_dd_BDTG.root"
file_ll = "/ceph/users/pmachnik/Lambdab/BDT/Jpsi_sWeights_ll_BDTG.root"
file_dd = "/ceph/users/pmachnik/Lambdab/BDT/Jpsi_sWeights_dd_BDTG.root"

# create dfs 
df_ll = r.RDataFrame("DecayTree", file_ll)
df_dd = r.RDataFrame("DecayTree", file_dd)
df_MC_ll = r.RDataFrame("DecayTree", file_MC_ll)
df_MC_dd = r.RDataFrame("DecayTree", file_MC_dd)

# Match the mass range of the MC and data
df_MC_ll = df_MC_ll.Filter("Lambdab_LOKI_MASS_LConstr > 5300 && Lambdab_LOKI_MASS_LConstr < 5940")
df_MC_dd = df_MC_dd.Filter("Lambdab_LOKI_MASS_LConstr > 5300 && Lambdab_LOKI_MASS_LConstr < 5940")


# plot the BDTG score for ll and dd
pf.plot_feature(
    dfs = [df_ll, df_MC_ll],
    normalize = True,
    styles = ["E", "E"],
    colors = ["black", "red"],
    marker_colors= ["black", "red"],
    weights = ["signal_yield_sw", "total_weights"],
    names = ["sWeighted Data", "Corrected MC"],
    markers = [20, 22],
    features = ["BDTG", "BDTG"],
    xlabel="BDTG score",
    ylabel="Normalized Events",
    xmin=-1,
    xmax=1,
    nbins=20,
    output_file="../plots/BDT/BDTG_score_ll.pdf",
)

pf.plot_feature(
    dfs = [df_dd, df_MC_dd],
    normalize = True,
    styles = ["E", "E"],
    colors = ["black", "red"],
    marker_colors= ["black", "red"],
    weights = ["signal_yield_sw", "total_weights"],
    names = ["sWeighted Data", "Corrected MC"],
    markers = [20, 22],
    features = ["BDTG", "BDTG"],
    xlabel="BDTG score",
    ylabel="Normalized Events",
    xmin=-1,
    xmax=1,
    nbins=20,
    output_file="../plots/BDT/BDTG_score_dd.pdf",
)

# plot the Lambab mass for ll data
pf.plot_feature(
    dfs = [df_ll],
    styles= ["HIST"],
    colors = ["red"],
    weights = ["signal_yield_sw"],
    names = ["sWeighted Data"],
    markers = [22],
    features = ["Lambdab_LOKI_MASS_AllConstr"],
    xlabel="m_{#Lambda_{b}} [MeV]",
    ylabel="sWeighted Events",
    xmin=5400,
    xmax=5800,
    nbins=50,
    output_file="../plots/BDT/Lambdab_mass_ll.pdf",
)
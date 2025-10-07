import ROOT as r
import yaml
from utils.filter import df_filter

# read the configuration files
with open("analysis_config.yml", "r") as analysis_config_file:
    config = yaml.safe_load(analysis_config_file)

file_MC_ll = "/ceph/users/pmachnik/Lambdab_Run3/BDT/MC_ll_BDTG.root"
file_MC_dd = "/ceph/users/pmachnik/Lambdab_Run3/BDT/MC_dd_BDTG.root"

# create dfs 
df_MC_ll = r.RDataFrame("DecayTree", file_MC_ll)
df_MC_dd = r.RDataFrame("DecayTree", file_MC_dd)

# calculate signal efficiency
count_ll = df_MC_ll.Count().GetValue()
count_dd = df_MC_dd.Count().GetValue()

df_MC_ll = df_filter(df_MC_ll, "BDT_selection_ll")
df_MC_dd = df_filter(df_MC_dd, "BDT_selection_dd")

count_ll_sel = df_MC_ll.Count().GetValue()
count_dd_sel = df_MC_dd.Count().GetValue()

# print results
print(f"LL BDT efficiency: {count_ll_sel / count_ll}")
print(f"DD BDT efficiency: {count_dd_sel / count_dd}")
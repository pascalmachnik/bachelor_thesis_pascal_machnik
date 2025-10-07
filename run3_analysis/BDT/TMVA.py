import ROOT as r
import argparse
import yaml

# Argparse setup
parser = argparse.ArgumentParser(description="Input and output file")
parser.add_argument("--input_signal", action="store", dest="input_signal", type=str, required=True)
parser.add_argument("--input_background", action="store", dest="input_background", type=str, required=True)
parser.add_argument("--features", action="store", nargs="+", dest="features", type=str, required=True)
args = parser.parse_args()

# read configuration file
with open("analysis_config.yml", "r") as config_file:
    config = yaml.safe_load(config_file)

r.TMVA.Tools.Instance()
r.EnableImplicitMT(16)

# Read input files
file1 = r.TFile(args.input_signal)
file2 = r.TFile(args.input_background)
tree_s = file1.Get("DecayTree")
tree_b = file2.Get("DecayTree")

# Set hyperparameters
track_type = args.input_signal[-7:-5]  # "ll" or "dd"
trees = config[f"{track_type}_trees"] #[200] 120, 160, 220, 240,
learning_rates = config[f"{track_type}_learning_rate"] # 0.12, 0.11  0.14, 0.16, 0.2, 0.3
depths = config[f"{track_type}_depth"]  # 3,4,5
# index = [9]

# for i in index: include to have a kFold
for N_TREES in trees:
    for learning_rate in learning_rates:
        for depth in depths:
            print("Number of Trees:", N_TREES, "Shrinkage:", learning_rate, "MaxDepth:", depth)
            # define output file
            fout = r.TFile(
                "TMVA_N{}_Sh{}_D{}_alldata_{}.root".format(N_TREES, learning_rate, depth, track_type), "RECREATE"
            )  # _kfold{}, i

            # define factory with options
            factory = r.TMVA.Factory(
                "TMVAClassification",
                fout,
                ":".join(
                    [
                        "!V",
                        "!Silent",
                        "Color",
                        "DrawProgressBar",
                        "Transformations=I;D;P;G,D",
                        "AnalysisType=Classification",
                    ]
                ),
            )

            # add discriminating variables for training
            dataloader = r.TMVA.DataLoader(
                "dataset_N{}_Sh{}_D{}_{}".format(N_TREES, learning_rate, depth,track_type)
            )  # _kfold{}, i

            for feature in args.features:
                dataloader.AddVariable(feature, feature, "", "F")
            

           #dataloader.AddVariable(
           #    "minIPCHI2_PPi:=min(log10(abs(Proton_IPCHI2_OWNPV)), log10(abs(Pion_IPCHI2_OWNPV)))",
           #    "log_minIPCHI2_PPi",
           #    "",
           #    "F",
           #)
           #dataloader.AddVariable(
           #    "L12_IPCHI2_OWNPV:=min(log10(abs(L1_IPCHI2_OWNPV)), log10(abs(L2_IPCHI2_OWNPV)))",
           #    "L12_minIPCHI2",
           #    "",
           #    "F",
           #)
           #dataloader.AddVariable("log10(abs(Jpsi_IPCHI2_OWNPV))", "log_Jpsi_IPCHI2", "", "F")
            # define kfold cut
            # sigCut_kf = r.TCut("(eventNumber%10)!={}".format(i))
            # bkgCut_kf = r.TCut("(eventNumber%10)!={}".format(i))

            # learning rate is scaled for filenames
            learning_rate *= 0.01

            # get number of entries in signal and background trees
            n_signal = tree_s.GetEntries()
            n_background = tree_b.GetEntries()

            # define weights
            mysigw = 1.0
            mybkgw = 1.0
            #n_signal / n_background  # z. B. 0.1

            # define signal and background trees
            dataloader.AddTree(tree_s, "Signal", mysigw)  # , sigCut_kf
            dataloader.AddTree(tree_b, "Background", mybkgw)  # bkgCut_kf

            # Add correction weights to MC signal
            #dataloader.SetSignalWeightExpression("total_weights")  # total_weights is the name of the weight variable in the signal tree

            # define additional cuts
            # dataloader.SetSignalWeightExpression("PT_weight*LDPAsy_weight*PIDK_map*(L0MuonWeight_Run2_v0 * Jpsi_L0MuonDecision_TOS + L0DiMuonWeight_Run2_v0 * Jpsi_L0DiMuonDecision_TOS - L0MuonWeight_Run2_v0 * Jpsi_L0MuonDecision_TOS * L0DiMuonWeight_Run2_v0 * Jpsi_L0DiMuonDecision_TOS)") #MC reweighting
            #finite_cut = " && ".join([f"TMath::Finite({feature})" for feature in args.features])
            #sigCut = r.TCut(finite_cut)
            #bkgCut = r.TCut(finite_cut)
            sigCut = r.TCut("1")
            bkgCut = r.TCut("1")
            # set options for trainings
            n_train_signal = int(0.7 * tree_s.GetEntries())
            n_train_bkg = int(0.7 * tree_b.GetEntries())
            dataloader.PrepareTrainingAndTestTree(
                sigCut,
                bkgCut,
                ":".join(
                    [
                        f"nTrain_Signal={n_train_signal}",
                        f"nTrain_Background={n_train_bkg}",
                        "SplitMode=Random",
                        "SplitSeed=68",
                        "NormMode=NumEvents",
                        "!V"
                    ]
                ),
            )
            factory.BookMethod(
                dataloader,
                r.TMVA.Types.kBDT,
                "BDTG",
                ":".join(
                    [
                        "!H",
                        "V",
                        "NTrees={}".format(N_TREES),
                        "NegWeightTreatment=PairNegWeightsGlobal",
                        "MaxDepth={}".format(depth),
                        "BoostType=Grad",
                        "Shrinkage={}".format(learning_rate),
                        "UseBaggedBoost",
                        "BaggedSampleFraction=0.55",
                        "SeparationType=CrossEntropy",
                    ]
                ),
            )
            factory.TrainAllMethods()
            factory.TestAllMethods()
            factory.EvaluateAllMethods()

# r.TMVA.TMVAGui("TMVA_N{}_Sh{}_D{}.root".format(N_TREES, learning_rate, depth))#_kfold{}, i

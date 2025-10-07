import uproot3 as uproot
import numpy as np
from Sampling import *
from lib.src.TCFit.parameter import Parameter
from lib.src.TCFit.formula import Formula
from lib.src.TCFit.basic_pdfs import Gaussian, Exponential, CrystalBall
from lib.src.TCFit.model import AddPDF
from lib.src.TCFit.model import Model
from lib.src.TCFit.fit import LikelihoodFit
from lib.src.TCFit.plot import Plot
import matplotlib.pyplot as plt

prepare = True
sampleSize = 10000
if prepare:

    bin_dict = {"P": [5000., 7500., 10000., 15000., 20000., 30000., 40000., 60000., 100000., 200000.],
                "ETA": [1.9, 2.4, 2.8, 3.0, 3.2, 3.6, 4.0, 4.5, 4.9],
                "nSPDHits": [0, 100, 200, 300, 400, 450],
                "nPVs": [0.5, 1.5, 2.5, 3.5, 4.5, 5.5]}
    obs = "J_psi_1S_M"
    var_list = [i for i in sorted(bin_dict.keys())]
    var_dict_2D = {
        "P-ETA": ["P", "ETA"],
    }
    var_list_glob = [
        "totCandidates",
        "nTracks",
        "nLongTracks",
        "nTTTracks",
        "nVeloTracks",
        "nPVs",
        "nSPDHits",
        "nVeloClusters",
        "nITClusters",
        "nOTClusters",
        "nTTClusters",
        "nRich1Hits",
        "nRich2Hits",
        "nMuonTracks",
        "matched",
        "Polarity",
        "runNumber"
    ]

    cuts = "Mother_ETA > 0"

    method = "Long"

    infile = uproot.open("../TrackingEfficiencyTool/trackEffTuple_MC_2015_25ns_Long_Sim09b.root")
    treePlus = infile["TrackEffTreePlus" + method]
    treeMinus = infile["TrackEffTreeMinus" + method]
    dfs_orig = [None] * 2
    dfs_orig[1] = treePlus.pandas.df()
    dfs_orig[0] = treeMinus.pandas.df()
    dfs = {}
    actual_entries = np.shape(dfs_orig[0])[0] + np.shape(dfs_orig[1])[0]

    print("Sample Size before generic cuts: {}".format(actual_entries))

    dfs[-1] = dfs_orig[0][select(dfs_orig[0], cuts)]
    dfs[1] = dfs_orig[1][select(dfs_orig[1], cuts)]

    sampleSize = np.shape(dfs[-1])[0] + np.shape(dfs[1])[0]
    print("Sample size after generic cuts: " + str(sampleSize))
    datasets = [None] * 3
    m = {}
    m[-1] = matching(dfs[-1])
    m[1] = matching(dfs[1])

    datasets[0] = np.concatenate([dfs[-1][m[-1]][obs].to_numpy(), dfs[1][m[1]][obs].to_numpy()])
    datasets[1] = np.concatenate([dfs[-1][~m[-1]][obs].to_numpy(), dfs[1][~m[1]][obs].to_numpy()])
    datasets[2] = np.concatenate([dfs[-1][obs].to_numpy(), dfs[1][obs].to_numpy()])

    var_mult = []
    var_rem = []
    for i in var_list:
        if i in var_list_glob:
            var_mult.append(i)
        else:
            var_rem.append(i)

    for key, varl in var_dict_2D.items():
        for var in varl:
            if var not in var_rem:
                var_rem.append(var)

    list_of_vars = var_rem
    list_of_vars.extend(var_mult)
    list_of_vars.extend(var_dict_2D.keys())

    ds_dict = make_datasets(dfs, list_of_vars, bin_dict, obs)

mass = Parameter("J_psi_1S_M", limits=(2900., 3300.), unit="MeV/c^{2}")
meanCB = Parameter("meanCB", 3100.0, limits=(3090.0, 3150.0))
meanP = Parameter("meanP", 3100.0, limits=(3090.0, 3150.0))
meanF = Parameter("meanF", 3100.0, limits=(3090.0, 3150.0))
sigma = Parameter("sigma", 30.0, limits=(20.0, 60.0))
sigmaP = Parameter("sigmaP", 30.0, limits=(20.0, 40.0))
sigmaF = Parameter("sigmaF", 30.0, limits=(20.0, 40.0))
sigma2 = Parameter("sigma2", 90.0, limits=(15.0, 300.0))
fracCB = Parameter("fracCB", 0.3, limits=(0.3, 0.7))
alpha = Parameter("alpha", 1.0, limits=(0., 15.0))
n = Parameter("n", 1.2, limits=(1.1, 15.0))
tau = Parameter("tau", 0.001, limits=(0, 0.1))
tauF = Parameter("tauF", 0.0, limits=(-0.01, 0.01))
signal_yield = Parameter("signal_yield", 0.5 * sampleSize, limits = (1, sampleSize))
background_yield = Parameter("background_yield", 0.01 * sampleSize, limits = (0.1, sampleSize))

myCB = CrystalBall("myCB", mass, meanCB, sigma, alpha, n)
myCB2 = CrystalBall("myCB2", mass, meanCB, sigma2, alpha, n)
myCBsum = AddPDF("myCBsum", [myCB, myCB2], coefficients=[fracCB])

myGauss = Gaussian("myGaussian", mass, meanCB, sigma)
myGauss2 = Gaussian("myGaussian2", mass, meanCB, sigma2)
myGaussSum = AddPDF("myGaussianSum", [myGauss, myGauss2], coefficients=[fracCB])

# some definitions for simultaneous fit
sigma2F = Parameter("sigma2F", 70.0, limits = (15.0, 300.0))
fracCBF = Parameter("fracCBF", 0.5, limits = (0., 1.))
myExpoP = Exponential("myExpoP", mass, tau)
myCBP = CrystalBall("myCBP", mass, meanP, sigma, alpha, n)
myCB2P = CrystalBall("myCB2P", mass, meanCB, sigma2, alpha, n)
myCBSumP = AddPDF("myCBsumP", [myCBP, myCB2P], coefficients=[fracCB])
myExpoF = Exponential("myExpoF", mass, tau)
myCBF = CrystalBall("myCBF", mass, meanF, sigma, alpha, n)
myCB2F = CrystalBall("myCB2F", mass, meanCB, sigma2F, alpha, n)
myCBSumF = AddPDF("myCBsumF", [myCBF, myCB2F], coefficients=[fracCBF])

myGaussP = Gaussian("myGaussianP", mass, meanP, sigmaP)
myGauss2P = Gaussian("myGaussian2P", mass, meanCB, sigma2)
myGaussSumP = AddPDF("myGaussianSumP", [myGaussP, myGauss2P], coefficients=[fracCB])
signal_yieldP = Parameter("signal_yieldP", 0.5 * sampleSize, limits = (0.1, sampleSize))
background_yieldP = Parameter("background_yieldP", 0.00001 * sampleSize, limits = (0.0, sampleSize ))

myGaussF = Gaussian("myGaussianF", mass, meanP, sigmaF)
myGauss2F = Gaussian("myGaussian2F", mass, meanCB, sigma2F)
myGaussSumF = AddPDF("myGaussianSumF", [myGaussF, myGauss2F], coefficients=[fracCBF])
signal_yieldF = Parameter("signal_yieldF", 0.5 * sampleSize, limits =(0.1, sampleSize))
background_yieldF = Parameter("background_yieldF", 0.00001 * sampleSize, limits= (0.0, sampleSize))

efficiency_sig = Parameter("efficiency_sig", 0.9, limits=(0.,1.))
# efficiency_bkg  = Parameter("efficiency_bkg",0.5, low = 0., high = 1.)
f_pass = Formula("f_pass", "@0*@1", [signal_yield, efficiency_sig])
f_fail = Formula("f_fail", "@0*(1.-@1)", [signal_yield, efficiency_sig])
# f_pass_bkg      = Formula("f_pass_bkg", "@0*@1", [background_yield,efficiency_bkg])
# f_fail_bkg      = Formula("f_fail_bkg", "@0*(1.-@1)", [background_yield,efficiency_bkg])

# totalShapeP = Model("totalShapeP",
#                    [myGaussSumP, myExpoP],
#                    [f_pass, background_yieldP])

# print(totalShapeP)
# totalShapeF = Model("totalShapeF",
#                    [myGaussSumF, myExpoF],
# 3                    [f_fail, background_yieldF])
print("Declaring total shape P")
totalShapeP = AddPDF("totalShapeP",
                    [myCBSumP, myExpoP],
                    [f_pass, background_yieldP])



print("test --------------------- ")
for par in totalShapeP.parameters:
    print(par)

for coef in totalShapeP.coefficients:
    print(coef)


print(totalShapeP)
totalShapeF = AddPDF("totalShapeF",
                    [myCBSumF, myExpoF],
                    [f_fail, background_yieldF])


totalShape = Model("sim fit", pdfs=[totalShapeP,totalShapeF])
print(ds_dict)

if prepare:
    data_set = []
    data_set.append(np.concatenate([ds_dict["ETA"][+1]["matched"][0], ds_dict["ETA"][-1]["matched"][0]]))
    data_set.append(np.concatenate([ds_dict["ETA"][+1]["failed"][0], ds_dict["ETA"][-1]["failed"][0]]))
    signal_yield.high = (len(data_set[0]) + len(data_set[1]))
    signal_yield.value = len(data_set[0]) + len(data_set[1]) * 0.9
    signal_yield.low = (len(data_set[0]) + len(data_set[1])) * 0.2
    efficiency_sig.value = len(data_set[0]) * 1. / signal_yield.value * 1.

    signal_yieldP.high = len(data_set[0])
    signal_yieldP.value = signal_yieldP.high / 1.1
    signal_yieldF.high = len(data_set[1])
    signal_yieldF.value = signal_yieldF.high / 1.1

    background_yieldP.value = 0.0
    background_yieldP.high = len(data_set[0])
    background_yieldF.value = 0.0
    background_yieldF.high = len(data_set[1])

    models = [totalShapeP, totalShapeF]

    print(data_set)
    print(len(data_set))
    print(totalShape)
    print(len(totalShape.pdfs))
    fit = LikelihoodFit(data=data_set, model=totalShape)
    minuit = fit.minimise(hesse=False, minos=False)
    minuit = fit.minimise(hesse=True, minos=True)
    # minuit = fit.minimize(hesse=True, minos=True)
    # minuit = fit.minimize(hesse=False, minos=True)

    print("is valid :: ", minuit.fmin.is_valid)
    for par in minuit.params:
        print(par.name, par.value, par.error, par.merror)

    nbins = 50
    plot = Plot(minuit=minuit, model=totalShape, data=data_set, limit=mass.limits, n_bins=nbins)
    plot.plot()
    plot.save("total")

    binwidth = ((mass.limits[1] - mass.limits[0]) / nbins)
    low, high = mass.limits

    fig, (ax1, ax2) = plt.subplots(1, 2)
    fig.suptitle('Simultaneous fit')

    ax1.hist(data_set[0], bins=nbins, range=[low, high], histtype='step', edgecolor='r', linewidth=3, label="pass")
    ax2.hist(data_set[1], bins=nbins, range=[low, high], histtype='step', edgecolor='b', linewidth=3, label="fail")

    x = np.linspace(low, high, 1000)
    yP = totalShapeP.evaluate(x) * binwidth
    yCBP = signal_yield.value * efficiency_sig.value * myCBP.evaluate(x) * binwidth
    yBP = background_yieldP.value * myExpoP.evaluate(x) * binwidth
    ax1.plot(x, yCBP, label='fit Passed')
    ax1.plot(x, yBP, label='background', color='g', linestyle='--')
    yF = totalShapeF.evaluate(x) * binwidth
    ax2.plot(x, yF, label='fit Failed')
    ax1.plot()
    ax2.plot()
    fig.savefig("prova_new.pdf")

    print(tau)
    print(f_pass.value)
    print(f_fail.value)

    print((len(data_set[0]) + len(data_set[1])))

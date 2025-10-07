# to run this script: lb-conda pidgen, next: python pid_correct.py

import argparse
import numpy as np

# import ROOT as r

##START OF CONFIG
# list of MC ntuples, format: (inputfile, outputfile, dataset) -> note that the outputfile should be without ".root" extension

parser = argparse.ArgumentParser(description="Input tuple and dataset")
parser.add_argument("--input", action="store", dest="input", type=str, required=True)
parser.add_argument("--output", action="store", dest="output", type=str, required=True)
parser.add_argument("--dataset", action="store", dest="dataset", type=str, required=True)
args = parser.parse_args()

dict = {
    "MD2015": "MagDown_2015",
    "MD2016": "MagDown_2016",
    "MD2017": "MagDown_2017",
    "MD2018": "MagDown_2018",
    "MU2015": "MagUp_2015",
    "MU2016": "MagUp_2016",
    "MU2017": "MagUp_2017",
    "MU2018": "MagUp_2018",
}

files = [(f"{args.input}", f"{args.output}"[:-5], f"{dict[args.dataset]}")]

# name of input tree
input_tree = "DecayTree"

# postfixes for Pt, Eta and Ntracks variables
ptvar = "PT"
etavar = "ETA"
ntrvar = "nTracks"
pvar = None
simversion = "Sim09"

friend = False  # if True, friend tree is written instead of copying whole tree

scale = (
    1,
    1,
    1.1,
)  # list of scale factors for each dimension of input data. in the example, nTracks is increased by 15%.

# list of (kernel, seeds) combinations for raw template smearing
kernels = [("default", [0])]  # three 'bootstrapped' templates with seeds 1...3

# Configuration dictionary for resampling, in the form {particle_name}:{pidvars}
# For each {particle_name}, {pidvars} is a dictionary in the form {ntuple_variable}:({sample}, {PID_var}, {kernels}),
#   where
#     {ntuple_variable} is the name of the corresponding ntuple PID variable without the particle name part
#                       (e.g. for "pi_PIDK" branch of particle "pi" it should be just "PIDK");
#     {sample} is the name of the calibration sample
#     {PID_var} is the string describing the PID variable template.
#     {kernels} is the list of kernels for template smearing (see example above)

config = {
    "L2": {"ProbNNmu": ("ProbNNmu", "mu_Jpsi2mumunopt_isMuon", "MC15TuneV1_ProbNNmu", kernels)},
    "L1": {"ProbNNmu": ("ProbNNmu", "mu_Jpsi2mumunopt_isMuon", "MC15TuneV1_ProbNNmu", kernels)},
}

##END OF CONFIG
output_tree = input_tree
treename = input_tree


from pidgen2.correct import correct

for input_file, output_file, dataset in files:

    input_location = f"{input_file}:{input_tree}"

    for track, subst in config.items():
        for var, (pidvar, sample, calibvar, kernel_list) in subst.items():

            # Create the list of input branches, depending on whether Eta or P variable is available
            if pvar is None:
                branches = f"{track}_{pidvar}:{track}_{ptvar}:{track}_{etavar}:{ntrvar}"
                eta_from_p = False
            else:
                branches = f"{track}_{pidvar}:{track}_{ptvar}:{track}_{pvar}:{ntrvar}"
                eta_from_p = True

            if friend:
                output_root_file = f"{output_file}_{track}_{var}.root"
            else:
                output_root_file = f"{output_file}.root"

            # Run resampling of a single variable in a single file
            correct(
                input=input_location,  # Input tuple
                simversion=simversion,  # Simulation version of the input MC file
                sample=sample,  # Calibration sample (e.g. mu_Jpsi2mumu)
                dataset=dataset,  # Dataset (e.g. "MagDown_2018")
                variable=calibvar,  # Calibration variable (e.g. "MC15TuneV1_ProbNNmu")
                branches=branches,  # List of resampling branches (typically, corresponding to Pt, Eta and Ntracks, e.g. "pt:eta:ntr")
                output=output_root_file,  # Output ROOT file name
                outtree=output_tree,  # Output tree name
                plot=True,  # If template needs to be created from scratch, produce control plots
                pidcorr=f"{track}_{var}_pidcorr",  # Name of the corrected PID branch
                stat=f"{track}_{var}_pidstat",  # Name of output branch with calibration statistics for each resampled event
                mcstat=f"{track}_{var}_pidmcstat",  # Name of output branch with MC statistics for each resampled event
                kernels=kernel_list,  # List of kernels and template seeds
                verbose=False,  # Print debugging information
                eta_from_p=eta_from_p,  # If eta needs to be calculated from p and pt
                friend=friend,  # If True, write output to friend trees
                library="ak",  # Library to handle ROOT files with uproot if friend=False, can be
                # "ak" (Awkward array), "np" (numpy) or "pd" (Pandas)
                step_size=50000,  # Chunk size when writing files with uproot if friend=False
                nan=-1000.0,  # Numerical value to substitute NaN, for regions w/o calibration data
                scale=scale,  # Scale factors for input data
                local_storage="/ceph/users/pmachnik/templates",  # Directory for local template storage, used when the template is not available in
                # the global storage on EOS.
                local_mc_storage="/ceph/users/pmachnik/templates",  # Directory for local template storage, used when the template is not available in
                # the global storage on EOS.
            )

            if not friend:
                # All subsequent calls use output file as input
                input_location = f"{output_root_file}:{output_tree}"

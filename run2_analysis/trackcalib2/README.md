# TrackCalib2

contacts: Maurice Pierre Morgenthaler (maurice.pierre.morgenthaler _at_ cern.ch) Giulia Frau (giulia.frau _at_ cern.ch) and Flavio Archilli (flavio.archilli _at_ cern.ch)

TrackCalib2, the successor of the the beloved [TrackCalib twiki page](https://twiki.cern.ch/twiki/bin/view/LHCb/TrackCalib), is a tool to allow custom variables, binning, etc. for Monte-Carlo tracking efficiency correction evaluation. TrackCalib2 package is now Python3 compatible and it doesn't have any dependency on other packages (e.g. ROOT or Urania). The principle of use is almost identical to the original TrackCalib.


If the default tracking efficiency correction tables, found at the [TrackCalib twiki page](https://twiki.cern.ch/twiki/bin/view/LHCb/TrackCalib "TrackCalib twiki page"), are not suitable for your analysis, this tool is the correct one for you!

## Prerequisites

All the necessary packages are listed together with the used versions in [requirements.yaml](https://gitlab.cern.ch/farchill/trackcalib2/-/blob/master/requirements.yaml). To install them, set up a conda environment called `trackcalib2` by running the following command:

```
conda env create -f requirements.yaml
```

Then, the library with the fitter needs to be initialized and downloaded:
```
git submodule init
git submodule update
```

Now you should be ready to go! 
To use the conda environment append `conda run -n trackcalib2` before the python command.

## Help
The manual with all the parameter details, options and default values are printed directly in the terminal in order to keep the manual as up-to-date as possible. Running with the -h or -help option is not going to print the details, instead, run
```
conda run -n trackcalib2 python trackcalib.py manual
```

## Running TrackCalib2
TrackCalib2 is divided into three sub-tasks that have to be invoked separately and one and only one is mandatory. The tasks are:

* `prepare`: only runs Prepare script, reducing the dataset to specified variables, weighting MC in multiplicity, and applying additional matching criteria
* `fit`: creates the efficiency fits in a given set of variables, mode, and method with an earlier created tuple
* `plot`: Plots the efficiencies from the fit, and creates data/MC correction tables in all 2D variables

Before running the whole thing, always test all the three steps with the -verbose and -max_entries smallNumber option! Otherwise you risk a lot of wasted computing time: not good for you, nor the collaboration, nor the planet.

### Prepare 

This step takes the tuples stored at EOS and creates smaller ROOT file(s) with the relevant information. If you do not need special preselection, skip this step and use the -official option at the fit part. The ROOT files are handled using uproot3.

An example, that creates the smaller ROOT files for MC, long Method, simulation version Sim09b for the year 2016:

```
conda run -n trackcalib2 python trackcalib.py prepare -year "2016_25ns" -mode "MC" -method "Long" -sim_ver "Sim09b" 
```

### Fit

Takes the tuple(s) created by the prepare step and performs the fits in the J/psi mass in the desired bins. The output is stored as a pickle (.pkl) file. If you wish to test eg different mass shapes of J/psi, use the option -subfolder to prevent overwriting the existing pickle file. This option does not store the fit output in ./results/year/ but in ./results/year/subfolder.

Always check the automatically created output and warnings files + the Jpsi mass fit plots created in ./results/year/subfolder/plots!

An example of fitting the Velo method on 2012 Data with simultaneous fit to the failed and matched J/psi:
```
conda run -n trackcalib2 python trackcalib.py fit -year "2012_50ns_strip" -v -sim_fit -method "Velo" -mode "Data" 
```

### Plot

Running the plot task creates pretty plots from your fit results. It takes either the default list of pickle files (ie Long, Velo, T methods in Data and MC for given year) and creates all plots automatically, or it takes a list of pickle files to be plotted. An example would be the direct comparison of 2011 and 2012 efficiencies:

```
 conda run -n trackcalib2 python trackcalib.py plot -year "2012_50ns_strip" -file_list "results/2011_50ns_strip/trackEff_Data_Long.pkl,results/2012_50ns_strip/trackEff_Data_Long.pkl" 
```

This way of passing the files is rather lengthy, however this ensures the most possible versatility of the code. If you pass files with Velo and T methods, the Combined method is automatically computed and saved as a pickle file, but only if the sample is one mode and one year! Similarly, Final method is calculated and saved, if possibly. Lastly, if you provide same samples in Data and MC, the ratios are calculated and saved.



## Useful links:
* [TrackCalib2 report with a short how-to](https://indico.cern.ch/event/1164387/contributions/4889759/attachments/2453815/4205270/2022_06_01_TrackCalib2.pdf)
* [TrackEff twiki](https://twiki.cern.ch/twiki/bin/viewauth/LHCbInternal/LHCbTrackingEfficiencies)
* [Outdated TrackCalib twiki](https://twiki.cern.ch/twiki/bin/view/LHCb/TrackCalib)
* [Run 1 TrackEff paper](http://iopscience.iop.org/1748-0221/10/02/P02007/)						
* [Chapter 10 in Renata's thesis](https://archiv.ub.uni-heidelberg.de/volltextserver/30967/)	

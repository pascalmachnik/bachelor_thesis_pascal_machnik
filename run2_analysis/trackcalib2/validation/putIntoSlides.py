from __future__ import print_function
import numpy as np 

from putIntoSlides_Utils import *

#Path to your results folder
FIT_PLOT_FOLDER  = "/home/renata/lxplus/afs_work/TrackCalib2/trackcalib2/results/"

#Path to your plots folder
EFF_PLOT_FOLDER  = "./figures/"

#1D variables with corresponding number of bins
DEF_VARS =  {"P": 5,
                "ETA": 2,
                "nSPDHits": 5,
                "nPVs": 5
            }

#1D variables with corresponding number of bins for fine binning
DEF_VARS_FINE = {"P": 9,
                "ETA": 8,
                "nSPDHits": 5,
                "nPVs": 5
                }

#2D variables with corresponding number of bins
DEF_VARS_2D = {"P-ETA" : [5,2]}

#Basic methods
METHODS = ["Long", "Velo", "T"]

#Methods including combined and final
METHODS_ALL = [ "Velo", "T", "Long", "Combined", "Final"]

#Dictionary with the years and year tags, latex is not happy with _, hence the tags
YEARS = {
    "2011_50ns_strip": "2011strip",
    "2012_50ns_strip": "2012strip",
    "2016_25ns_strip": "2016strip",
    "2017_25ns_strip": "2017strip",
    "2015_50ns":       "2015with50ns",
    "2015_25ns":       "2015with25ns",
    "2016_25ns":       "2016",
    "2017_25ns":       "2017",
    "2018_25ns":       "2018"  
}


def getFitPlotName(mode, method, var = "", bin = -1, tag = ""):          
    '''
        Returns the fit plot name. 
    '''
    if (tag != ""): tag = "_"+tag
    if type(var) == list: #2D binning
        if (len(bin) != len(var)):
            print(f"Number of bins {len(bin)} does not agee with the number of variables {len(var)}!")
        path =  f"{mode}_{method}_{var[0]}_bin_{bin[0]}_{var[1]}_bin_{bin[1]}{tag}"     
    else: #1D binning             
        if (var == ""): #Full dataset, no binning in any var
            path = f"{mode}_{method}_full_sample{tag}"
        else:
            if (bin==-1):
                print(f"Wrong variable {var} or bin {bin}!")
            path = f"{mode}_{method}_{var}_{bin}{tag}"    
    return path

def massFits1D(folder, year, mode, method, var, nBins, tag): 
    '''
        Dumps 1D mass fits into slides
    '''
    path = getFitPlotName(mode, method, var, 0, tag)
    path = path.replace("_0_","___LOOP___")

    listOfEntries = list(map(str,np.arange(0,nBins[0],1))) 
    defaultTitle = f"{year} {mode} {method} {var}" 

    outputName = f"texFiles/{year}/{year}_{mode}_{method}_{var}_{tag}" 
    figsPerSlide= 4 
    figsPerRow= 2
    make_figure_slides(folder, path, listOfEntries, defaultTitle, 
                        figsPerSlide, figsPerRow, ".pdf", outputName+".tex")


def massFits2D(folder, year, mode, method, var, nBins,tag): 
    '''
        Dumps 2D mass fits into slides
    '''
    path = getFitPlotName(mode, method, var, [0,0], tag)   
    path = path.replace(f"_{var[0]}_bin_0_{var[1]}_bin_0_","__LOOP__")


    listOfEntries = []
    for bin_x in range(0,nBins[0]):
        for bin_y in range(0,nBins[1]):
            listOfEntries.append(f"_{var[0]}_bin_{bin_x}_{var[1]}_bin_{bin_y}_")    
    defaultTitle = f"{year} {mode} {method} {var[0]} {var[1]}" 
    outputName = f"texFiles/{year}/{year}_{mode}_{method}_{var[0]}_{var[1]}_{tag}" 
    figsPerSlide= 4 
    figsPerRow= 2

    make_figure_slides(folder, path, listOfEntries, defaultTitle, 
                        figsPerSlide, figsPerRow, ".pdf", outputName+".tex")    

def efficiencies_MCvsData(year, yearTag, subfolder, var):
    '''
        Dumps 1D efficiency plots into slides
    '''
    folder = f"{EFF_PLOT_FOLDER}plots/{year}/{subfolder}/"
    plotName = f"trackEff_Data_MC_{year}___LOOP___{var}"

    listOfEntries = METHODS_ALL

    defaultTitle = f"Efficiency {yearTag} {var}" 

    outputName = f"texFiles/{yearTag}/{yearTag}_MCvsData_{var}" 
    figsPerSlide= 6
    figsPerRow= 3
    make_figure_slides(folder, plotName, listOfEntries, defaultTitle, 
                        figsPerSlide, figsPerRow, ".png", outputName+".tex")
    return


def ratio(year, yearTag, subfolder):
    '''
        Dumps efficiency ratio plots into slides
    '''
    folder = f"./figures/plots/{year}/{subfolder}/"
    plotName = f"trackEff_Ratio__{year}___LOOP___P-ETA"

    listOfEntries = METHODS_ALL

    defaultTitle = f"Ratio {yearTag}" 

    outputName = f"texFiles/{yearTag}/{yearTag}_Ratios" 
    figsPerSlide= 6
    figsPerRow= 3
    make_figure_slides(folder, plotName, listOfEntries, defaultTitle, 
                        figsPerSlide, figsPerRow, ".png", outputName+".tex")
    return


    
#Don't forget to use  \epstopdfsetup{outdir=./figures/out/2011/}  if you are loading the plots directly using ssh from lxplus

def massFits(subfolder, mode, simVer = ""):
    for year, yearTag in YEARS.items():
        for method in METHODS:
            for var, nBins in DEF_VARS.items():
                massFits1D(FIT_PLOT_FOLDER+f"{year}/{subfolder}/{simVer}/plots/",yearTag,mode,method,var,[nBins], "sim")
            for var, nBins in DEF_VARS_2D.items():
                vars = var.split("-")
                massFits2D(FIT_PLOT_FOLDER+f"{year}/{subfolder}/{simVer}/plots/",yearTag,mode,method,vars,nBins, "sim")
    return

def efficiencies(subfolder):
    for year, yearTag in YEARS.items():
        for var in DEF_VARS.keys():
            efficiencies_MCvsData(year,yearTag,subfolder,var)

def ratios(subfolder):
    for year, yearTag in YEARS.items():
        ratio(year,yearTag,subfolder)


massFits("","Data")
massFits("","MC","Sim09h")
efficiencies("")
ratios("")
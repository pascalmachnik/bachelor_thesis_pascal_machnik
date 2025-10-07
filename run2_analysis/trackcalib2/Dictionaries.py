import re
from Utilities import OutputColoring as OC
from typing import List
from collections import namedtuple

# code version
version = 0.6

# contacts and their emails
contacts = {
    "Flavio Archilli": "flavio.archilli@cern.ch",
    "Raphael Skuza": "raphael.skuza@cern.ch"
}

# authors
authors = ["F. Archilli", "M. Kolpin", "R. Kopecna", "R. Skuza"]

# methods available
methods = ["Long", "T", "Velo"]
modes = ["MC", "Data"]
polarities = ["MagUp", "MagDown"]

# available and default MC versions available for each datataking year
simVer = namedtuple("simVer", "sim_list sim_default")

# available years, available and  default sim versions
year_dict = {
    "2011_50ns_strip": simVer(["Sim09h"], "Sim09h"),
    "2012_50ns_strip": simVer(["Sim09h"], "Sim09h"),
    "2016_25ns_strip": simVer(["Sim09h"], "Sim09h"),
    "2017_25ns_strip": simVer(["Sim09h"], "Sim09h"),
    "2015_50ns": simVer(["Sim08h", "Sim09b"], "Sim09b"),
    "2015_25ns": simVer(["Sim09a", "Sim09b"], "Sim09b"),
    "2016_25ns": simVer(["Sim09a", "Sim09b", "Sim09d"], "Sim09b"),
    "2017_25ns": simVer(["Sim09h"], "Sim09h"),
    "2018_25ns": simVer(["Sim09h"], "Sim09h")
}

# In 2016 and 2017, the T method is missing
# There were two different bugs
# Use this list
year_missing_T_method = ["2016_25ns", "2017_25ns"]

### Fit dictionaries ###
# fit configuration parameters
fit_config = {
    "nbins": 30,
    "obs": {
        "name": "J_psi_1S_M",
        "bounds": {
            "Long": (2625.0, 3575.0),
            "Velo": (2925.0, 3275.0),
            "T": (2625.0, 3575.0),
            "T_Run1": (2700.0, 3500.0)
        }
    }
}

#Binning configuration
bin_dict_default = {
    "P": [5000., 10000., 20000., 40000., 100000., 200000.],
    "ETA": [1.9, 3.2, 4.9],
    "nSPDHits": [0, 100, 200, 300, 400, 450],
    "nPVs": [0.5, 1.5, 2.5, 3.5, 4.5, 5.5]
}

bin_dict_fine = {
    "P": [
        5000., 7500., 10000., 15000., 20000., 30000., 40000., 60000., 100000.,
        200000.
    ],
    "ETA": [1.9, 2.4, 2.8, 3.0, 3.2, 3.6, 4.0, 4.5, 4.9],
    "nSPDHits": [0, 100, 200, 300, 400, 450],
    "nPVs": [0.5, 1.5, 2.5, 3.5, 4.5, 5.5]
}

# default set of 2D variable combinations
var_dict_2D_default = {"P-ETA": ["P", "ETA"]}

#default variables lists to keep on tuple (split in global/particle vars)
var_list_part = [
    "P",
    "PT",
    "ETA",
    "M",
    "IPCHI2_OWNPV",
    "OWNPV_X",
    "OWNPV_Y",
    "OWNPV_Z",
    "IP_OWNPV",
    "TRACK_CHI2NDOF",  #no jpsi
    "TRACK_CHI2",  #no jpsi
    #    "BKGCAT",
    #    "Reconstructed",
    "Assoc",  #no jpsi
    "Overlap_TTFraction",  #no jpsi
    "Overlap_ITFraction",  #no jpsi
    "Overlap_OTFraction",  #no jpsi
    "Overlap_TFraction",  #no jpsi
    "Overlap_VeloFraction",  #no jpsi
    "Overlap_MuonFraction",  #no jpsi
    "Matched_GhostProb",  #no jpsi
    #no mu mc#    "TRACK_nVeloHits",#no jpsi   #no mu mc
    #no mu mc#    "TRACK_nVeloRHits",#no jpsi
    #no mu mc#    "TRACK_nVeloPhiHits",#no jpsi
    #no mu mc#    "TRACK_nTTHits",#no jpsi
    #no mu mc#    "TRACK_nITHits",#no jpsi
    #no mu mc#    "TRACK_nOTHits",#no jpsi
    #no mu mc#    "TRACK_nMuHits",
    "PIDmu"  #no jpsi
]

var_list_mother = [
    "P",
    "PT",
    "ETA",
    "M",
    "ENDVERTEX_CHI2",  #no mu
    "ENDVERTEX_X",  #no mu
    "ENDVERTEX_Y",  #no mu
    "ENDVERTEX_Z",  #no mu
    "IPCHI2_OWNPV",
    "OWNPV_X",
    "OWNPV_Y",
    "OWNPV_Z",
    "IP_OWNPV",
]

#If any other variable is needed by the user, it needs to be hardcoded here.
var_list_glob = [
    "totCandidates", "nTracks", "nLongTracks", "nTTracks", "nVeloTracks",
    "nPVs", "nSPDHits", "nVeloClusters", "nITClusters", "nOTClusters",
    "nTTClusters", "nRich1Hits", "nRich2Hits", "nMuonTracks", "Polarity",
    "runNumber"
    #    "matched"
]

var_list_L0 = ["L0MuonDecision_TOS", "L0MuonDecision_TIS"]

var_list_Hlt1 = [
    "Hlt1TrackMuonDecision_TOS",
    "Hlt1TrackMuonDecision_TIS",
    "Hlt1TrackAllL0Decision_TOS",
    "Hlt1TrackAllL0Decision_TIS",
    "Hlt1TrackMVADecision_TOS",
    "Hlt1TrackMVADecision_TIS"  #,
    #    "Hlt1TrackMuonMVADecision_TOS" ,
    #    "Hlt1TrackMuonMVADecision_TIS"
]

#Default uncertainty and default value used both for efficiencies and ratios
default_error = 0.05
default_value = 1.0

#Variables that should have log scale in the 2D plot
logVars = ["P"]

latex_names_dict = {
    "J_psi_1S_M": r"M$_\mu^{+}\mu^{-}$",
    "P": r"$p$ [MeV/c]",
    "PT": r"$P_T$ [MeV/c]",
    "ETA": r"$\eta$",
    "M": r"$M$ [MeV/c$^2$]",
    "IPCHI2_OWNPV": r"IP $\chi^2$",
    "OWNPV_X": r"$X_{PV}$",
    "OWNPV_Y": r"$Y_{PV}$",
    "OWNPV_Z": r"$Z_{PV}$",
    "IP_OWNPV": r"$IP$",
    "TRACK_CHI2NDOF": r"Track $\chi^2$",
    "TRACK_CHI2": r"Track $\chi^2$",
    "BKGCAT": r"Bkg Category",
    "Matched_GhostProb": r"Ghost Prob Matched",
    "TRACK_nVeloHits": r"Track $\mathrm{N_{VELO~Hits}}$",
    "TRACK_nVeloRHits": r"Track $\mathrm{N_{VELO~R~Hits}}$",
    "TRACK_nVeloPhiHits": r"Track $\mathrm{N_{VELO~\phi~Hits}}$",
    "TRACK_nTTHits": r"Track $\mathrm{N_{TT~Hits}}$",
    "TRACK_nITHits": r"Track $\mathrm{N_{IT~Hits}}$",
    "TRACK_nOTHits": r"Track $\mathrm{N_{OT~Hits}}$",
    "TRACK_nMuHits": r"Track $\mathrm{N_{Muon~Hits}}$",
    "PIDmu": r"Muon PID",
    "totCandidates": r"$\mathrm{N_{Tot~Candidates}}$",
    "nTracks": r"$\mathrm{N_{All~Tracks}}$",
    "nLongTracks": r"$\mathrm{N_{Long~Tracks}}$",
    "nTTracks": r"$\mathrm{N_{T~Tracks}}$",
    "nVeloTracks": r"$\mathrm{N_{VELO~Tracks}}$",
    "nPVs": r"$\mathrm{N_{PV}}$",
    "nSPDHits": r"$\mathrm{N_{SPD~Hits}}$",
    "nVeloClusters": r"$\mathrm{N_{VELO~clusters}}$",
    "nITClusters": r"$\mathrm{N_{IT~clusters}}$",
    "nOTClusters": r"$\mathrm{N_{OT~clusters}}$",
    "nTTClusters": r"$\mathrm{N_{TT~clusters}}$",
    "nRich1Hits": r"$\mathrm{N_{RICH1~hits}}$",
    "nRich2Hits": r"$\mathrm{N_{RICH2~hits}}$",
    "nMuonTracks": r"$\mathrm{N_{Muon~Tracks}}$",
    "Polarity": r"Polarity",
}


# Set the nSPDHits based on the trigger response
def get_nSPDHits_list(year):
    '''
        Returns the list of nSPD hits cuts. These cuts are included in the Muon Trigger,
        hence they are different for Run I and Run II. 
        For Run I, it returns [0,600,900], for Run II [0,450,900]
    '''
    nSPDHitsList = [0]
    if (int(year.split("_")[0]) < 2014):
        nSPDHitsList.append(600)  #In Run I, the cut was at 600
    else:
        nSPDHitsList.append(450)  #In RunII, it was 450
    nSPDHitsList.append(900)

    return nSPDHitsList


def get_nSPDHits_cut(year, high_multi):
    '''
       Returns the stirng with the corresponding nSPDHits cut
    '''
    #get_nSPDHits_list returns [0,450 or 600, 900],
    # so cut list is [0, 450 or 600] or for high multi [450 or 600, 900]
    cut_list = get_nSPDHits_list(year)[int(high_multi):int(high_multi + 2)]
    return f"nSPDHits>{cut_list[0]},nSPDHits<{cut_list[1]}"


# define basic cuts; crucial for 2015, in 2016 all moved to the trigger
# (here a bit tighter cuts), except for totCandidates

# TODO add DLL cuts at some point as a cross-check,
# check for cuts for different years at some point

# WARNING for users: custom cuts overwrite these cuts!


def get_Tmethod_mass_cuts(year):
    if (year == "2011_50ns_strip" or year == "2012_50ns_strip"):
        return [2700, 3500]
    else:
        return [2600, 3600]


def get_default_cuts_dict(year, high_multi):
    nSPDHits_cut = get_nSPDHits_cut(year, high_multi)
    mass_T_cut = get_Tmethod_mass_cuts(year)
    default_cuts = {
        'Long': {
            'Probe': 'P>5000,PT>1000',
            'Tag': 'P>10000,PT>1300,TRACK_CHI2NDOF<2',
            'Mother': 'ENDVERTEX_CHI2<5,IP_OWNPV<0.8,M<3600,M>2600,PT>1000',
            'Global': 'totCandidates==1,' + nSPDHits_cut
        },
        'T': {
            'Probe': 'P>5000,PT>500,TRACK_CHI2NDOF<5',
            'Tag': 'P>7000,PT>500,TRACK_CHI2NDOF<5,IP_OWNPV>0.2',
            'Mother':
            f'ENDVERTEX_CHI2<5,M<{mass_T_cut[1]},M>{mass_T_cut[0]},PT>500',
            'Global': 'totCandidates==1,' + nSPDHits_cut
        },
        'Velo': {
            'Probe': 'P>5000,PT>700,TRACK_CHI2NDOF<10',
            'Tag': 'P>5000,PT>700,TRACK_CHI2NDOF<5,IP_OWNPV>0.5',
            'Mother': 'ENDVERTEX_CHI2<5,M<3300,M>2900',
            'Global': 'totCandidates==1,' + nSPDHits_cut
        }
    }
    return default_cuts


def get_very_tight_cuts(year, high_multi):
    nSPDHits_cut = get_nSPDHits_cut(year, high_multi)
    very_tight_cuts = {
        'Long': {
            'Probe': 'P>5000,PT>500',
            'Tag': 'P>10000,PT>1300,TRACK_CHI2NDOF<2,PIDmu>2',
            'Mother': 'ENDVERTEX_CHI2<5,IP_OWNPV<0.8,M<3600,M>2600,PT>1000',
            'Global': 'totCandidates==1,' + nSPDHits_cut
        },
        'T': {
            'Probe': 'P>5000,PT>500,TRACK_CHI2NDOF<5',
            'Tag': 'P>7000,PT>500,TRACK_CHI2NDOF<3,IP_OWNPV>0.2',
            'Mother': 'ENDVERTEX_CHI2<5,M<3600,M>2600,PT>500',
            'Global': 'totCandidates==1,' + nSPDHits_cut
        },
        'Velo': {
            'Probe': 'P>5000,PT>700,TRACK_CHI2NDOF<10',
            'Tag': 'P>5000,PT>700,TRACK_CHI2NDOF<5,IP_OWNPV>0.5,PIDmu>-2',
            'Mother': 'ENDVERTEX_CHI2<5,M<3300,M>2900',
            'Global': 'totCandidates==1,' + nSPDHits_cut
        }
    }
    return very_tight_cuts


# Define the default mathcing criteria that is used
default_matching_criteria = {'Long':\
                                 {'Probe': "Overlap_TTFraction > 0.4,\
                                            Overlap_MuonFraction > 0.4,\
                                            Assoc == True"
                                    },
                             'Velo':\
                                 {'Probe': 'Overlap_TTFraction > 0.4,\
                                            Overlap_TFraction > 0.4,\
                                            Assoc == True'
                                    },
                             'T':\
                                 {'Probe': 'Overlap_VeloFraction > 0.4,\
                                            Overlap_MuonFraction > 0.4, \
                                            Assoc == True'
                                    }
                             }


def get_var_list(list, type="part"):

    var_list = []
    if type.lower() == "global":
        for var in list:
            var_list.append(var)

    else:

        if type == "part":
            part_list = ['muplus', 'muminus']

        elif type == "mother":
            part_list = ['J_psi_1S']

        else:
            part_list = [type]

        for part in part_list:
            for var in list:
                var_list.append("%s_%s" % (part, var))

    return var_list


def get_default_var_list():

    var_list = []
    var_list += get_var_list(var_list_part)
    var_list += get_var_list(var_list_mother, "mother")
    var_list += get_var_list(var_list_L0, "mother")
    var_list += get_var_list(var_list_Hlt1)
    var_list += get_var_list(var_list_glob, "Global")

    return var_list


def translate_dict_cut(cut_dict, Probe="muplus"):

    Tag = "muminus"
    if Probe == "muminus":
        Tag = "muplus"
    Mother = "J_psi_1S"
    Global = ""

    meth_cuts = cut_dict

    cut = ''
    first = True
    for key, value in meth_cuts.items():

        l = value.replace(" ", "").split(',')

        if len(l) > 0:
            if key != "Global":
                tmp = "%s_%s" % (eval(key), l[0])
            else:
                tmp = "%s" % (l[0])
            if len(l) > 1:
                for item in l[1:]:
                    if key != "Global":
                        tmp += " && %s_%s" % (eval(key), item)
                    else:
                        tmp += " && %s" % (item)

            if first:
                cut += tmp
                first = False
            else:
                cut += " && %s" % tmp

    return cut


def translate_dictmeth_cut(cut_dict, method="Long", probe="muplus"):

    meth_cuts = cut_dict[method]
    return translate_dict_cut(meth_cuts, probe)


def translate_substring_dict(string):

    d = {}
    tmp = string.replace(" ", "").split(";")

    if len(tmp) > 0:
        for target in tmp:
            line = target.split(":")
            if len(line) > 1:
                d[line[0]] = line[1]

    return d


def translate_string_dict(string):

    d = {}
    tmp = re.split(r':\{|\};', string.replace(" ", ""))
    # remove trailing }
    tmp[-1] = tmp[-1].replace("}", "")
    
    if len(tmp) > 1:
        for key, value in zip(tmp[::2], tmp[1::2]):
            d[key] = translate_substring_dict(value)

    else:
        for method in methods:
            d[method] = translate_substring_dict(string)

    return d


def translate_dict_cutdict(dict, probe="muplus"):

    cuts = {}
    for method in methods:
        cuts[method] = translate_dictmeth_cut(dict, method, probe)

    return cuts


def get_matching_criteria(string, probe="muplus"):

    d = translate_substring_dict(string)
    cuts = translate_dict_cut(d, probe)
    if len(cuts) > 0:
        return cuts
    else:
        return []


def strip_var_from_criteria(string, probe):

    d = translate_substring_dict(string)

    var_list = []
    tmp_ass = {"Probe": probe, "Global": "global"}

    for key, value in d.items():
        tmp = re.split(r'\<|\<=|\>|\>=|==|\&\&|\|\|', value.replace(" ", ""))
        var_list += get_var_list(tmp[0::2], tmp_ass[key])

    return var_list


def create_bin_dictionary(binning, auto_binning, fine_binning, verbose = False):

    bin_dict = {}
    if (fine_binning):
        bin_dict = bin_dict_fine
    else:
        bin_dict = bin_dict_default


# check if custom binning is defined, overwrite respective defaults/equidistant binning
    if binning != "":
        bin_custom = binning
        bin_custom = bin_custom.split(';')

        for var_string in bin_custom:
            tmpvar = var_string.split(':')
            tmpbin = tmpvar[1].split(',')
            tmpbinlist = [float(i) for i in tmpbin]

            # fill dictionary
            bin_dict[tmpvar[0]] = tmpbinlist
        if verbose:
            OC.get_info_text("Custom binning scheme used: " + str(bin_dict))

    elif auto_binning != "":
        bin_custom = auto_binning
        bin_custom = bin_custom.split(';')

        for var_string in bin_custom:
            tmpvar = var_string.split(':')
            tmpbin = tmpvar[1].split(',')
            if (len(tmpbin) != 3):
                OC.get_warning_text(
                    "-autoBinning needs var:NoBins,min,max as an input! -autoBinning will be ignored now."
                )
            else:
                # to keep the code simple, create a dictionary of binsInFit
                bin_bins = int(tmpbin[0])
                # define bin step as  (bin_max - bin_min() / bin_bins
                bin_step = (float(tmpbin[2]) - float(tmpbin[1])) / bin_bins
                # Fill dictionary from bin_min to bin_max (that's why xrange(bin_bins + 1))
                tmpbinlist = [
                    float(float(tmpbin[0]) + i * bin_step)
                    for i in range(bin_bins + 1)
                ]
                bin_dict[tmpvar[0]] = tmpbinlist
                if verbose:
                    OC.get_info_text("Automatic binning scheme used: " +
                                     str(bin_dict))

    return bin_dict


def one_bin_bin_dict(binDict, whichBin):
    '''
        If whichBin not empty string, return a dictionary with just one bin depending on whichBin
    '''
    if (whichBin == ""):
        return binDict
    else:
        bin_dict = {}
        dict_list = whichBin.split(";")
        for one_bin in dict_list:
            var_bin_list = one_bin.split(":")
            var_list = var_bin_list[0].split("-")
            bin_list = var_bin_list[1].split("-")
            for i in range(len(var_list)):
                bin_dict[var_list[i]] = [
                    binDict[var_list[i]][int(bin_list[i])],
                    binDict[var_list[i]][int(bin_list[i]) + 1]
                ]
        OC.get_debug_text(bin_dict)
        return bin_dict


def create_var_list(variables, variables_2D, bin_dict, verbose):

    # check if binning for all variables available, else throw a warning and remove it
    new_variables = []
    for var in variables:
        if var not in bin_dict.keys():
            OC.get_warning_text(
                f'No binning given for custom variable {var}, it will be ignored!'
            )
        else:
            new_variables.append(var)
    variables = new_variables

    new_variables_2D = {}
    for key, item in variables_2D.items():
        if (item[0] in bin_dict.keys()) and (item[1] in bin_dict.keys()):
            new_variables_2D[key] = item
    variables_2D = new_variables_2D

    # split the list into global and remaining variables
    # faster due to shared info for mu+ and mu
    var_mult = []
    var_rem = []
    for i in variables:
        if i in var_list_glob:
            var_mult.append(i)
        else:
            var_rem.append(i)
    # also add missing variables from 2D set
    for key, varList in variables_2D.items():
        for var in varList:
            if var not in var_rem:
                var_rem.append(var)
    # concatenate all the variables 1D and 2D
    list_of_vars = var_rem
    list_of_vars.extend(var_mult)
    list_of_vars.extend(variables_2D.keys())

    return list_of_vars


def make_list_branches(list_of_vars: List[str]) -> List[str]:

    __branches = []
    for var in list_of_vars:
        if '-' not in var:
            if var in var_list_glob:
                __branches.append(var)
            elif var in var_list_part:
                __branches.append("muplus_" + var)
                __branches.append("muminus_" + var)
            else:
                __branches.append(var)

    return __branches

import sys
import numpy as np

import Utils
import Dictionaries
from Utilities import OutputColoring as OC
from Paths import FileUtils, Paths, adapt_path
import uproot
import re
from numba import njit
import math
from mergedeep import merge #pip3 install mergedeep



# ------------------------------------------------------
# -- Prepare data, so it can be used more efficiently later
# ------------------------------------------------------
# TODO: Add additional cuts for Run2 stripped data
#

# set charge strings 
charges = ["Plus", "Minus"]


class Weight:
    def __init__(self, __data, __mc):
        __scaled_data, __scaled_data_err = self.__normalise(__data)
        __scaled_mc, __scaled_mc_err = self.__normalise(__mc)
        
        self.ratio, self.error = self.__make_ratio(__scaled_data,__scaled_data_err, __scaled_mc, __scaled_mc_err)
        self.weights = self.__set_weights()
        self.n_bins = len(self.ratio)

    def __normalise(self, __hist):
        __new_hist = np.array(__hist, dtype=np.float64)
        __err = np.array(__hist, dtype=np.float64)
        total = sum(__hist)
        if total > 0:
            __new_hist = __hist.astype(np.float64) / total
            __err = np.abs(__hist) ** 0.5 / total
        else:
            raise ValueError(OC.get_error_text("The total number of entries in the histogram used to\
                obtain weights is <=0! Abort."))
        return __new_hist, __err

    def __make_ratio(self, __h1, __e1, __h2, __e2):
        __h = np.array(__h1, np.float64)
        __e = np.array(__h1, np.float64)

        if len(__h1) == len(__h2):
            __h = np.divide(__h1, __h2, where=__h2 != 0)

            __rel_e1 = np.divide(__e1, __h1, where=__h1 != 0.0, out=np.zeros_like(__e1))
            __rel_e2 = np.divide(__e2, __h2, where=__h2 != 0.0, out=np.zeros_like(__e2))

            cond = __rel_e1 * __rel_e2 > 0

            __e = np.where(cond, np.abs(__h) * (__rel_e1 ** 2 + __rel_e2 ** 2) ** 0.5, 0)

        else:
            ValueError(OC.get_error_text(f"The data and MC histograms cannot be divided as they have different number of bins!"))        
        return __h, __e

    def __set_weights(self): 
        __w = np.array(self.ratio, float)
        for i in range(len(self.ratio)):
            __w[i] = 1.0 if (self.ratio[i] == 0 or self.error[i]/self.ratio[i] >= 0.3 or self.ratio[i] >= 5 or self.ratio[i] < 1e-100) else self.ratio[i]
        return __w

    def __get_weight(self, index: int):
        #Index is an array of indices pointing to correct weight
        # this returns weight[index]: setting index to 1 will return the 2nd weight bin value
        # but we want to return 1
        for idx in index:
            if idx >= self.n_bins: 
                OC.get_debug_text(f"Index {idx}>= nBins {self.n_bins}")
                return 1.0
            if idx < 0:
                OC.get_debug_text(f"Index {idx}< 0")
                return 1.0

        return self.weights[index]        

    def __getitem__(self, item: int):
        return self.__get_weight(item)

#TODO: Document the way of creating these files
def data_files_for_reweighting(year, method, charge: str = "Plus"):
    return {adapt_path(Paths.getReweigtingFilePath(year)):f"Tuple{method}{charge}/DecayTree"}

def __find_bin(data, ranges, nbins):
    if data < ranges[0]:
        return -1
    elif data >= ranges[1]:
        return nbins
    else:
        return int(math.floor(nbins * (data - ranges[0]) / (ranges[1] - ranges[0])))


def __get_bins(data, ranges, nbins):
    bins = np.ones(len(data), int)
    for index, data_point in enumerate(data):
         bins[index] = __find_bin(data_point, ranges, nbins)
    return bins


@njit(parallel=True)
def __fill_hist(__hist, data, lo, hi, nbins):
    if data >=lo and data<=hi:
        __hist[__find_bin(data, lo, hi, nbins)] += 1


def __get_eta(x_1, y_1, z_1, x_2, y_2, z_2):
    x = x_2 - x_1
    y = y_2 - y_1
    z = z_2 - z_1
    mag = np.sqrt(x * x + y * y + z * z)
    cos_theta = z / mag
    return np.where(cos_theta * cos_theta < 1, -0.5 * np.log((1.0 - cos_theta) / (1.0 + cos_theta)),
                    np.where(z == 0, 0,
                             np.where(z > 0, 1e11, -1e11)))


def check_var(part: str, probe: str = 'muminus'):
    tag = 'muminus'
    if probe == tag:
        tag = 'muplus'

    if part == "Probe":
        return f"{probe}_"
    elif part == "Tag":
        return f"{tag}_"
    elif part == "Mother":
        return "J_psi_1S_"
    else:
        return ''


def make_dict_operations(op_dict, probe: str = 'muplus'):
    cut_dict = {}
    for key, value in op_dict.items():
        base = check_var(key, probe)
        value = value.replace(" ","") #I am not sure why it is needed, but it is needed
        for cut in value.split(","):
            objects = re.split(r'(==|>=|<=|!=|!|&|\)|\(|\||>|<)', cut)
            varName = f"{base}{objects[0]}"
            if (varName) in cut_dict.keys(): #Make sure when multiple cuts are required, such as mass windows, take them into account
                cut_dict[varName].append(objects[1:])
            else: cut_dict[varName] = [objects[1:]]
    return cut_dict

def make_dict_methods_operations(method_list, op_dict, probe: str = 'muplus'):
    cut_dict = {}    
    for method in method_list:
        if method in op_dict.keys():
            cut_dict[method] = make_dict_operations(op_dict[method], probe)
        else:
            raise ValueError(OC.get_error_text(f"{method} not in dictionay with selections"))

    return cut_dict

def apply_selection(chunk, report, cut_dict):
    '''
        Parse the cut_dict into True/False numpy array
    '''
    selection = np.ones(report.stop - report.start, np.bool_)
    for var, cuts in cut_dict.items():
        var = var.replace(" ","")
        for cut in cuts:
            cut[1] = cut[1].replace(" ","")            
            if cut[0] == "==":
                if cut[1].lower() in ["true", "false"]:
                    selection = np.logical_and(selection, chunk[var] == bool(cut[1]))
                else:
                    selection = np.logical_and(selection, chunk[var] == float(cut[1]))
            elif cut[0] == ">=":
                selection = np.logical_and(selection, chunk[var] >= float(cut[1]))
            elif cut[0] == "<=":
                selection = np.logical_and(selection, chunk[var] <= float(cut[1]))
            elif cut[0] == "!=" or cut[0] == "!":
                if cut[1].lower() in ["true", "false"]:
                    selection = np.logical_and(selection, chunk[var] is not bool(cut[1]))
                else:
                    selection = np.logical_and(selection, chunk[var] != float(cut[1]))
            elif cut[0] == "&":
                selection = np.logical_and(selection, chunk[var] and chunk[cut[1]])
            elif cut[0] == "|":
                selection = np.logical_and(selection, chunk[var] or chunk[cut[1]])
            elif cut[0] == ">":
                selection = np.logical_and(selection, chunk[var] > float(cut[1]))
            elif cut[0] == "<":
                selection = np.logical_and(selection, chunk[var] < float(cut[1]))
            else:
                print(f"{cut[0]} not implemented!")

    return selection

def list_of_trees(year, mode, method_list, wg_production, **kwargs):

    sim_ver = "" 
    if "sim_ver" in kwargs.keys():
        sim_ver = kwargs["sim_ver"]

    test = False
    if "test" in kwargs.keys():
        test = kwargs["test"]

    trees = {}

    files = Paths.getPrepareInput(year,mode,sim_ver,wg_production,test)


    for method in method_list:
        tree_locations = [f"Tuple{method}Plus/DecayTree",
                          f"Tuple{method}Minus/DecayTree"]

        trees[method] = {}
        trees[method][0] = {f"{adapt_path(Paths.EOS_ROOT+file)}" : f"{tree_locations[0]}" for file in files}
        trees[method][1] = {f"{adapt_path(Paths.EOS_ROOT+file)}" : f"{tree_locations[1]}" for file in files}
        
        if len(trees[method][0]) == 0 and len(trees[method][0]) == 0:
            OC.get_warning_text(f"No files found for {method} method")            
            OC.get_info_text(f"Tried {Paths.EOS_ROOT}{files}")
        # Some files in 2017 do not have T method, which causes a crash
        # So here check for the existence of the tuple 
        if (year in Dictionaries.year_missing_T_method and method == "T"): 
            for path in list(trees[method][0].keys()):
                file = uproot.open(adapt_path(path))
                if "TupleTPlus;1" not in file.keys(): del trees[method][0][path]
                if "TupleTMinus;1" not in file.keys(): del trees[method][1][path]
                                
    return trees

def Prepare(mode, opts):
    # weight ranges
    year = opts.year
    method_list = opts.method
    verbose = opts.verbose
    max_entries = opts.max_entries


    # start printout
    Utils.multiple_lines_title('Welcome to TrackEffReducer', 'Running on %s' % mode.ljust(4))

    # default variables lists to keep on tuple (split in global/particle vars)
    var_list = Dictionaries.get_default_var_list()

    trees = list_of_trees(year=year,
                          mode=mode,
                          method_list=opts.method,
                          wg_production=opts.WGProduct,
                          sim_ver=opts.sim_ver,
                          verbose=verbose,
                          test=opts.test)

    # dictionary of branches in the original TTree
    # this is automatically generated from the list of variables `var_list` and
    # inspecting one TTree in the list of `trees`

    tmp_key = list(trees[method_list[0]][0])[0]
    tmp_filename = adapt_path(f"{tmp_key}:{trees[method_list[0]][0][tmp_key]}")
    tmp_file = uproot.open(tmp_filename)
    branch_orig_dict = {}

    for branch in var_list:
        if branch in tmp_file.keys():
            branch_orig_dict[branch] = str(tmp_file[branch].interpretation.from_dtype)
        else:
            OC.get_info_text(f"Branch {branch} not present in the original production files!")

    # new TTree contains additional branches
    branch_dict = branch_orig_dict.copy()
    branch_dict.update({"Mother_ETA": '>f8', "weight": '>f8', "matched": '|b1'})

    ### START OF CUSTOM CUT READ-IN ###

    if opts.cuts != "":
        # cut dictionary for custom cuts, split by type
        opts.cuts = Utils.insert_Global_cut(opts.cuts,Dictionaries.get_nSPDHits_cut(year,opts.high_multi))
        custom_cut_dict = Dictionaries.translate_string_dict(opts.cuts)
        cuts = [make_dict_methods_operations(method_list, custom_cut_dict, 'muplus'),
                make_dict_methods_operations(method_list, custom_cut_dict, 'muminus')]
    else:
        # default cut dictionary 
        cuts = [make_dict_methods_operations(method_list, Dictionaries.get_default_cuts_dict(year,opts.high_multi), 'muplus'),
                make_dict_methods_operations(method_list, Dictionaries.get_default_cuts_dict(year,opts.high_multi), 'muminus')]

    if opts.verbose:
        print()
        OC.get_info_text(f"Positive track cuts string used:")
        print(*[str(method) + ': ' + str(cut) for method,cut in cuts[0].items()], sep = "\n\n")
        print()
        OC.get_info_text(f"Negative track cuts string used:")
        print(*[str(method) + ': ' + str(cut) for method,cut in cuts[1].items()], sep = "\n\n")
        print()

    ### START OF MATCHING CRITERIA READ-IN ###
    # cut dictionary for custom cuts, split by type
    # add extra custom matching requirements from option, either on probe track or globally
    # example '-matchCrit "Probe:Matched_GhostProb<0.3;Global:nSPDHits<400" '
    # prepare the list of default matching criteria
    matching_criteria = [make_dict_methods_operations(method_list, Dictionaries.default_matching_criteria, 'muplus'),
                         make_dict_methods_operations(method_list, Dictionaries.default_matching_criteria, 'muminus')]

    if opts.match_crit != "":
        extra_matching_dict = Dictionaries.translate_string_dict(opts.match_crit)
         
        extra_matching = [make_dict_methods_operations(method_list, extra_matching_dict, 'muplus'),
                          make_dict_methods_operations(method_list, extra_matching_dict, 'muminus')]

        # extend default matching criteria
        for charge in [0,1]:            
            merge(matching_criteria[charge],matching_criteria[charge],extra_matching[charge])
                  
    if verbose:
        OC.get_info_text(f"Positive particle matching criteria used:")
        print(*[str(method) + ': ' + str(cut) for method,cut in matching_criteria[0].items()], sep = "\n\n")
        print()
        OC.get_info_text(f"Negative particle matching criteria used:")
        print(*[str(method) + ': ' + str(cut) for method,cut in matching_criteria[1].items()], sep = "\n\n")
        print()


    ### END OF MATCHING CRITERIA READ-IN ###

    # list of min, max range for some standard weighting vars #Increased, just in case
    nSPDHits_edges = Dictionaries.get_nSPDHits_list(year)
    weight_range = {
        "nSPDHits": (nSPDHits_edges[opts.high_multi], nSPDHits_edges[opts.high_multi+1]),
        "nTracks": (0, 600),
        "nLongTracks": (0, 2000),
        "nVeloTracks": (0, 400),
        "nTTracks": (0, 6000)
    }

    OC.get_debug_text(f"weight_range[nSPDHits] {weight_range['nSPDHits']}")
    weight_nBins = 200

    for method in method_list:

        cache = uproot.LRUArrayCache("200 MB")

        # produce immediately the multiplicity histogram
        __weights = [None, None]
        if mode == "MC":
            mc_arrays = [None, None]
            data_arrays = [None, None]
            __mc_hist = [None, None]
            __data_hist = [None, None]

            for charge in (0, 1):
                if (opts.test):
                    small_tree_dict = {k:  trees[method][charge][k] for k in list( trees[method][charge])[:10]}
                    mc_arrays[charge] = uproot.dask(small_tree_dict, filter_name=opts.weight_var, step_size=10*max_entries, cache=cache, max_num_elements = max_entries )
                else:
                    OC.get_info_text("Loading the full MC sample to get the weights. .  .  .")
                    mc_arrays[charge] = uproot.dask(trees[method][charge], filter_name=opts.weight_var, step_size=10000, cache=cache, max_num_elements = 1000000)
                    OC.get_info_text("Done!")

                OC.get_debug_text(f"Got {len(mc_arrays[charge].compute()[opts.weight_var])} MC events for weighting.")

                __mc_hist[charge], __bins = np.histogram(mc_arrays[charge].compute()[opts.weight_var].to_numpy(),
                                                         bins=weight_nBins,
                                                         range=weight_range[opts.weight_var])

                data_arrays[charge] = uproot.dask(data_files_for_reweighting(opts.year, method, charges[charge]),
                                                  filter_name=opts.weight_var, step_size=10000, cache=cache)

                OC.get_debug_text(f"Got {len(data_arrays[charge].compute()[opts.weight_var])} Data events for weighting.")

                __data_hist[charge], __bins = np.histogram(data_arrays[charge].compute()[opts.weight_var].to_numpy(),
                                                           bins=weight_nBins,
                                                           range=weight_range[opts.weight_var])                                                                      
                __weights[charge] = Weight(__data_hist[charge], __mc_hist[charge])

        
        #info_message("Entries in the sample: " + str(actual_entries))
        if opts.max_entries != -1 or verbose:
            OC.get_info_text(f"Maximum number of entries requested: {opts.max_entries}")

        #create files
        outputFilePath = Paths.getPrepareOutput(year=year,
                                                mode=mode,
                                                method=method,
                                                official=opts.store_official,
                                                wg_production=opts.WGProduct,
                                                sim_ver=opts.sim_ver,
                                                match_crit=opts.match_crit,
                                                high_multi=opts.high_multi,
                                                verbose=verbose)
        if not FileUtils.CheckIfCanOverwrite(outputFilePath,opts.force):
            sys.exit(1)                       

        file = uproot.recreate(outputFilePath)
        for charge in (0,1):
            file.mktree(f"TrackEffTree{charges[charge]}{method}",branch_dict)
            for chunk, report in uproot.iterate(trees[method][charge], branch_orig_dict,
                                                report=True, cache=cache):

                selection = apply_selection(chunk, report, cuts[charge][method])

                chunk['Mother_ETA'] = __get_eta(chunk['J_psi_1S_OWNPV_X'],
                                                chunk['J_psi_1S_OWNPV_Y'],
                                                chunk['J_psi_1S_OWNPV_Z'],
                                                chunk['J_psi_1S_ENDVERTEX_X'],
                                                chunk['J_psi_1S_ENDVERTEX_Y'],
                                                chunk['J_psi_1S_ENDVERTEX_Z'])   

                chunk['matched'] = apply_selection(chunk, report, matching_criteria[charge][method])
                if (mode == "Data"):
                    chunk['weight'] = np.ones(len(chunk['matched']), int) #Data have weight 1
                else:
                    chunk['weight'] = __weights[charge][__get_bins(chunk[opts.weight_var],
                                                               weight_range[opts.weight_var],
                                                               weight_nBins)]
                file[f"TrackEffTree{charges[charge]}{method}"].extend(chunk[branch_dict.keys(), selection])

                if opts.max_entries >0 and report.stop >= opts.max_entries:
                    break

            if opts.max_entries > 0:
                OC.get_info_text(f"{opts.max_entries} requested, {report.stop} processed for method {method} and charge {charges[charge]}!")

            OC.get_info_text(f"{report.stop} events processed for method {method} and charge {charges[charge]}!")


    OC.get_ok_text("Preparation finished.")

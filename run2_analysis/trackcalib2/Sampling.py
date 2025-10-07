from typing import Dict, List, Tuple, Union

import numpy as np

from Dictionaries import var_list_glob
import pandas
import pandas as pd
import re


# Build varname for not global variables
def varname(var: str, glob: bool, charge: int):
    if glob:
        return var
    else:
        if charge == -1:
            return "muminus_%s" % var
        elif charge == 1:
            return "muplus_%s" % var


# Helper for pandas' DataFrame - bin definition
def sel(df, var: str, low: float, high: float):
    return np.logical_and(df[var] > low, df[var] <= high)

def add_sel(selection_1, selection_2):
    return np.logical_and(selection_1, selection_2)

# Helper for pandas' DataFrame - matching function
def matching(df, matching_var: str = "matched"):
    return (df[matching_var] == 1)


def subsample(df: dict, bin_dict: dict, variables: List[str], obs: str, index: Union[int, List[int]], charge: int = 0,
              matching_var: str = "matched") -> Tuple[List, List]:

    datasets = [[], []]

    if charge == 0:
        __charges = [-1, 1]
    else:
        __charges = [charge]

    for __charge in __charges:
        m = matching(df[__charge], matching_var)

        b = True
        for i, variable in enumerate(variables):
            glob = True if (variable in var_list_glob) else False
            var_name = varname(variable, glob, __charge)

            lo = bin_dict[variable][index[i]]
            hi = bin_dict[variable][index[i] + 1]
            b = add_sel(b, sel(df[__charge], var_name, lo, hi))

        datasets[0].extend(df[__charge][m & b][obs].to_numpy())
        datasets[1].extend(df[__charge][~m & b][obs].to_numpy())

    return datasets[0], datasets[1]


# core of the sampler
def sampling(df, ds_list: dict, var: str, obs: str, bin_dict: dict, charge: int, matching_var: str = "matched"):
    # define matching criteria
    m = matching(df[charge], matching_var)

    # do 1d splitting
    if "-" not in var:

        # if variable is a global variable (e.g. nTracks) the suffix is muplus or muminus is not added
        glob = True if (var in var_list_glob) else False

        # variable name with the correct suffix
        varn = varname(var, glob, charge)

        # dictionary with datasamples
        ds_list[var][charge] = {}
        ds_list[var][charge]["matched"] = []
        ds_list[var][charge]["failed"] = []
        ds_list[var][charge]["total"] = []
        for bin in range(len(bin_dict[var]) - 1):
            lo = bin_dict[var][bin]
            hi = bin_dict[var][bin + 1]
            b = sel(df[charge], varn, lo, hi)

            ds_list[var][charge]["matched"].append((df[charge][m & b][obs]).to_numpy())
            ds_list[var][charge]["failed"].append((df[charge][~m & b][obs]).to_numpy())
            ds_list[var][charge]["total"].append((df[charge][b][obs]).to_numpy())

    # do 2d splitting
    else:
        variables = var.split("-")
        iglob = True if (variables[0] in var_list_glob) else False
        jglob = True if (variables[1] in var_list_glob) else False
        ivarn = varname(variables[0], iglob, charge)
        jvarn = varname(variables[1], jglob, charge)
        ds_list[var][charge] = {}
        ds_list[var][charge]["matched"] = []
        ds_list[var][charge]["failed"] = []
        ds_list[var][charge]["total"] = []

        for ibin in range(len(bin_dict[variables[0]]) - 1):

            ilo = bin_dict[variables[0]][ibin]
            ihi = bin_dict[variables[0]][ibin + 1]
            ib = sel(df[charge], ivarn, ilo, ihi)

            for jbin in range(len(bin_dict[variables[1]]) - 1):
                jlo = bin_dict[variables[1]][jbin]
                jhi = bin_dict[variables[1]][jbin + 1]
                jb = sel(df[charge], jvarn, jlo, jhi)

                ds_list[var][charge]["matched"].append((df[charge][m & ib & jb][obs]).to_numpy())
                ds_list[var][charge]["failed"].append((df[charge][~m & ib & jb][obs]).to_numpy())
                ds_list[var][charge]["total"].append((df[charge][ib & jb][obs]).to_numpy())


# Produce 2 (positive and negative charges) dictionaries with subsamples for each bin
def make_datasets(df: Dict[int, pd.DataFrame], list_vars: List[str], bin_dict: dict, obs: str,
                  matching_var: str = "matched"):
    ds_dict = {}
    for var in list_vars:
        ds_dict[var] = {}
        for charge in [-1, 1]:
            sampling(df, ds_dict, var, obs, bin_dict, charge, matching_var)

    return ds_dict


# Helper to apply selection in ROOT fashioned way to a pandas DataFrame
def select(df: pandas.DataFrame, selection):
    selection_list = re.split('([==|&|!=|>|<|>=|<=|\||\(|\)|\ ])', selection)
    new = ""
    for term in selection_list:
        if term in (df.columns.tolist()):
            new += "df[\"%s\"]" % term
        else:
            new += term
    return eval(new)


def array_select(_array, selection):
    """selection for awkward arrays"""
    selection_list = re.split('([==|&|!=|>|<|>=|<=|\||\(|\)|\ ])', selection)
    new = ""
    for term in selection_list:
        if term in (_array.fields):
            new += "_array[\"%s\"]" % term
        else:
            new += term
    return eval(new)

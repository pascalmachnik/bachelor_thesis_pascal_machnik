import yaml


def df_define(df, define):
    """
    Apply a definition to a DataFrame based on a configuration file.

    Parameters:
    df (RDataFrame): The DataFrame to define.
    define (str): The name of the definition to apply, which corresponds to a key in the configuration file.

    Returns:
    RDataFrame: The defined DataFrame.
    """
    # read the configuration file
    with open("analysis_config.yml", "r") as config_file:
        config = yaml.safe_load(config_file)

    # get the definition from the config
    filter = config[define]

    # apply the definition to the DataFrame
    df = df.Define(define, filter)

    return df


def BDT_define(df, corrected=False):
    """
    Define BDT variables in the DataFrame.

    Parameters:
    df (RDataFrame): The DataFrame to define BDT variables in.

    Returns:
    RDataFrame: The DataFrame with BDT variables defined.
    """
    # read the configuration file
    with open("analysis_config.yml", "r") as bdt_config_file:
        bdt_config = yaml.safe_load(bdt_config_file)

    definitions = {
        "q2": "(pow(Lambdab_PE-Lambda_PE,2) - pow(Lambdab_PX-Lambda_PX,2) - pow(Lambdab_PY-Lambda_PY,2) - pow(Lambdab_PZ-Lambda_PZ,2))/1000000",
        "log10_Lambdab_pT": "log10(abs(Lambdab_PT))",
        "log10_Proton_pT": "log10(abs(Proton_PT))",
        "log10_Pion_pT": "log10(abs(Pion_PT))",
        "log10_L1_pT": "log10(abs(L1_PT))",
        "ProbNNmu_max": "max(L1_ProbNNmu, L2_ProbNNmu)",
        "ProbNNmu_min": "min(L1_ProbNNmu, L2_ProbNNmu)",
        "chi2ip_max": "max(L1_IPCHI2_OWNPV, L2_IPCHI2_OWNPV)",
        "chi2ip_min": "min(L1_IPCHI2_OWNPV, L2_IPCHI2_OWNPV)",
        "chi2ip_min_pi_p": "min(Pion_IPCHI2_OWNPV, Proton_IPCHI2_OWNPV)",
        "log10_chi2ip_max": "log10(chi2ip_max)",
        "log10_chi2ip_min": "log10(chi2ip_min)",
        "log10_chi2ip_min_pi_p": "log10(chi2ip_min_pi_p)",
        "log10_Lambda_pT": "log10(abs(Lambda_PT))",
        "log10_Lambdab_DIRA_OWNPV": "log10(Lambdab_DIRA_OWNPV)",
        "log10_ProbNNmu_max": "log10(ProbNNmu_max)",
        "log10_ProbNNmu_min": "log10(ProbNNmu_min)",
    }

    if corrected:
        definitions["ProbNNmu_max"] = "max(L1_ProbNNmu_pidcorr_default, L2_ProbNNmu_pidcorr_default)"
        definitions["ProbNNmu_min"] = "min(L1_ProbNNmu_pidcorr_default, L2_ProbNNmu_pidcorr_default)"

    # apply the definitions to the DataFrame
    for name, expression in definitions.items():
        df = df.Define(name, expression)

    return df

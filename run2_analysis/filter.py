import yaml


def df_filter(df, filter):
    """
    Apply a filter to a DataFrame based on a configuration file.

    Parameters:
    df (RDataFrame): The DataFrame to filter.
    filter (str): The name of the filter to apply, which corresponds to a key in the configuration file.

    Returns:
    RDataFrame: The filtered DataFrame.
    """
    # read the configuration file
    with open("analysis_config.yml", "r") as config_file:
        config = yaml.safe_load(config_file)

    # get the filter from the config
    filter = config[filter]

    # apply the filter to the DataFrame
    df = df.Filter(filter)

    return df

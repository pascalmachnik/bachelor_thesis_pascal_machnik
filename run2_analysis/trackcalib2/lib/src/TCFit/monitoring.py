from abc import abstractmethod
from typing import List, Tuple, Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ..TCFit.parameter import Parameter
import seaborn as sns

##########################################################################################################
#                                   Interface Monitor Implementation                                     #
##########################################################################################################


class IFMonitor:
    """For debugging it might be useful to monitor the development of 
    the parameters when being updated. For this purpose, this monitor is
    build
    
    Methods (abstract)
    ------------------
        display_parameter_evolution(epochs : int)
            The recorded parameters are being plotted on a rectangular grid
            and displayed to the user
    """


    @abstractmethod
    def display_parameter_evolution(self, epochs: int) -> None:
        """When the parameter update finished, this method will plot the recorded
        parameters over the epochs it took to converge to a solution

        Parameters
        ----------
            epochs: int
                `epochs` is the number of iterations which were needed
                to optimise the parameters
        """
        pass


##########################################################################################################
#                                           Monitor Implementation                                       #
##########################################################################################################


class MonitorLikeLihood(IFMonitor):
    """For debugging it might be useful to monitor the development of 
    the parameters when being updated. For this purpose, this monitor is
    build

    Parameters
    ----------
        parameters : list[utils.parameter.Parameter]
            The list of parameters which should be recorded and plotted afterwards
        monitor_on : bool
            If set to `False`, the plots will not be created, if set to `True`, the plots
            corresponding to the parameters are displayed to the user
        theme      : str, default = darkgrid
            The seaborn theme in order to influence the look of the plots

    Attributes
    ----------
        loss_metrics     : dict[str, floats]
            The following loss metrics are stored: the log-likelihood,
            the sum of normalisations and the loss metric
        monitor_elements : int
            The number of parameters which have to be monitored

    Methods:
    --------
        display_parameter_evolution(epochs : int)
            The recorded parameters are being plotted on a rectangular grid
            and displayed to the user
    """


    def __init__(self, parameters: List[Parameter], monitor_on: bool, theme: str = "darkgrid"):
        self.parameters            = {parameter.name: [parameter.value] for parameter in parameters}
        self.loss_metrics          = {"NLL": [np.nan], "Norm": [np.nan], "Loss": [np.nan]}
        self.monitor_elements      = len(self.parameters) + 3
        self.monitor_on            = monitor_on

        sns.set_theme(style = theme)


    def display_parameter_evolution(self, epochs: int) -> None:
        """When the parameter update finished, this method will plot the recorded
        parameters over the epochs it took to converge to a solution

        Parameters
        ----------
            epochs: int
                `epochs` is the number of iterations which were needed
                to optimise the parameters
        """

        if self.monitor_on:
            dataset, column_names = self.__prepare_dataset(epochs)
            grid_row_dim = int(self.monitor_elements / 2)
            if self.monitor_elements % 2 != 0:
                grid_row_dim = int((self.monitor_elements + 1) / 2)
            grid_col_dim = 2
            plt.rc("axes", labelsize = 10)
            _, axes = plt.subplots(nrows = grid_row_dim, ncols = grid_col_dim, figsize = (12, 6), sharex = "all")
            count   = 0
            for row in range(grid_row_dim):
                for col in range(grid_col_dim):
                    if count >= self.monitor_elements:
                        break
                    sns.lineplot(x = "epochs", y = f"{column_names[count]}", data = dataset, ax = axes[row, col])
                    count += 1
            plt.suptitle("Evolution of updated parameters and loss metrics:", fontsize = 20)
            plt.show()


    def __prepare_dataset(self, epochs: int) -> Tuple[pd.DataFrame, Sequence[str]]:
        """In order to plot the data, a dataset has to be created, the dataset
        consists of the epoch sequence and the parameter values with their respective
        labels

        Parameters
        ----------
            epochs: int
                `epochs` is the number of iterations which were needed
                to optimise the parameters

        Returns
        -------
            pd.DataFrame, array-like[str]
                A pandas dataframe is returned which stores all the information as
                the parameter values, the epochs, ... . The array contains the labels
                of the columns
        """
        parameter_names  = [name for name in self.parameters.keys()]
        parameter_values = [value for value in self.parameters.values()]
        self.loss_metrics["NLL"].append(np.nan)
        self.loss_metrics["Norm"].append(np.nan)
        self.loss_metrics["Loss"].append(np.nan)
        epochs = range(epochs)
        data   = []
        for index in range(len(epochs)):
            parameter_sequence = []
            for subindex in range(len(parameter_values)):
                parameter_sequence.append(parameter_values[subindex][index])
            data.append([epochs[index], *parameter_sequence, self.loss_metrics["NLL"][index],
                         self.loss_metrics["Norm"][index], self.loss_metrics["Loss"][index]])
        dataset = pd.DataFrame(data = data, columns = ["epochs", *parameter_names, "NLL", "Sum of Norm.", "Loss"])

        return dataset, [*parameter_names, "NLL", "Sum of Norm.", "Loss"]


    def __call__(self, name, value):
        if name in ["NLL", "Norm", "Loss"]:
            self.loss_metrics[name].append(value)
            return
        self.parameters[name].append(value)
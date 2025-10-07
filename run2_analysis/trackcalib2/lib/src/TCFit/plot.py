import pickle
from abc import abstractmethod
from statistics import mode
from typing import Union

import matplotlib.patches as mpl_patches
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib import ticker
import mplhep
import numpy as np
from ..TCFit.model import AddPDF, Model
from ..TCFit.stats import eval_unc
from ..TCFit.latex_format import latex_obs_names
from uncertainties import ufloat


v_eval_unc = np.vectorize(eval_unc)

class IF_LikelihoodPlot:
    @abstractmethod
    def plot(self):
        pass


class Plot(IF_LikelihoodPlot):
    __width = 14
    __height = 8
    __linewidth = 3
    __nplots = 1
    __maxcols = 4
    __fig = None
    __leg_size = 24

    def __init__(self, minuit: object, model: Model, data: Union[np.ndarray, list],
                 limit: Union[list, tuple], n_bins: int):
        self.limit = limit
        self.n_bins = n_bins
        self.data = data
        self.model = model
        self.minuit = minuit
        self.bin_width = 1.

    def __setup(self):
        plt.style.use(["default", mplhep.style.LHCb2])

        self.x_min, self.x_max = self.limit[0], self.limit[1]

        self.bin_width = (self.x_max - self.x_min) / self.n_bins

        self.value_y = {}
        self.value_y_error = {}
        self.value_x = {}
        edges = {}
        for i in range(len(self.data)):
            self.value_y[i], edges[i] = np.histogram(self.data[i], bins=self.n_bins, range=self.limit)
            self.value_x[i] = 0.5 * (edges[i][1:] + edges[i][:-1])
            self.value_y_error[i] = v_eval_unc(self.value_y[i])


    def addText(axis, text, pos = [0.25,0.25], size = __leg_size):
        # transform is needed in order to be able to add the text relative to the box coordinates
        # instead of using the x- and y- axis values as coordinates
        axis.text(pos[0],pos[1], text, style = 'normal', size = size, transform=axis.transAxes)
        return
    
    def format_latex(self, name):
        '''
        Turn the name of the parameter into a latex string that looks better on the plots.
        If the name is not in the dictionary, return the name as is.
        '''
        if name in latex_obs_names:
            return latex_obs_names[name]
        else:
            return name
        
    def round_value(self, value, error, fixed = False):
        if fixed:
            return str(ufloat(value,error).n) + " (fixed)"
        else:
            return str(ufloat(value,error)).replace( "+/-"," $\\pm$ ")

    def calculate_pull(self, model_y, data_y, data_y_unc):
        pulls = np.zeros(len(model_y))
        pull_uncs = np.zeros((2, len(model_y)))
        for i in range(len(model_y)):
            # Case 1: Model value is smaller than data value
            # --> Relevant uncertainty is the lower one
            if model_y[i] < data_y[i]:
                pulls[i] = (data_y[i] - model_y[i]) / data_y_unc[0][i]
                pull_uncs[0][i] = 1
                pull_uncs[1][i] = data_y_unc[1][i]/data_y_unc[0][i]
            # Case 2: Model value is larger than data value
            # --> Relevant uncertainty is the upper one
            else:
                pulls[i] = (data_y[i] - model_y[i]) / data_y_unc[1][i]
                pull_uncs[0][i] = data_y_unc[0][i]/data_y_unc[1][i]
                pull_uncs[1][i] = 1
                
        return pulls, pull_uncs

    def plot(self, show_components: bool = True, show_fit_result: bool = True, muon_charge: int = 0,
            tag = None, mode = None, method = None, simFit = False, binnedFit = False):
        """ The given model and data are plotted """
        self.__setup()

        x = np.linspace(self.x_min, self.x_max, self.n_bins * 10)
        y = np.zeros((len(self.model.pdfs), x.shape[0])) 
        y_pull = np.zeros((len(self.model.pdfs), self.value_x[0].shape[0]))
        pull = np.zeros((len(self.model.pdfs), self.value_x[0].shape[0]))
        pull_unc = np.zeros((len(self.model.pdfs), 2, self.value_x[0].shape[0]))

        n_rows = (int(len(self.model.pdfs) / self.__maxcols) + 1)
        n_cols = (self.__maxcols if len(self.model.pdfs) > self.__maxcols else len(self.model.pdfs))
        # build a plot with multiple subplots
        self.__fig, axs = plt.subplots(nrows=n_rows*2,
                                       ncols=n_cols,
                                       figsize=(self.__width * n_cols, self.__height * n_rows*1.33),
                                       gridspec_kw={'height_ratios': [1, 0.33]*n_rows})
        self.__fig.tight_layout()
        # add more space between the subplots
        self.__fig.subplots_adjust(wspace=0.15)
        self.__fig.subplots_adjust(hspace=0)
        
        # Ensure axs is always a 2D array for consistent indexing
        if n_cols == 1:
            axs = np.expand_dims(axs, axis=1)
        
        # Remove ticks from the upper subplots
        for ax in axs[0, :]:
            ax.tick_params(axis='both', which='both', labelbottom=False)  
      
        for i in range(n_cols):
            axs[0][i].sharex(axs[1, i])

        if not isinstance(axs, np.ndarray):
            axs = np.array([axs])

        # set axs ranges
        for ax in axs[0][:]:
            ax.set_xlim(self.limit)
            #ax.set_ylim(0, 1.1 * max([max(self.value_y[i]) for i in range(len(self.data))]))
        for ax in axs[1][:]:
            ax.set_xlim(self.limit)
            
        # set help lines for the pull plots
        for ax in axs[1][:]:
            ax.axhline(-5, color='black', lw=0.5, ls='-.')
            ax.axhline(-3, color='black', lw=0.5, ls='--')
            ax.axhline( 0, color='black', lw=0.5)
            ax.axhline( 3, color='black', lw=0.5, ls='--')
            ax.axhline( 5, color='black', lw=0.5, ls='-.')
        # for each category
        for index in range(len(self.model.pdfs)):
            y[index] = self.model.pdfs[index].evaluate(x) * self.bin_width
            y_pull[index] = self.model.pdfs[index].evaluate(self.value_x[index]) * self.bin_width
            pull[index], pull_unc[index] = self.calculate_pull(
                model_y=y_pull[index], 
                data_y=self.value_y[index], 
                data_y_unc=self.value_y_error[index])

            sub_y = {}
            
            if hasattr(self.model.pdfs[index],'pdfs'):
                for shape in self.model.pdfs[index].pdfs:
                    sub_y[shape.name] = shape.evaluate(x) * self.bin_width
                    if isinstance(shape, AddPDF):
                        for sub_shape in shape.pdfs:
                            sub_y[sub_shape.name] = sub_shape.evaluate(x) * self.bin_width

            # build the text string which contains all the fit parameters
            fit_info = []
            fit_info_main = []
            
            # in simFit want on the pass side the efficiency and the number of passed signal and background events first
            # on the fail side the efficiency and the number of failed signal and background events first
            if simFit:
                if "f_pass" in self.model.pdfs[index].get_coefficient_names:
                    rounded_value = self.round_value(self.minuit.values["efficiency_sig"], 
                                                     self.minuit.errors["efficiency_sig"])
                    fit_info_main.append(f"{self.format_latex("efficiency_sig")} = {rounded_value}\n")
                    
                    for param in ["f_pass", "f_pass_bkg"]:
                        rounded_value = self.round_value(self.model.pdfs[index].get_coefficient(param).value,
                                                         self.model.pdfs[index].get_coefficient(param).error)
                        fit_info_main.append(f"{self.format_latex(param)} = {rounded_value}\n")
                        
                elif "f_fail" in self.model.pdfs[index].get_coefficient_names:
                    rounded_value = self.round_value(self.minuit.values["efficiency_bkg"], 
                                                     self.minuit.errors["efficiency_bkg"])
                    fit_info_main.append(f"{self.format_latex("efficiency_bkg")} = {rounded_value}\n")
                    
                    for param in ["f_fail", "f_fail_bkg"]:
                        rounded_value = self.round_value(self.model.pdfs[index].get_coefficient(param).value,
                                                         self.model.pdfs[index].get_coefficient(param).error)
                        fit_info_main.append(f"{self.format_latex(param)} = {rounded_value}\n")
                else:
                    raise ValueError("Simultaneous fit parameters not found")
            else:
                for param in ["signal_yield", "background_yield"]:
                    rounded_value = self.round_value(self.minuit.values[param], self.minuit.errors[param])
                    fit_info_main.append(f"{self.format_latex(param)} = {rounded_value}\n")
            
            for name, value, error in zip([param.name  for param in self.model.pdfs[index].parameters], 
                                          [param.value for param in self.model.pdfs[index].parameters], 
                                          [param.error for param in self.model.pdfs[index].parameters]):
                if name in ["efficiency_sig", "signal_yield", "efficiency_bkg", "background_yield"]:
                    continue
                if simFit and index == 1 and name[-2:] != "_F":
                    continue
                #TODO: maybe at some point fix integers to be integers, but whatever at this point
                rounded_value = self.round_value(value, error, self.minuit.fixed[name])
                #rounded_value = "${:L}$".format(ufloat(value,error))
                fit_info.append(f"{self.format_latex(name)} = {rounded_value}\n")

            # plot full pdf
            axs[0][index].plot(x, y[index], color='b', lw=self.__linewidth, ls='-', marker=None)
            
            # plot sub-components of the pdf
            if show_components:
                for shape_name, shape_values in sub_y.items():
                    axs[0][index].plot(x, shape_values, lw=self.__linewidth, ls='--')

            # plot data
            axs[0][index].errorbar(self.value_x[index], self.value_y[index], 
                                yerr=self.value_y_error[index], xerr = self.bin_width/2.0, 
                                fmt=".", color='k')
            
            # plot pull
            axs[1][index].errorbar(self.value_x[index], pull[index], 
                                yerr=pull_unc[index], xerr = self.bin_width/2.0, 
                                fmt=".", color='k')
            axs[1][index].set_ylabel = r"Pull [$\sigma$]"

            #Make the y-axis tick labels reasonably looking
            # axs[0][index].yaxis.set_major_formatter(mticker.FormatStrFormatter('%.1e'))
            axs[0][index].ticklabel_format(axis="y", style="sci", scilimits=(0,0), useMathText=True)
            #axs[0][index].tick_params(axis='y', which='major', labelsize=18)

            #Add the tags
            y_pos = 0.9
            if (tag != None):    
                Plot.addText(axs[0][index], tag[index], [0.05,y_pos], 33)
                y_pos = y_pos - 0.07
            if (mode != None):   
                Plot.addText(axs[0][index], mode, [0.05,y_pos], 30)
                y_pos = y_pos - 0.07
            if (method != None): 
                Plot.addText(axs[0][index], method + " method", [0.05,y_pos], 30)
                y_pos = y_pos - 0.07
            if (muon_charge == 1):
                Plot.addText(axs[0][index], "Muon charge plus", [0.05,y_pos], 30)
                y_pos = y_pos - 0.07
            elif (muon_charge == -1):
                Plot.addText(axs[0][index], "Muon charge minus", [0.05,y_pos], 30)
                y_pos = y_pos - 0.07
            if (binnedFit):      
                Plot.addText(axs[0][index],"BinnedFit", [0.05,y_pos], 30)
                y_pos = y_pos - 0.04

            #Add axis labels
            x_title = self.model.observable.latex
            if self.model.observable.unit is not None:
                x_title = f"{self.model.observable.latex} [{self.model.observable.unit}]"
            axs[1][index].set_xlabel(x_title)
            y_label = f"Entries"
            if self.model.observable.unit is not None:
                _bin_width = round(self.bin_width) if self.bin_width.is_integer else round(self._bin_width,1)
                y_label = f"Entries / ({_bin_width} {self.model.observable.unit})"
            axs[0][index].set_ylabel(y_label)

            handles = [mpl_patches.Rectangle((0.0, 0), 1, 1, fc="white", ec="white",
                                         lw=0, alpha=0)] * len(fit_info_main+fit_info) 

            if show_fit_result:
                axs[0][index].legend(handles, fit_info_main+fit_info,
                            loc = "upper right", 
                            bbox_to_anchor = [0.8,0.5,0.2,0.5], #x0,y0, width, height
                            fontsize=self.__leg_size,
                            labelspacing = -0.25
                            )
        
        return self.__fig

    def save(self, name: str = "default"):
        self.__fig.savefig(f"./{name}")

    def delete(self):
        plt.close(self.__fig)

    def flush(self, name: str = "default"):
        self.save(name)
        self.delete()

    def __del__(self):
        plt.cla()
        plt.close(self.__fig)

    def dump(self, file):
        pickle.dump(self.__fig, file)



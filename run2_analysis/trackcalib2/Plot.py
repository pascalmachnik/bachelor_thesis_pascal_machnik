import matplotlib.pyplot as plt
from PlotDesign import colorlist, PlotDesign

from Utilities import OutputColoring as OC
from Utils import readPickleFile, binEdges

import Dictionaries
from PlotMaker import singlePlot, singlePlotUtils

from Paths import Paths, FileUtils
import numpy as np
import pandas as pd
import os
from PlotHelpers import PlotHelpers

######################################################################################
#                                    Plotter                                         #
######################################################################################

width = 10
height = 6  #TODO move to design


class PlotReader:

    @staticmethod
    def getInfoAndEfficiency(filename, setBinDict_list):
        '''
            Returns the dictionary with the fit information and 
            the dictionary with the efficiency plot information
            from the pickle file
        '''
        OC.get_info_text(f"Opening file {filename}")

        #Check if there are any values to be set by hand:
        var_list = []
        binValue_list = []
        change_dict = {}
        for dict in setBinDict_list:
            if os.path.abspath(dict["file"]) == os.path.abspath(
                    filename):  #in to prevent issues with ./ and stuff
                var_list.append(dict["var"])
                binValue_list.append([dict["bin"], dict["value"]])
        #Init the dicitonary
        for var in var_list:
            change_dict[var] = []

        for i, var in enumerate(var_list):
            change_dict[var].append(binValue_list[i])

        #Read from the file
        items = list(readPickleFile(filename))
        if (items[0]['mode'] == "Data"):
            items[0]["sim_ver"] = ""
        return items[0], items[1], change_dict

    @staticmethod
    def loadEfficiency(pkl_dict, set_values_dict):
        # The data input looks like this
        #{'P': {0: {'boundaries': [5000.0, 10000.0],
        #            'efficiency': [-1, -1, -1]
        #            },
        #        1: {'boundaries': [10000.0, 20000.0],
        #            'efficiency': [0.9768466856190094, 0.008403095254462034, -0.008403095254462034]
        #            },
        #        ...
        #    }
        #}
        #P-ETA {0: {'boundaries': ([5000.0, 10000.0], [1.9, 3.2]), 'efficiency': [0.9739229024943311, 0.00545530847345922, 0.006624957663235576]},
        #       1: {'boundaries': ([5000.0, 10000.0], [3.2, 4.9]), 'efficiency': [0.9090909090909091, 0.09090909090909094, 0.19610029414485486]},
        #       ...
        #       }

        plotDict = {}
        for var, bins in pkl_dict.items():
            bin_id = []
            center_list = []
            boundaries_list = []
            binHalfWidth_list = []
            eff_list = []
            eff_err_list = []

            for bin, val in bins.items():

                bin_id.append(bin)
                x, x_err = PlotHelpers.GetBinCenterWithErr(val['boundaries'])
                boundaries_list.append(val['boundaries'])
                center_list.append(x)
                binHalfWidth_list.append(x_err)

                tmpVal = PlotHelpers.protectNegatives(val['efficiency'])

                eff_list.append(tmpVal[0])
                eff_err_list.append(
                    [tmpVal[2], tmpVal[1]]
                )  #Python is idiotic and needs the errors the other way around
                #low first, upper second

            #Flatten and simplify the bounaries list to only contain the bin edges
            x_edges = PlotHelpers.getBinEdges(boundaries_list)

            # Bin
            if ("-" in var):
                #In the case of a 2D variable, remove the 2D counting
                center_list = [
                    list(dict.fromkeys(np.array(center_list)[:, 0])),
                    list(dict.fromkeys(np.array(center_list)[:, 1]))
                ]

                #Overwrite the histogram based on the desired value
                if (var in set_values_dict.keys()):
                    for bin in set_values_dict[var]:
                        nBin = bin[0][1] + bin[0][0] * (len(x_edges[1]) - 1)
                        eff_list[nBin] = bin[1][0]
                        eff_err_list[nBin] = [bin[1][1], bin[1][2]]

            else:  #Overwrite the histogram based on the desired value
                if (var in set_values_dict.keys()):
                    for bin in set_values_dict[var]:
                        eff_list[bin[0][0]] = bin[1][0]
                        eff_err_list[bin[0][0]] = [bin[1][1], bin[1][2]]

            plotDict[var] = singlePlot(
                nBins=len(bin_id),
                center_list=center_list,
                binHalfWidth_list=binHalfWidth_list,
                x_edges=x_edges,
                eff_list=eff_list,
                eff_err_list=np.array(eff_err_list).T.tolist())
        return plotDict

    @staticmethod
    def isAllowed(listOfDicts, keyName, keyList):
        '''
            Checks whether there is T and Velo method, with the rest of the 
            dictionary being the same.
            If so, return the indices of the T and Velo samples.
        '''

        # List of allowed differences between the two samples in order to calculate
        # the Combined method, the Final method and the ratio between Data and MC
        diff_opts_allowed = ['binned_fit', 'sim_fit', 'fit_model']
        if (keyName == 'mode'): diff_opts_allowed.append('sim_ver')

        idx_list = []
        for i in range(0, len(listOfDicts)):
            for j in range(i + 1, len(listOfDicts)):
                d1_extra, d2_extra, diff, same = PlotHelpers.compareDict(
                    listOfDicts[i], listOfDicts[j])
                #OC.get_debug_text(f"d1_extra {d1_extra}")
                #OC.get_debug_text(f"d2_extra {d2_extra}")
                #OC.get_debug_text(f"diff {diff}")
                diff_opts = diff.keys()
                if ((set(diff_opts) - set(diff_opts_allowed)) != {keyName}):
                    continue
                if (sorted(diff[keyName]) == sorted(keyList)):
                    idx_list.append([i, j])
                else:
                    continue

        return idx_list


class Plotter():

    def __init__(self, listOfFiles, setBinDict_list):
        '''
            Initialize the plotter: save the lists of samples to be plotted 
        '''

        self._graphs = []
        self._varList = []

        for filename in listOfFiles:
            info_dict, eff_dict, set_values_dict = PlotReader.getInfoAndEfficiency(
                filename, setBinDict_list)
            plot_dict = PlotReader.loadEfficiency(eff_dict, set_values_dict)
            self._graphs.append({
                'info': info_dict,
                'plot': plot_dict
            })  # TODO: not really needed to have the info with the plots,
            # at some point remove

        # Now find only the common variables
        self._varList = set.intersection(
            *map(set, [d['plot'] for d in self._graphs]))

        # Make a list of all the infos
        self._allInfos = [d['info'] for d in self._graphs]

        # Convert the graphs into a numpy array for easier handling
        self._graphs = np.array(self._graphs)

        # If possible, add combined method
        for idx in PlotReader.isAllowed(self._allInfos, 'method',
                                        ["Velo", "T"]):
            self.addNew(idx, singlePlotUtils.getCombinedPlot, 'method',
                        "Combined")

        # If possible, add final method
        for idx in PlotReader.isAllowed(self._allInfos, 'method',
                                        ["Long", "Combined"]):
            self.addNew(idx, singlePlotUtils.getFinalPlot, 'method', "Final")

        # If possible, add ratio of Data/MC
        for idx in PlotReader.isAllowed(self._allInfos, 'mode',
                                        ["Data", "MC"]):
            self.addNew(idx, singlePlotUtils.getRatioPlot, 'mode', "Ratio")

        return

    def addNew(self, idx, function, keyName, newKeyValue):
        '''
            Function that saves the combined/final method, ratios and relative differences (useful for systematical studies).
            Input is the indices typically obtained by the isAllowed function,
            checking whether is is allowed to calculate the new efficiencies/ratios.
            The item corresponding to keyName will be overwritten with a newKeyValue.
        '''

        # Check whether the order for the ratio is correct (doesn't matter
        # for the combined and final methods)
        if (newKeyValue == "Ratio"
                and self._graphs[idx[0]]['info']['mode'] != "Data"):
            OC.get_info_text(
                f"Attempted to calculate the ratio of {self._graphs[idx[0]]['info']['mode']}/{self._graphs[idx[1]]['info']['mode']}. Flipping the order."
            )
            idx.reverse()
        tmpDict = {}
        for var in self._varList:
            newPlot = function(self._graphs[idx[0]]['plot'][var],
                               self._graphs[idx[1]]['plot'][var])
            tmpDict[var] = newPlot
        tmpInfo = dict(self._graphs[idx[0]]['info'])
        tmpInfo[keyName] = newKeyValue
        singlePlotUtils.savePlots(
            tmpDict, tmpInfo,
            Paths.getFitOutput(tmpInfo, tmpInfo['mode'], tmpInfo['method']))
        self._graphs = np.append(self._graphs, {
            'info': tmpInfo,
            'plot': tmpDict
        })

        self._allInfos.append(tmpInfo)
        return

    def createSampleTag(self, idx_list):

        # First make sure to use only the desired stuff
        info_list = [d['info'] for d in self._graphs[idx_list]]

        mergedDicts = (pd.DataFrame(info_list).groupby(
            ['year']).agg(set).reset_index().to_dict('records'))
        tagKeyList = ['year', 'mode']
        print(mergedDicts)
        for mergedDict in mergedDicts:  #It returns a list for different years
            for key, item in mergedDict.items():
                if (len(item) > 1): tagKeyList.append(key)
        tagKeyList = list(dict.fromkeys(tagKeyList))  #Remove duplicates

        return tagKeyList

    def make1DPlot(self, var, idx_list, tagList,
                   opts):  #opts needed for the file name

        OC.get_info_text(f"Plotting variable {var} in 1D.")

        #Check the colorlist is long enough for the required number of plots
        PlotDesign.checkNumberOfPlots(len(idx_list))

        #First make sure to plot only the desired stuff
        graph_list = np.array(self._graphs[idx_list])

        #Get the plt figure
        fig, axs = plt.subplots(figsize=(width, height))

        #Always include the method name
        methodName = PlotDesign.getMethodName(tagList, graph_list)
        modeName = []

        for i, graph in enumerate(graph_list):
            label = ""
            for item in tagList:
                if (item in graph['info']):
                    toAdd = graph['info'][item]
                    print(toAdd)
                    if (toAdd == False): continue
                    if (toAdd == True): toAdd = item
                    label = label + f"{toAdd} "
            #Check the maximal errors and set the values acordingly
            #yes, handling the numpy arrays/lists is all over the place now
            eff, eff_err = singlePlotUtils.checkError(
                np.array(graph['plot'][var].eff_list),
                np.array(graph['plot'][var].eff_err_list), opts.max_error_1D)

            #TODO: move the axis to PlotDesign
            axs.errorbar(graph['plot'][var].center_list,
                         eff,
                         xerr=graph['plot'][var].binHalfWidth_list,
                         yerr=np.abs(eff_err),
                         fmt="o",
                         color=colorlist[i],
                         label=label,
                         barsabove=True,
                         elinewidth=1.2,
                         capsize=3,
                         markersize=5)
            if ('method' in tagList):
                methodName.append(graph['info']['method'])
            modeName.append(graph['info']['mode'])
            PlotDesign.design1D(axs, methodName)
            PlotDesign.setLabels(Dictionaries.latex_names_dict[var],
                                 r"$\varepsilon$", axs)

        fig_name = Paths.getPlotPath(opts=opts,
                                     year=opts.year,
                                     method="_".join(methodName),
                                     tag="_".join(modeName),
                                     var=var)

        if (opts.verbose): OC.get_info_text(f"Saving into {fig_name}")
        fig.savefig(fig_name, bbox_inches='tight')
        plt.close()
        return

    def make2DPlot(self, var, idx_list, taglist, opts):
        '''
            Plots a 2D plot from singlePlot
        '''

        OC.get_info_text(f"Plotting variable {var} in 2D.")

        #First make sure to plot only the desired stuff
        graph_list = np.array(self._graphs[idx_list])

        var_x, var_y = var.split('-')

        for graph_dict in graph_list:
            graph = graph_dict['plot']
            fig, axs = plt.subplots(figsize=(width, height))
            PlotDesign.setLogScale(var_x, var_y)

            #Check the maximal errors and set the values acordingly
            eff, eff_err = singlePlotUtils.checkError(
                np.array(graph[var].eff_list),
                np.array(graph[var].eff_err_list), opts.max_error)

            if (graph_dict['info']['mode'] == "Ratio"):
                eff, eff_err = singlePlotUtils.checkDeviation(
                    np.array(graph[var].eff_list),
                    np.array(graph[var].eff_err_list), opts.max_deviation)

            xBins = len(graph[var].x_edges[0]) - 1
            eff_arr = PlotHelpers.get2Deff(eff, xBins)
            eff_err_arr = PlotHelpers.get2DeffErr(eff_err, xBins)

            im = axs.pcolormesh(graph[var].x_edges[0],
                                graph[var].x_edges[1],
                                eff_arr,
                                edgecolors='k',
                                linewidths=0.5,
                                cmap='rainbow',
                                vmin=0.85,
                                vmax=1.10)

            PlotDesign.addText2D(graph[var].center_list[0],
                                 graph[var].center_list[1], eff_arr,
                                 eff_err_arr)
            PlotDesign.setLabels(Dictionaries.latex_names_dict[var_x],
                                 Dictionaries.latex_names_dict[var_y], axs)
            colBar = fig.colorbar(im, pad=0.005)
            PlotDesign.design2D(axs, colBar, graph_dict['info']['method'],
                                graph_dict['info']['mode'], "efficiency")
            fig_name = Paths.getPlotPath(opts=opts,
                                         year = opts.year,
                                         method = graph_dict['info']['method'],
                                         tag = graph_dict['info']['mode']+\
                                               "_"+graph_dict['info']['sim_ver'],
                                         var =var)
            if (opts.verbose): OC.get_info_text(f"Saving into {fig_name}")
            fig.savefig(fig_name, bbox_inches='tight')
            plt.close()
        return

    def makePlot(self, idx_list, opts):
        PlotDesign.setFontFamily()
        tagList = self.createSampleTag(idx_list)
        for var in self._varList:
            if ("-" in var): self.make2DPlot(var, idx_list, tagList, opts)
            else: self.make1DPlot(var, idx_list, tagList, opts)

        return

    def PlotMCvsData(self, opts):
        OC.get_info_text("Plotting Data vs MC.")
        for idx in PlotReader.isAllowed([d['info'] for d in self._graphs],
                                        'mode', ["Data", "MC"]):
            idx.reverse()
            self.makePlot(idx, opts)
        return

    def PlotOnly(self, key, item, opts):
        for idx in range(0, len(self._graphs)):
            if (self._graphs[idx]['info'][key] == item):
                self.makePlot([idx], opts)
        return

    def PlotOnlyMethod(self, method, opts):
        self.PlotOnly('method', method, opts)
        return

    def PlotRatios(self, opts):
        OC.get_info_text("Plotting the ratios of efficiency in Data/MC.")
        self.PlotOnly("mode", "Ratio", opts)
        return

    def PlotMethodsComparison(self, methods, opts):
        for idx in PlotReader.isAllowed([d['info'] for d in self._graphs],
                                        'method', methods):
            self.makePlot(idx, opts)
        return

    def PlotLongFinalComparison(self, opts):
        self.PlotMethodsComparison(["Long", "Final"], opts)
        return


def Plot(opts):

    FileUtils.CreateFolder(Paths.getPlotFolder(opts), opts.verbose)

    if (opts.file_list) == None:
        file_list = []
        for method in opts.method:
            file_list.append(Paths.getFitOutput(opts, "MC", method))
            file_list.append(Paths.getFitOutput(opts, "Data", method))

        pltr = Plotter(file_list, opts.set_bins)
        pltr.PlotMCvsData(opts)
        pltr.PlotRatios(opts)

    else:
        for filenames in opts.file_list:
            if opts.set_bins == None:
                pltr = Plotter(filenames, Dictionaries.create_bin_dictionary("","",""))
            else:
                pltr = Plotter(filenames, opts.set_bins)
            idx_list = list(range(0, len(pltr._graphs)))
            pltr.makePlot(idx_list, opts)
    return

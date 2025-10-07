from Utilities import OutputColoring as OC
import Dictionaries
import numpy as np
from collections import namedtuple
import pickle as pkl
from Utils import binEdges
import Stats


#Named tuple containing the information needed to plot one file
singlePlot = namedtuple("singlePlot", 
                    "nBins " #Number of bins
                    +"center_list "  #List of bin centers
                    +"binHalfWidth_list " #List of half-widths of the bins
                    +"x_edges "  #List of bin edges
                    +"eff_list "  #Efficiency values
                    +"eff_err_list" #Errors                    
                    )
 
class singlePlotUtils: 
    def checkError(val, err, max_error):
        if int(max_error) == 0: return val, err
        mask_plus = err > max_error
        mask_minus = err < -max_error
        val[np.where(mask_minus[0]+mask_plus[1])] = Dictionaries.default_value
        err[np.where(mask_plus)] = Dictionaries.default_error
        err[np.where(mask_minus)] = -Dictionaries.default_error
        return val, err

    def checkDeviation(val, err, max_deviation):
        if int(max_deviation) == 0: return val, err
        mask = abs(val-Dictionaries.default_value) > max_deviation
        val[np.where(mask)] = Dictionaries.default_value
        err[:,np.where(mask)] = Dictionaries.default_error
        err[0] = -err[0]
        return val, err

    def checkBinning(plot1, plot2):
        '''
            Checks whether the binning in plot1 agrees with plot2. If yes, returns True.
        ''' 
        if not isinstance(plot1,singlePlot):
            OC.get_error_text("You cannot compare calss SinglePlot with anything else!")
            return NotImplemented 

        for key in ['center_list','binHalfWidth_list','x_edges']:
            if (getattr(plot1,key) != getattr(plot2,key)):                 
                OC.get_warning_text(f"Binning error: {key} does not agree!")
                OC.get_warning_text(f"{getattr(plot1,key)} != {getattr(plot2,key)}")
                return False
        return True

    def getNewPlot(plot1, newEff, newEffErr):
        '''
            Checks whether the binning in plot1 agrees with plot2. If yes, returns True.
        ''' 
        newPlot = plot1._replace(eff_list = newEff, eff_err_list = newEffErr)        
        return newPlot

 
    @staticmethod
    def getFinalPlot(plot_Long, plot_Comb):
        #Check the binning agrees
        if not singlePlotUtils.checkBinning(plot_Long, plot_Comb):
            raise Exception(OC.get_error_text("Your LONG and COMBINED plots have different binnings. Abort."))

        #Calculate the new efficiency and errors
        final_eff, final_eff_err = Stats.getFinalMethod(plot_Long.eff_list, plot_Long.eff_err_list, 
                                                                   plot_Comb.eff_list, plot_Comb.eff_err_list)
        #Return the new named tuple
        return singlePlotUtils.getNewPlot(plot_Long,final_eff,final_eff_err)


    @staticmethod
    def getCombinedPlot(plot_Velo, plot_T):
        #Check the binning agrees
        if not singlePlotUtils.checkBinning(plot_Velo, plot_T):
            raise Exception(OC.get_error_text("Your VELO and T plots have different binnings. Abort."))        

        #Calculate the new efficiency and errors
        comb_eff, comb_eff_err = Stats.getCombinedMethod(plot_Velo.eff_list, plot_Velo.eff_err_list, 
                                                                   plot_T.eff_list, plot_T.eff_err_list)

        #Return the new named tuple
        return singlePlotUtils.getNewPlot(plot_Velo,comb_eff,comb_eff_err)
    
    @staticmethod
    def getRatioPlot(plot_Data, plot_MC):
         #Check the binning agrees
        if not singlePlotUtils.checkBinning(plot_Data, plot_MC): 
            raise Exception(OC.get_error_text("Your DATA and MC plots have different binnings. Abort."))

        #Calculate the new efficiency and errors
        ratio_val, ratio_err = Stats.getRatio(plot_Data.eff_list, plot_Data.eff_err_list, 
                                                        plot_MC.eff_list, plot_MC.eff_err_list)

        #Return the new named tuple
        return singlePlotUtils.getNewPlot(plot_Data,ratio_val,ratio_err)



    @staticmethod
    def getRelativeDiffPlot(first, second):
         #Check the binning agrees
        if not singlePlotUtils.checkBinning(first, second): 
            raise Exception(OC.get_error_text("Your DATA and MC plots have different binnings. Abort."))

        #Calculate the new efficiency and errors
        diff_val, diff_err = Stats.getRelativeDiff(first.eff_list, first.eff_err_list, second.eff_list, second.eff_err_list)

        #Return the new named tuple
        return singlePlotUtils.getNewPlot(first,diff_val,diff_err)


    

    @staticmethod
    def savePlots(plot_dict, info_dict, output_path):
        '''
            Takes the plot dicitonary and the info dictionary and saves them into a pkl file
            defined by the output_path in the same format as the Fit.py
        '''
        outFile = open(output_path,"wb") 
        results_binned = {}
        for var, plot in plot_dict.items(): #This is stupid, I hate it, but it is what it is
            results_binned[var] = {}
            for i in range(0,plot.nBins):      #So I have to somehow recalculate the bin order      
                if ("-" in var):
                    x_len = len(plot.x_edges[0])-1                    
                    y_len = len(plot.x_edges[1])-1
                    j = i/y_len+(i%y_len)*x_len
                    j = int(j)
                else:
                    j = i                     
                results_binned[var][i] = {}             
                results_binned[var][i]["efficiency"] = [plot.eff_list[j], plot.eff_err_list[1][j], plot.eff_err_list[0][j]]
                results_binned[var][i]["boundaries"] = binEdges.getBoundaries(j,plot.x_edges)                
        pkl.dump(info_dict,outFile)
        pkl.dump(results_binned,outFile)
        outFile.close()
        return


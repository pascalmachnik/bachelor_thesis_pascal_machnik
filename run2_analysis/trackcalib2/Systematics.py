from Fit import Fit
from Paths import Paths, FileUtils
from Plot import Plotter
from PlotMaker import singlePlotUtils
from Utilities import OutputColoring as OC

class Systematics:
    
    def __init__(self, opts):
        self.options = opts
        return

    def fit_for_systematics(self, sysName, optionName, optionTags_list, optionValues_list):

        '''
            Performs the fits used to check the difference between samples with 
            different optionValues in the optionName. Saves into a subfolder with
            systematics/sysName/optionTags. 
        '''

        #Check the lenght of all the lists
        list_len = len(optionTags_list) 
        if (list_len != len(optionValues_list)):
            raise Exception(OC.get_error_text(f"Your input of optionTags {optionTags_list} has different lenght to optilValues {optionValues_list}. Abort."))
        

        for method in self.options.method:            
            for mode in self.options.mode:          
                for i in range(0, list_len):      
                    self.options.subfolder = f"systematics/{sysName}/{optionTags_list[i]}"
                    vars(self.options)[optionName] = optionValues_list[i]
                    Fit(method=method, mode=mode, opts=self.options)

        return 


    def plot_for_systematics(self, sysName, optionName, optionTags_list, optionValues_list):
        '''
            Makes the plots used to check the difference between samples with 
            different optionValues in the optionName. Saves into a subfolder with
            systematics/sysName/. It has to be performed separately to the fit due to the subparsers in argParse.
        '''

        file_list = []
        idx_list = []
        counter = 0 
        list_len = len(optionTags_list) 

        #For now, only use 2 different options
        if (list_len != 2):
            raise Exception(OC.get_error_text(f"Your input of optionTags {optionTags_list} does not contian exactly two values. Abort."))
        #Check the lenght of all the lists
        if (list_len != len(optionValues_list) ):
            raise Exception(OC.get_error_text(f"Your input of optionTags {optionTags_list} has different lenght to optilValues {optionValues_list}. Abort."))
        
        for method in self.options.method:            
            for mode in self.options.mode:
                self.options.subfolder = f"systematics/{sysName}/{optionTags_list[0]}"
                vars(self.options)[optionName] = optionValues_list[0]
                file_list.append(Paths.getFitOutput(self.options, mode, method)) 

                self.options.subfolder = f"systematics/{sysName}/{optionTags_list[1]}"
                vars(self.options)[optionName] = optionValues_list[1]
                file_list.append(Paths.getFitOutput(self.options, mode, method)) 
                
                idx_list.append([counter,counter+1])
                counter = counter +2
        
        self.options.subfolder = f"systematics/{sysName}/"
        FileUtils.CreateFolder(Paths.getPlotFolder(self.options), self.options.verbose)
        pltr = Plotter(file_list)
        for idx in idx_list:
            pltr.makePlot(idx,self.options)
            pltr.addNew(idx,singlePlotUtils.getRelativeDiffPlot, 'method', 'diff')

        pltr.PlotOnly('method','diff',self.options)
        return


    # Check the difference in effiicencies between fit with linear and with exponential background 
    def fit_linearVsExponential(self):
        self.fit_for_systematics("line_vs_expo","simple_bkg",["expo","line"], [False,True])
        return

    def plot_linearVsExponential(self):
        self.plot_for_systematics("line_vs_expo","simple_bkg",["expo","line"], [False,True])
        return
    
    #Check the difference in effiicencies between fit with double Gauss and double CB signal model 
    def fit_gaussVsCB(self):
        self.fit_for_systematics("gauss_vs_cb","simple_sig",["cb","gauss"], [False,True])
        return

    def plot_gaussVsCB(self):
        self.plot_for_systematics("gauss_vs_cb","simple_sig",["cb","gauss"], [False,True])
        return

    
    #Check the difference in effiicencies between magUp and magDown samples
    def fit_polarity(self):
        self.fit_for_systematics("polarity","polarity",["MagUp","MagDown"], ["MagUp","MagDown"])
        return

    def plot_polarity(self):
        self.plot_for_systematics("polarity","polarity",["MagUp","MagDown"], ["MagUp","MagDown"])
        return

    # Among others, for the final systematical uncertainty on the tables, the overlapVariation,
    # the motherEta and the double counting happening in the combined method should be taken into account
        
def runSystematics(opts):
       
    myErr  = Systematics(opts)
    if (opts.task == "fit"):    
        myErr.fit_linearVsExponential()
        myErr.fit_gaussVsCB()
        myErr.fit_polarity()
    elif (opts.task == "plot"): 
        myErr.plot_linearVsExponential()
        myErr.plot_gaussVsCB()
        myErr.plot_polarity()
    return
    
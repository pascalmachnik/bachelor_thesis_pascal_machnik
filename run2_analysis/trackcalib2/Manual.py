# Manual to be displayed by TrackCalib2
# Use option -man to get the manual 


from ArgParser import parse_args
from Utilities import OutputColoring as OC
import sys
from shutil import get_terminal_size
from collections import deque #For substracting a dictionary
from prettytable import PrettyTable
from textwrap import fill

references = {

    "TrackCalib Twiki": "https://twiki.cern.ch/twiki/bin/view/LHCb/TrackCalib",
    "TrackEff Twiki": "https://twiki.cern.ch/twiki/bin/view/LHCbInternal/LHCbTrackingEfficiencies",
    "Crystal Ball Wiki": "https://en.wikipedia.org/wiki/Crystal_Ball_function",
}

def printTasks():
    print (OC.TRKCLB + " is divided into three parts, executed separately using the following arguments:")

    print ("{0:10}{1}".format("prepare","Creates reduced tuple for each method with required variables. Provides weights to MC in multiplicity (in the given variable) and applies additional matching criteria when specified."))

    print ("{0:10}{1}".format("fit","Performs the fits to Jpsi mass using the previously created tuples. The mass is binned in the given set of variables."))

    print ("{0:10}{1}".format("plot","Plots the efficiencies vs the given variables, as obtained from the fit. It creates data/MC correction tables in all 2D variables."))
    return 

class Manual():

    def __init__(self,termSize):
        self._parent_dict, self._prepare_dict, self._fit_dict, self._plot_dict = parse_args(sys.argv[1:],True)
        #The parent dictionary is included in the prepare, fit and plot dictionaries, hence we want to get rid of them
        consume = deque(maxlen=0).extend
        consume(self._prepare_dict.pop(key,None) for key in self._parent_dict)
        consume(self._fit_dict.pop(key,None) for key in self._parent_dict)
        consume(self._plot_dict.pop(key,None) for key in self._parent_dict)
        #Remove the help option as it 
        self.terminal_size = termSize
        return
    
    def formatHelp(self,srt_Help):
        helpList = srt_Help.split("\n")
        return helpList
    
    def printArgDict(self,dict):        #Print the argument in a nice pretty table 
        table = PrettyTable(["Option","Default value","Usage"], columns_width = self.terminal_size/2, autowrap = True)
        for arg,details in dict.items():     
            helpList = self.formatHelp(details['help'])            
            table.add_row([', '.join(details['name']),details['default'],fill(helpList[0])])
            if (len(helpList)>0):
                for helpLine in helpList[1:]:
                    table.add_row(["","",fill(helpLine)])
        table.align["Option"] = "l"
        table.align["Default value"] = "c"
        table.align["Usage"] = "l"
        table.padding_width = 1
        print (table)
        return 
    
    def printGlobal(self):
        self.printArgDict(self._parent_dict)
        return
    
    def printPrepare(self):
        self.printArgDict(self._prepare_dict)
        return

    def printFit(self):
        self.printArgDict(self._fit_dict)
        return

    def printPlot(self):
        self.printArgDict(self._plot_dict)
        return

    @staticmethod
    def printReferences():
        for key,item in references.items():
            print ("{0:25} {1}".format(OC.bold_text(key),item))
        return        



def printManual():
    print ('''     

        _|      _|                                               _|  
        _|_|  _|_|     _|_|_|   _|_|_|     _|    _|     _|_|_|   _|  
        _|  _|  _|   _|    _|   _|    _|   _|    _|   _|    _|   _|  
        _|      _|   _|    _|   _|    _|   _|    _|   _|    _|   _|  
        _|      _|     _|_|_|   _|    _|     _|_|_|     _|_|_|   _|  
                                                                                    
    ''')
    terminal_size = get_terminal_size().columns
    print (f"{OC.DOT} {OC.TRKCLB} is a tool creating tracking efficiency correction tables. It allows for user-defined track-quality cuts, binning and variable combinations. It also allows checking the 1D and 2D efficiency dependencies for variables of choice. ")
    print (f"{OC.DOT} The name and the idea were shamelessly stolen from PIDCalib.")
    print (f"{OC.DOT} {OC.TRKCLB} is a successor of the beloved TrackCalib. ")
    print ("\n")
    printTasks()
    print ("\n Additional arguments can be specified for each of the three steps. The arguments that can (or have to be) specified for all three steps are: ")    
    manual = Manual(terminal_size)
    manual.printGlobal()
    print (f"\n The optional arguments applicable only to the {OC.bold_text('prepare')} part are:\n") 
    manual.printPrepare()
    print (f"\n The optional arguments applicable only to the {OC.bold_text('fit')} part are:\n") 
    manual.printFit()
    print (f"\n The optional arguments applicable only to the {OC.bold_text('plot')} part are:\n") 
    manual.printPlot()
    print ("\n\nFor more details and more information, see the references:\n ")
    manual.printReferences()    
    print("\n\nSidenote: If the table prints with linebreaks, increase the terminal width for better readability.")
    return   


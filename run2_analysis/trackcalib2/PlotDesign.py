import matplotlib.pyplot as plt
from Utilities import OutputColoring as OC
from Dictionaries import logVars

#Define default sizes

class Sizes:
    AXIS_LABEL_SIZE = 24
    AXIS_TICKLABEL_SIZE = 20
    LEGEND_SIZE = 16
    TITLE_SIZE = 20
    LHCBTAG_SIZE = 26
    METHOD_SIZE = 24
    VALUE_TEXT_SIZE = 16


colorlist = ["#000000", "#940D0D", #black, red
             "#00478F","#0A3900",  #blue, green
             "#9E8900","#A72564",  #yellow, magenta
             "#007B7D","#707070",  #teal, gray
             "#C76000","#6313B9",  #orange, purple
]

class PlotDesign:

    @staticmethod
    def checkNumberOfPlots(nPlots):
        try:
            colorlist[nPlots-1]
        except IndexError:
            OC.get_error_text("Too many graphs in one plot requested!")
            OC.get_error_text(f"The current possible maximal number of graphs is {len(colorlist)}.")
            OC.get_error_text("Either plot less graphs or extend the list of colors in PlotDesign.py.")
            exit()    

    @staticmethod
    def setFontFamily():
        plt.rcParams['text.usetex'] = True
        plt.rcParams["font.family"] = "Times New Roman" 


    @staticmethod
    def setLogScale(var_x,var_y):
        if (var_x in logVars): plt.xscale('log')
        if (var_y in logVars): plt.yscale('log') 
        pass

    @staticmethod
    def addLHCbtag(x_pos,y_pos,axis):
        plt.text(x_pos, y_pos, f"LHCb preliminary", transform = axis.transAxes, fontweight = "normal", fontsize = Sizes.LHCBTAG_SIZE)

    @staticmethod
    def getMethodName(tagList,graph_list):
        methodName = []
        if ('method' not in tagList): 
            methodName = [graph_list[0]['info']['method']]
        return methodName

    @staticmethod   
    def design1D(axis, method):           
        plt.setp(axis.get_xticklabels(), fontsize=Sizes.AXIS_TICKLABEL_SIZE)
        plt.setp(axis.get_yticklabels(), fontsize=Sizes.AXIS_TICKLABEL_SIZE)
        axis.set_title(f"", fontsize = Sizes.TITLE_SIZE)        
        plt.ylim(0.85, 1.05)
        PlotDesign.addLHCbtag(0.60,0.925,axis)
        if len(method) == 1:
            plt.text(0.05, 0.925, f"{method[0]} method",transform = axis.transAxes, 
                fontweight = "normal", fontsize = Sizes.METHOD_SIZE)
        plt.tick_params(axis = 'both', which = 'both')   
        plt.legend(prop={'size': Sizes.LEGEND_SIZE}, loc = 'lower right', 
                    frameon = True, edgecolor = "white", facecolor = "white")
        return

    @staticmethod    
    def addText2D(x_centers, y_centers, eff_arr, eff_err_arr): #I miss ROOT :(   
        for y in range(len(y_centers)):
            for x in range(len(x_centers)):                
                plt.text(x_centers[x], y_centers[y],
                        r"$\bf{"+str(f" {eff_arr[y,x]:.3f}")+r"}$"+f"\n\ $+$ {eff_err_arr[1][y,x]:.3f}\n $-$ {-eff_err_arr[0][y,x]:.3f}",
                        horizontalalignment='center', verticalalignment='center',
                        fontsize  = Sizes.VALUE_TEXT_SIZE, weight="bold"
                 ) 
        return

    @staticmethod     
    def design2D(axis, colorBar, method, mode, tag):        
        colorBar.ax.tick_params(labelsize=Sizes.AXIS_TICKLABEL_SIZE-2) 
        plt.setp(axis.get_xticklabels(), fontsize=Sizes.AXIS_TICKLABEL_SIZE)
        plt.setp(axis.get_yticklabels(), fontsize=Sizes.AXIS_TICKLABEL_SIZE)
        PlotDesign.addLHCbtag(0.70,1.05,axis)
        plt.text(-0.1, 1.05, f"{method} method {tag} {mode}", transform = axis.transAxes, 
                 fontweight = "normal", fontsize = Sizes.METHOD_SIZE)
        return
    

    @staticmethod
    def setLabels(x_label, y_label,axis):
        axis.set_xlabel(x_label, fontsize = Sizes.AXIS_LABEL_SIZE, loc = "right")
        axis.set_ylabel(y_label, fontsize = Sizes.AXIS_LABEL_SIZE, loc = "top")
        return
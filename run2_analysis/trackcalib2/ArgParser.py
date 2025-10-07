import argparse
from distutils.log import error
import sys
import os
from typing import Dict

import Dictionaries
from Utils import argument_dict
from Utilities import OutputColoring as OC
from numba import int64


class ShowArgumentsParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n\n' % message)
        self.print_usage(sys.stderr)
        try:
            sys.stderr.write('\n' + self.description)
        finally:
            print("no value inserted!")
        sys.exit(2)


def parse_args(args, forManual = False):

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        prog=os.path.basename(sys.argv[0]),
        description=("""
  Welcome to TrackCalib2. 

  For a full help, including the list of arguments, use examples and all available samples, do: 
            
            
      'python {0} manual'
            
  """).format(os.path.basename(sys.argv[0]))
    )

    # parent parser for common commands
    parent_parser = argparse.ArgumentParser(add_help=False) 

    # add the positional arguments
    parent_parser.add_argument("-v",
                               "--verbose",
                               default=False,
                               action="store_true",
                               help="Increase output verbosity.")

    parent_parser.add_argument("-systematics",
                               "--systematics",
                               default=False,
                               action="store_true",
                               help="Allows to run basic systematic tests. Use together with the fit option to perform the necessary fits and with plot option to create corresponding plots. It cannot be performed in one step!")

    subparser = parser.add_subparsers(title="tasks")
    parser_man = subparser.add_parser("manual", parents=[parent_parser], add_help=True) 
    parser_man.set_defaults(task = "manual")

    # prepare, fit, plot
    parent_parser.add_argument('-mode',
                               metavar='<mode>',
                               default=["MC","Data"],
                               action=CheckMode,
                               type=str,
                               help="""By using this argument you can choose between running TrackCalib2 only on Data or MC. If you desire to use both, leave this argument out. Use as -mode "MC" or  -mode "Data". """)

    # prepare, fit, plot
    parent_parser.add_argument('-year',
                               metavar='<year>',
                               required=True,
                               type = str,
                               action=CheckYear,
                               help=OC.WARNING +"Required argument!"+OC.ENDC+" Sets the year and conditions for given data-taking period. The full list of available samples is "+ '\n'.join(Dictionaries.year_dict.keys()))

    # prepare, fit, plot
    parent_parser.add_argument('-sim_ver',
                               metavar='<sim_ver>',
                               default=None,
                               action=CheckSimVersion,
                               type=str,
                               help="Sets the Sim version for MC. If none is selected, the default version for each year is automatically assigned. The available and default versions for each year are:\n"
                               + "\n".join("{0:20} {2:10} {1}".format(yr,'{}'.format(sim.sim_list),sim.sim_default) for yr,sim in Dictionaries.year_dict.items())
                               )

    # prepare, fit, plot
    parent_parser.add_argument('-method',
                               metavar='<method>',
                               default=Dictionaries.methods,
                               action=CheckMethod,
                               type=str,
                               help="Sets the tracking efficiency method. Choose any (combination) of " + ", ".join(Dictionaries.methods)+ ". The methods have to be separated by commas. If none is explicitly assigned, the default is using all methods.")

    # prepare, fit, plot
    parent_parser.add_argument('-variables',
                               metavar='<variables>',
                               default=list(Dictionaries.bin_dict_default.keys()),
                               type=str,
                               action=CheckVarList,
                               help="""Sets the variables versus which are the tracking efficiencies obtained. The variables have to be separated by commas. Example of use is "P,PT,ETA".""") 

    # fit, plot
    parent_parser.add_argument('-variables_2D', 
                               metavar='<variables_2D>',
                               default=Dictionaries.var_dict_2D_default,
                               type=str, 
                               action=CheckVarDict2D,
                               help=f"""Sets the 2D variables for which the efficiency dependency is evaluated. Minus connects the variables within a set, comma separates sets. An example of use is -variables_2D "PT-ETA, PT-P". If you do not want to create the ratios, run with -variables_2D "". """) #The default for the manual prinitng is specified in Utils, argument_dict

    # fit, plot
    parent_parser.add_argument('-polarity',
                               metavar='<polarity>',
                               default="",
                               action=CheckPolarity,
                               type=str,
                               help="Sets the mode to be split by polarity. Possible options are "+ ' or '.join(Dictionaries.polarities)+". If not used, both polarities are used.")

    # fit, plot
    parent_parser.add_argument('-subfolder',
                               metavar='<subfolder>',
                               default="",
                               type=str,
                               help="Set the name for the subfolder used to save the fit results. Especially helpful when performing tests or systematic studies. The fit results are then saved in results/year/subfolder/ and plot results are saved in plots/year/subfolder.")

    # prepare, fit
    parent_parser.add_argument('-max_entries',
                               metavar='<maxEntries>',
                               dest='max_entries',
                               default=-1,
                               type=int64,
                               help="Sets the maximum number of events used per method and charge.") 

    parent_parser.add_argument('-match_crit',
                                metavar='<match_crit>',
                                dest='match_crit',
                                default="",
                                type=str, #TODO: add check action
                                help="""Adds custom criteria to match the probe track (when is an event considered efficient).  Usage:\n-matchCrit 'Probe:Matched_GhostProb<0.3;Global:nSPDHits<400'
                                """)
    parent_parser.add_argument('-high_multi',
                                action=CheckMultiplicity,
                                nargs=0,
                                default = False,
                                help="""When using the stripped samples, events with higher nSPDHits are available. Use -high_multi to study these events. Available only in the stripped samples!
                                """)
                                #TODO: test it    

    parser_prepare = subparser.add_parser("prepare", parents=[parent_parser], add_help=False)
    parser_prepare.set_defaults(task = "prepare")

    parser_fit = subparser.add_parser("fit", parents=[parent_parser], add_help=False) 
    parser_fit.set_defaults(task = "fit")

    parser_plot = subparser.add_parser("plot", parents=[parent_parser], add_help=False) 
    parser_plot.set_defaults(task = "plot")
     
    # prepare
    parser_prepare.add_argument('-store_official', #FOR EXPERTS ONLY! IF YOU ARE NOT AN EXPERT, BE PREPARED FOR THE WRATH OF KHAN
                                action="store_true",
                                default=False,
                                help="FOR EXPERTS ONLY! IF YOU ARE NOT AN EXPERT, BE PREPARED FOR THE WRATH OF KHAN!") #TODO: test it

    # prepare
    parser_prepare.add_argument("-t",
                                "--test",
                                dest='test',
                                action="store_true",
                                default=False,
                                help="Only for experts! Unless you want to break things, please do not use this option.") 

    # prepare
    parser_prepare.add_argument("-f",
                                "--force",
                                action="store_true",
                                default=False,
                                help="Overwrite previously created datasets.")

     # prepare
    parser_prepare.add_argument('-weight_var',
                                metavar='<weight_var>',
                                dest='weight_var',
                                default="nSPDHits",
                                type = str,
                                help="Defines the variable used to obtain weights for MC.")

    # prepare
    parser_prepare.add_argument('-cuts',
                                metavar='<cuts>',
                                default="",
                                type=str,
                                help="""Sets the additional cuts to be applied to the tuple. The cuts can be applied on probe, tag, Jpsi or the global event. 
                                Adds custom criteria to match the probe track (when is an event considered efficient).  Usage:\n -cuts "Long:{Probe:Matched_GhostProb<0.3;Global:nSPDHits<400};Velo:{Probe:Matched_GhostProb<0.3}" 
                                The lists of available variables are:"""
                                +OC.BOLD+"\nJpsi:\n"+OC.ENDC
                                +", ".join(Dictionaries.var_list_mother)
                                +OC.BOLD+"\nTag/Probe:\n"+OC.ENDC
                                +", ".join(Dictionaries.var_list_part)
                                +OC.BOLD+"\nGlobal:\n"+OC.ENDC
                                +", ".join(Dictionaries.var_list_glob)
                                ) 
                                #TODO: possibly add the HLT decisions, but only possibly
 
    # fit
    parser_fit.add_argument('-official', 
                            action="store_true",
                            dest='official',
                            default=False,
                            help="Use official ntuples in the eos official location")

    # fit
    parser_fit.add_argument('-binning',
                            metavar='<binning>',
                            default="",
                            dest='binning',
                            type=str,
                            action = CheckBinning,
                            help=f"""Sets customized binning scheme for a given variable. The bin borders are separated by a comma, the variables are separated by a semicolon. Use example would be\n -binning "ETA:1.5,2.5,5;P:5000,50000,100000,200000". The default is {Dictionaries.bin_dict_default}. \nFor more options on binning, see also -auto_binning and -fine_binning.""" #The default printed by the manual is specified in Utils, argument_dict
                            )
    # fit
    parser_fit.add_argument('-just_one_bin',
                            metavar='<just_one_bin>',
                            default="",
                            dest='just_one_bin',
                            type=str,
                            help=f"""Allows you to plot just one bin in given variable(s). The bin counts starts from 0! The variables are separated by a semicolon. Use example would be\n -just_one_bin "ETA:0;P:4;ETA-P:1-3"."""
                            )

    # fit
    parser_fit.add_argument('-auto_binning',
                            metavar='<auto_binning>',
                            dest='auto_binning',
                            default="",
                            type=str,
                            help="""Sets equidistant binning scheme for a given variable in the format <var>:<Nbins><,Low>,<High>. The variables are separated by a semicolon. Use example would be\n -auto_binning "ETA:3,1,4;P:5,50000,175000" resulting in bins for ETA:1,2,3,4 and for P: 50000,75000,100000,125000,150000,175000. If not specified, the default binning will be used."""
                            )

    # fit
    parser_fit.add_argument('-fine_binning',
                            dest='fine_binning',
                            action="store_true",
                            default=False,
                            help="Sets the default binning scheme to be more fine binned. The used binning will be: "
                            +str(Dictionaries.bin_dict_fine)+". If not specified, the default binning scheme will be used."
                           )

    # fit
    parser_fit.add_argument('-simple_sig',
                            dest='simple_sig',
                            action="store_true",
                            default=False,
                            help="Use simplified Jpsi mass signal model. When False, the fit shape Double Crystal Ball (for more details see references) for signal. When using this argument, the fit shape is two Gaussians for signal. This option should be turned on for samples with low statistics.")
    # fit
    parser_fit.add_argument('-simple_bkg',
                            dest='simple_bkg',
                            action="store_true",
                            default=False,
                            help="Use simplified Jpsi mass background model. When False, the fit shape is an exponential. When using this argument, the fit shape is set to a linear polynomial. This option should be turned on for samples with low statistics.") 
    # fit
    parser_fit.add_argument('-binned_fit',
                            dest='binned_fit',
                            action="store_true",
                            default=False,
                            help="Use binned fit. It is advised to use this option only for samples with enough statistics.") 

    # fit
    parser_fit.add_argument('-bins_in_fit',
                            metavar='<binsInFit>',
                            dest='bins_in_fit',
                            default=100,
                            type = int,
                            help="Sets the number of bins in the binned fit (= number of bins in the Jpsi mass). ")

    # fit
    parser_fit.add_argument('-sim_fit',
                            dest='sim_fit',
                            action="store_true",
                            default=False, 
                            help="Perform a simultaneous fit the the matched and failed events. When not used, the parameter performs two independent fits to the matched and failed events. This is recommended for small samples. For larger samples, use this argument in order to perform a simultaneous fit to all and matched events.")

    # fit
    parser_fit.add_argument('-use_weights',
                            dest='use_weights',
                            action="store_true",
                            default=False,
                            help="Use multiplicity weights in the fit (to be implemented!)") #TODO

    # fit
    parser_fit.add_argument('-input_filename',
                            dest='input_filename',
                            default="",
                            type=str,
                            help="(FOR EXPERT USER ONLY) Use a specific file for testing")
    # fit
    parser_fit.add_argument('-no_full_fit',
                            dest='no_full_fit',
                            action="store_true",
                            default=False,
                            help="Don't fit the full dataset (no binning applied, use the whole dataset in the fit to Jpsi).")
    
    # fit
    parser_fit.add_argument('-muon_charge',
                            dest='muon_charge',
                            default=0,
                            type=int,
                            help="Define which muon charge is used for the fit. 0 means both charges are used, 1 means only positive muon charge is used and -1 means only negative muon charge is used. The default is 0.")
    
    # plot
    parser_plot.add_argument('-max_deviation',
                        metavar='<max_deviation>',
                        default=0.03,
                        type = float,
                        help=f"Sets the maximum deviation of efficiency RATIO from one. If the error is larger than this limit, the corresponding efficiency ratio bin value is set to {Dictionaries.default_value}{OC.PM}{Dictionaries.default_error}. For no limit set to zero.")

    # plot
    parser_plot.add_argument('-max_error',
                             metavar='<max_error>',
                             default=0.03,
                             type = float,
                             help=f"Sets the maximum possible error for tracking efficiency RATIOS. If the error is larger than this limit, the corresponding efficiency ratio bin value is set to {Dictionaries.default_value}{OC.PM}{Dictionaries.default_error}. For no limit set to zero.")

    # plot
    parser_plot.add_argument('-max_error_1D',
                             metavar='<max_error_1D>',
                             default=0.1,
                             type = float,
                             help=f"Sets the maximum possible error for tracking efficiencies. If the error is larger than this limit, the corresponding efficiency ratio bin value is set to {Dictionaries.default_value}{OC.PM}{Dictionaries.default_error}. For no limit set to zero.")
    
    # plot
    parser_plot.add_argument('-set_bins',
                        metavar='<set_bins>',
                        default="",
                        type=str,
                        action=SetBinValue,
                        help=f"""Sets the given efficiency bins to given value and error. It is useful when you have low statistics in only one bin and hence want to use a value from a dedicated fit. Use as file_name:var:bin:value+plus uncertainty-uncertainty, separate variables by semicolon ;. Example of usage: 'results/2018_25ns/trackEff_Data_Long.pkl:P:0:0.99+0.1-0.2;results/2018_25ns/Sim09h/trackEff_MC_T.pkl:ETA:1:0.95+0.2-0.15' """
                        )

    # plot
    parser_plot.add_argument('-file_list',
                        metavar='<file_list>',
                        default=None,
                        type=str,
                        action = CheckPlotFiles,
                        help="List of pickle (.pkl) files that should be plotted in one plot. All variables stored in the pickle file will be plotted. The paths are separated by the comma. In the case of producing multiple plots, separate the paths by ';'. The default is to plot MC and Data efficiency in the given year: \n result/year/trackEff_Data_Long.pkl, results/year/sim_ver/trackEff_MC_Long.pkl; result/year/trackEff_Data_T.pkl, results/year/sim_ver/trackEff_MC_T.pkl; result/year/trackEff_Data_Velo.pkl, results/year/sim_ver/trackEff_MC_Velo.pkl; result/year/trackEff_Data_Combined.pkl, results/year/sim_ver/trackEff_MC_Combined.pkl; result/year/trackEff_Data_Final.pkl, results/year/sim_ver/trackEff_MC_Final.pkl") 

    
    if (forManual):
        parent_dict = argument_dict(parent_parser) #vars(parent_parser.parse_args(args))
        prepare_dict = argument_dict(parser_prepare)
        fit_dict = argument_dict(parser_fit)
        plot_dict = argument_dict(parser_plot)
        return parent_dict, prepare_dict, fit_dict, plot_dict

    options = parser.parse_args(args)
    
    #Skipping setting up simVer and WGproduct if task is manual
    if (options.task == "manual"): return options

    #Don't do full fit when fitting only one bin
    if (options.task == "fit" and options.just_one_bin != ""): options.no_full_fit = True
    
    # Set the WGProduct tag to True except for 2015
    # Legacy issue: in 2016, the plot production was done automatically and it was needed
    #               to test whether 2016 local and WG productions agree
    #               So only 2015 was not GP, the rest was, hence set to False if 2015, else True
    options.WGProduct = False if ('2015' in options.year) else True

    #Set default SimVersion if it is not already set
    if (options.sim_ver == None):
        options.sim_ver = Dictionaries.year_dict[options.year].sim_default
   

    return options


class CheckYear(argparse.Action):
    """Check if the year selected is among the available years"""
    def __call__(self, parser, namespace, values, option_string=None):
        #If manual to be printed, do not require the year to be properly initialized
        pass
        if values in Dictionaries.year_dict.keys():
            setattr(namespace, self.dest, values)
        else:
            #Year not found, printout the error message
            OC.get_error_text(f"The year {values} is not available!")
            OC.get_info_text("Available years are:")
            OC.get_info_text("\n".join(Dictionaries.year_dict.keys()))
            OC.get_error_text("Abort.")
            sys.exit(1)


class CheckSimVersion(argparse.Action):
    """Check if the simulation version is among all the available simulation"""
    def __call__(self, parser, namespace, values, option_string=None):
        if values in Dictionaries.year_dict[namespace.year].sim_list:
            if (namespace.verbose): print("Setting the simulation version to", values)
            setattr(namespace, self.dest, values)
        else:
            OC.get_error_text("%s not present!"%values)
            OC.get_info_text("Available MC versions are:")
            OC.get_info_text("\n".join(Dictionaries.year_dict[namespace.year].sim_list))
            sys.exit(1)



class CheckMethod(argparse.Action):
    """Check if the method/s selected is/are among the available methods"""
    def __call__(self, parser, namespace, values, option_string=None):
        values = values.split(",")
        __temp_values = []
        for value in values:
            if value in Dictionaries.methods:
                __temp_values.append(value)
            else:
                parser.error(OC.get_error_text("Invalid method %s given!"%value))
        if len(__temp_values) == 0:
            parser.error(OC.get_error_text("A problem appeared with the methods selected"))
        else:
            if (namespace.verbose): OC.get_info_text(f"Setting the methods to {__temp_values}")
            setattr(namespace, self.dest, __temp_values)
        return


class CheckVarList(argparse.Action):
    """Check if the 1D variables are sensible"""
    def __call__(self, parser, namespace, values, option_string=None):
        
        OC.get_warning_text("Custom cuts OVERWRITE the default cuts, they are not added!")

        if (values == ""):
            setattr(namespace, self.dest, [''])
            return

        #Split the input string
        var_list = values.split(",")
        for var in var_list:
            if var not in Dictionaries.var_list_part + Dictionaries.var_list_mother + Dictionaries.var_list_glob:                               
                parser.error(OC.get_error_text(f"Invalid variable {var} given."))
        if namespace.verbose: OC.get_info_text("Custom set of variables requested: "+ ','.join(var_list))
        setattr(namespace, self.dest, var_list)
        return

class CheckVarDict2D(argparse.Action):
    """Check if the 2D variables are sensible"""
    def __call__(self, parser, namespace, values, option_string=None):
        
        var_dict_2D = {}
        if (values == ""):
            OC.get_info_text("Turning off the 2D ratios calculation.")
            setattr(namespace, self.dest, {})
            return
        for string_2D in values.split(','):
            tmpvar = string_2D.split('-')
            #Check the formatting
            if (len(tmpvar)!=2): 
                parser.error(OC.get_error_text(f"""Invalid formatting of the 2D variable! In case of producing two tables, first one VAR1 vs VAR2 and the second one VAR2 vs VAR3, it has to be formatted as follows:\n -variables_2D "VAR1-VAR2,VAR3-VAR4".\n For more help see the README or run trackcalib.py manual."""))
            #Check the variables are meaningful
            for var in tmpvar:
                if var not in Dictionaries.var_list_part + Dictionaries.var_list_mother + Dictionaries.var_list_glob:
                    parser.error(OC.get_error_text(f"Invalid variable {var} given in -variables_2D."))
            #Check the variables are not the same
            if (tmpvar[0]==tmpvar[1]): 
                parser.error(OC.get_error_text(f"""Creating the efficiency ratio tables of VAR1 vs VAR2 is not meaningful. For creating the tables in eg. PT and ETA run with -variables_2D "PT-ETA". """))

            #If all good, pass it to the dictionary
            var_dict_2D[string_2D] = [tmpvar[0], tmpvar[1]]
  
        if namespace.verbose: OC.get_info_text("Custom set of 2D variables requested: " + ','.join(var_dict_2D))
        
        setattr(namespace, self.dest, var_dict_2D)
        return

class CheckPolarity(argparse.Action):
    """Check if the polarity selected is among the available polarities"""
    def __call__(self, parser, namespace, values, option_string=None):
        if values in Dictionaries.polarities:
            if (namespace.verbose): print ("Running only on %s polarity." % values)
            setattr(namespace, self.dest, values)
        else:
            parser.error(OC.get_error_text('Invalid polarity %s given, must be "MagUp" or "MagDown"!' % values))

class CheckBinning(argparse.Action):
    """Check if the custom binning has increasing bin borders"""
    def __call__(self, parser, namespace, values, option_string=None):
        tmpDict  = Dictionaries.create_bin_dictionary(values, "", False, namespace.verbose)        
        for var,bins in tmpDict.items():                    
            for i in range(len(bins)-1):
                if (bins[i]>=bins[i+1]): 
                    parser.error("Wrong binning for variable %s! The bin edges %d %d are not in an increasing order!" % (var,bins[i],bins[i+1])) #TODO: test                    
        setattr(namespace, self.dest, values)
        return

class CheckMode(argparse.Action):
    """Check if the mode/s selected is/are among the available modes"""
    def __call__(self, parser, namespace, values, option_string=None):        
        values = values.split(",")
        __temp_values = []
        if len(values) != 0:
            for value in values:
                if value in Dictionaries.modes:
                    __temp_values.append(value)
                else:
                    parser.error(OC.get_error_text('Invalid mode %s given, must be "Data" or "MC"!' % values))
            if len(__temp_values) == 0:
                parser.error(OC.get_error_text("A problem appeared with the modes selected"))
            else:
                setattr(namespace, self.dest, __temp_values)
        else:
            setattr(namespace, self.dest, Dictionaries.modes)


class CheckMultiplicity(argparse.Action):
    """Check if turning the high multiplicity study on makes sense"""
    def __call__(self, parser, namespace, values, option_string=None):        
        if "strip" in namespace.year.split("_"):
            setattr(namespace, self.dest, True)            
        else:
            OC.get_warning_text("High multiplicity studies are not available for Turbo samples! Setting -high_multi back to False.")
            setattr(namespace, self.dest, False)
        return


class CheckPlotFiles(argparse.Action):
    """Check if the selected files exist and make a list out of them"""
    def __call__(self, parser, namespace, values, option_string=None):
        file_list = []
        #Create the 2D list
        for lst in values.split(";"):
            file_list.append(lst.split(','))
        #Check if all the members are existing
        for i in range(len(file_list)):
            for file in file_list[i]:
                #If the file does not exist, drop this list of files (and continue with the rest of the plots)
                if not (os.path.exists(file)):
                    OC.get_error_text(f"File {file} does not exist!")
                    OC.get_error_text(f"Skipping plotting the files {file_list[i]}")
                    file_list.pop(i)
                    break            
        setattr(namespace, self.dest, file_list)


class SetBinValue(argparse.Action):
    """Parse the input string to be able to easily set value for given bin"""
    def __call__(self, parser, namespace, values, option_string=None):   
        info_dict_list = [] 
        full_list = values.split(";") 
        for info in full_list:
            info_list = info.split(":") 
            if (len(info_list) != 4):
                parser.error(OC.get_error_text(f"Something went wrong when trying to parse {info}! \n Use the format file:var:bin:eff+err-err"))
                
            #Now add the errors, check which one is upper and which one is lowe
            plus = info_list[3].find("+")
            minus = info_list[3].find("-")
            val = []
            if (plus > minus):            
                val = [float(info_list[3][0:minus]),
                       -float(info_list[3][minus+1:plus]),
                       float(info_list[3][plus+1:])
                    ]
            else:     
                val = [float(info_list[3][0:plus]),
                       -float(info_list[3][minus+1:]),
                       float(info_list[3][plus+1:minus])
                    ]
            bin = [int(b) for b in info_list[2].split("-")]

            info_dict_list.append({
                "file": info_list[0],
                "var": info_list[1],
                "bin": bin,
                "value": val,
            })            
        
        print(info_dict_list)
        setattr(namespace, self.dest, info_dict_list)


def check_file(options):
    """(FOR EXPERT USER ONLY) Check if the file used has some specific filename structure"""
    pass


import os
from XRootD import client

from pathlib import Path
from Utilities import OutputColoring as OC
from itertools import zip_longest

lfn_dict = {
    "Data" : {
        "2011_50ns_strip": {"MagDown": "00115408",
                            "MagUp": "00115398",
                        },
        "2012_50ns_strip": {"MagDown": "00115396",
                            "MagUp": "00115402",
                            },
        "2016_25ns_strip": {"MagDown": "00115404",
                            "MagUp": "00115406",
                            },
        "2017_25ns_strip": {"MagDown": "00115410",
                            "MagUp": "00115400",
                            },
        "2016_25ns": {"MagDown": "00066675",
                      "MagUp": "00066676",
                     },
        "2017_25ns": {"MagDown": "00079274",
                      "MagUp": "00079273",
                     },
        "2018_25ns": {"MagDown": "00083388",
                      "MagUp": "00083387",
                     }
    },
    "MC": {        
        "2011_50ns_strip": {"MagDown": "00115604",
                            "MagUp": "00115596",
                            },
        "2012_50ns_strip": {"MagDown": "00115594",
                            "MagUp": "00115600",
                            },
        "2016_25ns_strip": {"MagDown": "00115598",
                            "MagUp": "00115602",
                            },
        "2017_25ns_strip": {"MagDown": "00096208",
                            "MagUp": "00096209",
                            },
        "2017_25ns": {"MagDown": "00096208",
                      "MagUp": "00096209",
                     },
        "2018_25ns": {"MagDown": "00096206",
                      "MagUp": "00096207",
                     }     
        } #2016 MC is not in WG
}

class FileUtils():
    '''
        Class handling moving files around, creating folders, etc        
    '''
    
    @staticmethod
    def CreateFolder(path,verbose):
        '''
        Helper function to create folders enriched by proper info messages
        '''
        if not os.path.exists(path):
            if verbose: OC.get_info_text("Creating new folder for plots: "+ path)
            os.makedirs(path)
        else:
            if (verbose): OC.get_info_text(f"Folder {path} already exists")
        return

    def CheckIfCanOverwrite(filePath, force):
        '''
            Checks whether file already exists. If force is True, the existing file will be overwritten
        '''
        if os.path.isfile(filePath) and force is False:
            OC.get_error_text(f"File {filePath} exists.\nPlease move/delete this file or set --force/-f flag")
            return False    
            sys.exit(1)
        elif os.path.isfile(filePath) and force is True:
            OC.get_warning_text(f"File {filePath} exists.\n File will be overwritten!")
            return True
        else:
            OC.get_info_text(f"File {filePath} will be created")
            return True



class Paths():
    '''
    Class to get the correct paths for used files
    '''
    #Set eos locations
    
    #Location for the default prepared tuples
    EOS_LOC_OUT = "/eos/lhcb/wg/TrackingAlignment/TrackCalib/ToolTuples/" 

    #Location of the input files for the prepare.py    
    EOS_LOC_IN = "/eos/lhcb/wg/TrackingAlignment/WGP/"

    #Root tag
    EOS_ROOT = "root://eoslhcb.cern.ch/" 
    
    #Define the folders where to find the prepared tuples, fit files and the plots
    FOLDER_DICT = {
    "prepare": "./tuples",
    "fit": "./results",
    "plot": "./plots",
    }

    #Main filename
    MAIN_NAME = "trackEff"
    
    @staticmethod
    def interleave(*iterables):
        '''
        Zips together multiple iterables, continuing with the remaining elements if one iterable is shorter than the others.
        e.g. interleave([1, 2, 3, 4], [5, 6]) -> [1, 5, 2, 6, 3, 4]. This is used to interleave the MagDown and MagUp files.
        In case of low requested number of entries, there will be both MagDown and MagUp entries prepared.
        '''
        # Use zip_longest to handle iterables of different lengths
        zipped = zip_longest(*iterables, fillvalue=None)
        # Flatten the zipped pairs, filtering out the fillvalue
        return (item for group in zipped for item in group if item is not None)
    
    @staticmethod
    def getReweigtingFilePath(year): 
        return Paths.EOS_ROOT + Paths.EOS_LOC_OUT\
                + f"Data/{year}/Tuple_Data.root"        

    @staticmethod
    def getPrepareInputExtension(year,mode):     
        extension = "_1.trkcalib.root"
        info_year = year.split("_")
        # check if the is the ntuples coming from the stripping or turbo
        if "strip" in info_year:
            extension = "_1.tuple_" + info_year[0]
            if (mode == "Data"):
                if (info_year[0] != "2018"): extension += "_data_stripping.root" 
            else:
                if (info_year[0] == "2017" or info_year[0] == "2018"):
                    raise SystemExit(OC.get_error_text("No available MC for 2017 and 2018 stripping. Sorry, abort."))
                else: extension += "_mc_{1}_{0}.root"
        return extension

    @staticmethod
    def getPrepareInputData(eos_path, year, polarity, extension):

        location = f"{eos_path}Data_{year}_{polarity}/"                  
        filename_mask = lfn_dict["Data"][year][polarity] + f"_*{extension}"        
        
        return Path(location).glob(filename_mask)

    @staticmethod
    def getPrepareInputMC(eos_path,year, polarity, sim_ver,extension):

        location = f"{eos_path}MC/{year}/{sim_ver}/{polarity}/"
        info_year = year.split("_")
        if ("strip" in info_year) and info_year[0] in ["2016", "2017", "2018"]:
            location = f"{eos_path}MC/{year}/{sim_ver}/{polarity}/"        

        extension = extension.format(sim_ver.lower(), polarity[3:].lower())  
        filename_mask = lfn_dict["MC"][year][polarity] + f"_*{extension}"

        return Path(location).glob(filename_mask)

    @staticmethod
    def getPrepareInput(year, mode, sim_ver, wg_production, test):     
        '''   
            Returns the stripped/triggered ROOT files
        '''

        #If test, just quickly return a location
        if (False):#(test): 
            location = "../tuples/"
            filename_mask = "*_*_1.trkcalib.root"
            files = [str(file) for file in Path(location).glob(filename_mask)]
            OC.get_debug_text("Loading only small test sample...")
            return files
 
        #Get the extension based on whether strip/Turbo and Data/MC
        extension = Paths.getPrepareInputExtension(year,mode)

        if (mode == "Data"):
            if (wg_production):
                paths = Paths.interleave(
                    Paths.getPrepareInputData(Paths.EOS_LOC_IN,year,"MagDown",extension),
                    Paths.getPrepareInputData(Paths.EOS_LOC_IN,year,"MagUp",extension)
                )
            else:
                paths = [f"{ Paths.EOS_LOC_OUT}{mode}/{year}/Tuple_{mode}_MagDown.root",
                         f"{ Paths.EOS_LOC_OUT}{mode}/{year}/Tuple_{mode}_MagUp.root"]        
        else: #if mode=="MC"
            eos_path = Paths.EOS_LOC_OUT #Different to Data!
            if (wg_production and year != '2016_25ns'):
                paths = Paths.interleave(
                    Paths.getPrepareInputMC(eos_path,year,"MagDown",sim_ver,extension),
                    Paths.getPrepareInputMC(eos_path,year,"MagUp",sim_ver,extension)
                )     
            else:
                paths = [f"{ Paths.EOS_LOC_OUT}{mode}/{year}/{sim_ver}/Tuple_{mode}_MagDown.root",
                         f"{ Paths.EOS_LOC_OUT}{mode}/{year}/{sim_ver}/Tuple_{mode}_MagUp.root"]

        return [str(file) for file in paths]

    @staticmethod
    def getPrepareOutputFolder(official,verbose):
        '''
            Returns the folder where the prepare part output is stored.
        '''
        if (official):
            path = Paths.EOS_LOC_OUT +  Paths.FOLDER_DICT["prepare"]
        else:
            path = Paths.FOLDER_DICT["prepare"]

        FileUtils.CreateFolder(path,verbose)
       
        return path

    @staticmethod
    def getPrepareOutput(mode, method, year, **kwargs):
        '''
            Returns the path to the output of the prepare part. It is used by the fit as input files.
        '''
      
        verbose = False
        if ("verbose" in kwargs.keys()): verbose = kwargs["verbose"]
      
        match_crit = ""
        if "match_crit" in kwargs.keys(): match_crit = kwargs['match_crit']
      
        sim_ver = ""
        if "sim_ver" in kwargs.keys(): sim_ver = kwargs["sim_ver"]
      
        WGProduct = False    
        if ("wg_production" in kwargs.keys()): WGProduct = kwargs["wg_production"]

        official = False
        if ("official" in kwargs.keys()): official = kwargs["official"]

        high_multi = False
        if ("high_multi" in kwargs.keys()): high_multi = kwargs["high_multi"]
        
        path = Paths.getPrepareOutputFolder(official,verbose)

        path =  f"{path}/trackEffTuple_{mode}_{year}_{method}" \
                + ((f"_{sim_ver}") if mode == "MC" else "") \
                + ("_tight" if match_crit != "" else "") \
                + ("_high_mult" if high_multi else "") \
                + ("_WG" if WGProduct else "") \
                + ".root"
                #Be careful with the brackets, they are crucial!
        return path

    @staticmethod        
    def getFitOutputFolder(opts,mode):
        '''
            Returns the folder where the fit output files are stored
        '''
        if (type(opts) != dict):  #This allows the user to pass either options or a dedicated dictionary
            #Careful here, it should be argparse.Namespace, buuuut if something goes wrong, it can be anything
            opts_dict = vars(opts)
        else: opts_dict = opts
        if ('match_crit' not in opts_dict.keys()):
            opts_dict['match_crit'] = ""

        path = Paths.FOLDER_DICT["fit"] + "/" \
                + opts_dict['year'] \
                + (("_" + opts_dict['polarity']) if opts_dict['polarity'] != "" else "") \
                + ("_tight" if opts_dict['match_crit'] != "" else "") \
                + "/" + opts_dict['subfolder'] + "/" \
                + ((opts_dict['sim_ver']) if mode == "MC" else "")  \
                +"/" 
        # This protects double slash in case subfolder is empty,
        # but also protects the case when user puts 'subfolder/' as an option
        path = path.replace("//","/")
        return path


    @staticmethod
    def getStatusFile(opts, mode, method):
        '''
            Returns the path to the fit status file, where the basic fit information is saved
        '''
        path = Paths.getFitOutputFolder(opts,mode)\
                + f"{Paths.MAIN_NAME}_{mode}_{method}_fitStatus.out"
        return path

    @staticmethod
    def getWarningsFile(opts, mode, method):
        '''
            Returns the path to the fit warnings file, where important warnings and errors are stored. Existence of such a file is beneficial for debugging: the warnings/errors are concentrated in one place for quick problem identification.
        '''
        path = Paths.getFitOutputFolder(opts,mode) \
                + f"{Paths.MAIN_NAME}_{mode}_{method}_WARNINGS.out"
        return path
            

    @staticmethod        
    def getFitOutput(opts, mode, method):
        path = Paths.getFitOutputFolder(opts,mode)  \
             + f"{Paths.MAIN_NAME}_{mode}_{method}.pkl"
        return path

    @staticmethod        
    def getFitPlotFolder(opts,mode):
        '''
            Returns the folder where the fit control plots are stored
        '''        
        path = Paths.getFitOutputFolder(opts,mode) + "/plots/"
        #Create the folder if it is not there
        FileUtils.CreateFolder(path,opts.verbose)
        return path

    @staticmethod       
    def getFitPlotName(mode, method, opts, var = "", bin = -1):
        path = Paths.getFitPlotFolder(opts,mode)
       
        if type(var) == list: #2D binning
            if (len(bin) != len(var)):
                OC.get_error_text(f"Number of bins {len(bin)} does not agee with the number of variables {len(var)}!")
            path = path + f"{mode}_{method}_{var[0]}_bin_{bin[0]}_{var[1]}_bin_{bin[1]}"      
        else: #1D binning     
            #Else is not super safe here, but just for the path should be fine               
            if (var == ""): #Full dataset, no binning in any var
                path = path + f"{mode}_{method}_full_sample"
            else:
                if (bin==-1):
                    OC.get_error_text(f"Wrong variable {var} or bin {bin}!")
                path = path + f"{mode}_{method}_{var}_{bin}"    
        return path

    @staticmethod        
    def getPlotFolder(opts):
        '''
            Returns the folder where the final plots are stored
        '''   
        path = Paths.FOLDER_DICT["plot"] + "/" \
                + opts.year \
                + (("_" + opts.polarity) if opts.polarity != "" else "") \
                + ("_tight" if opts.match_crit != "" else "") \
                +"/" \
                + opts.subfolder + "/"
        # This protects double slash in case subfolder is empty,
        # but also protects the case when user puts 'subfolder/' as an option
        path = path.replace("//","/")
        
        return path


    @staticmethod        
    def getPlotPath(opts, year, method, tag, var):

        path = Paths.getPlotFolder(opts)  \
             + f"{Paths.MAIN_NAME}_{tag}_{year}_{method}_{var}.png"
        return path


    @staticmethod        
    def getPlotOutput(opts):
        path = ""
        return path

def adapt_path(url_path):
    if url_path.startswith("root://"):
        if os.path.exists("/eos") and "LBAP_TOKENS_FILE" not in os.environ:
            xrdurl = client.URL(url_path)
            url_path = xrdurl.path
    return url_path
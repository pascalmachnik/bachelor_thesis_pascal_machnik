import numpy as np
import gc
import uproot
import dask
import pickle as pkl
from typing import Union, List
import resource
from collections import namedtuple

from lib.src.TCFit.basic_pdfs import Gaussian, Exponential, CrystalBall, Polynomial
from lib.src.TCFit.model import Model, AddPDF
from lib.src.TCFit.fit import LikelihoodFit
from lib.src.TCFit.plot import Plot
from lib.src.TCFit.parameter import Parameter, ParameterDatabase
from lib.src.TCFit.formula import Formula
from lib.src.TCFit.stats import chi_square

from Utils import multiple_lines_title,\
    standard_text,\
    binEdges
from Utilities import OutputColoring as OC
from Paths import FileUtils, Paths, adapt_path
from fitMonitorFiles import monitorFiles
from Stats import calc_eff_and_error

from Sampling import array_select, matching, subsample
from Dictionaries import fit_config, \
    create_bin_dictionary, \
    one_bin_bin_dict, \
    create_var_list, \
    make_list_branches 


#set a named tuple for passing the tags for the fit control plots
pltTag = namedtuple("pltTag", "tag mode method muon_charge")
  
class fitInfo():
    def __init__(self,  **kwargs) -> None:        
        for key, value in kwargs.items():
            setattr(self, key, value)
        pass

    def dumpFitInfo(self,filename):
        pkl.dump(self.__dict__,filename)        
        

def using(point=""):
    usage = resource.getrusage(resource.RUSAGE_SELF)
    return '''%s: usertime=%s systime=%s mem=%s mb
           ''' % (point, usage[0], usage[1],
                  (usage[2] * resource.getpagesize()) / 1000000.0)

def define_model(method,
                 year: str,
                 simple_sig: bool = False,
                 simple_bkg: bool = False,
                 sample_size: int = 1000,
                 sim_fit: bool = True):
    """Define the mass model for the invariant mass fit.
    The invariant mass range depends on the method used and the year considered.
    The signal PDF can be described as sum of gaussians (simple_sig == True) or
    a sum of Crystal Balls (simple_sig == False).
    The background PDF can be described as a linear linr (simple_bkg == True) or
    as an exponential (simple_bkg == False).
    This function can also build the simultaneous pdf for fitting matched and
    unmatched events simultaneously (sim_fit == True).
    sample_size gives a rough estimante of signal and background yields constraints.
    Also returns the name of the signal and background models.
    """

    total_shape_names = {
        'signal' : "",
        'background': ""
    }
    
    # set invariant mass as observable; the mass range depends on the method used for the efficiency evaluation
    obs = fit_config["obs"]# mass = Parameter(obs['name'], limits=obs['bounds'][method], unit="MeV/c$^{2}$") 
    mass = Parameter(obs['name'], limits=obs['bounds'][method], unit="MeV/c$^{2}$", latex = r"m$_{\mu^{+}\mu^{-}}$")
    #TODO: if there is time, add the latex and unit parameters to the rest of the parameters too
    # Then also modify the plotting in the TCFit  #TODO
    if method == "T" and (year == "2011" or year == "2012"):
        mass.limit = obs['bounds']["T_Run1"]

    # set signal and background yields parameters
    signal_yield = Parameter("signal_yield",
                             value=0.5 * sample_size,
                             limits=(0.0, sample_size))
    background_yield = Parameter("background_yield",
                                 value=0.5 * sample_size,
                                 limits=(0.0, sample_size))

    # Signal parameters
    mean = Parameter("mean", value=3094.0, limits=(3090.0, 3115.0))
    sigma = Parameter("sigma", value=25.0, limits=(12.0, 45.0))
    sigma2 = Parameter("sigma2", value=70.0, limits=(20.0, 150.0))
    alpha = Parameter("alpha", value=1.0, limits=(0.1, 15.0))
    n = Parameter("n", value=1.0, limits=(1.0, 35.0))
    fraction = Parameter("frac", value=0.8, limits=(0.5, 1.))
    
    # Background parameters
    tau = Parameter("tau", value=0.0005, limits=(-0.005, 0.01))
    c0 = Parameter("c0", 1, (0,10000))
    c1 = Parameter("c1", -1, (-100,0))

    if simple_sig:
        gauss1 = Gaussian("mySingleGauss1", mass, mean, sigma)
        gauss2 = Gaussian("mySingleGauss2", mass, mean, sigma2)
        signal = AddPDF("myGaussSum", [gauss1, gauss2], [fraction])
        total_shape_names['signal'] = "GaussSum"
    else:
        crystal_ball1 = CrystalBall("myCB", mass, mean, sigma, alpha, n)
        crystal_ball2 = CrystalBall("myCB2", mass, mean, sigma2, alpha, n)
        signal = AddPDF("myCBsum",
                        pdfs=[crystal_ball1, crystal_ball2],
                        coefficients=[fraction])
        total_shape_names['signal'] = "CBsum"
    if simple_bkg:
        background = Polynomial("myPoly",mass,[c0, c1]) 
    else:
        background = Exponential("myExpo", mass, tau)

    # Total PDF as an extended sum od signal and background PDFs
    sum_pdfs = AddPDF('total_shape', pdfs=[signal, background],
                      coefficients=[signal_yield, background_yield])
    total_shape_names['background'] = "Expo"                 
    total_shape = Model("totalShape", pdfs=[sum_pdfs])

    # if sim_fit == True (simultaneous fit on matched and unmatched categories)
    # additional shape for unmatched events is defined
    if sim_fit:
        # some additional definitions for simultaneous fit
        sigma2_f = Parameter("sigma2F", value=70.0, limits=(20.0, 150.0))
        fraction_f = Parameter("fracF", value=0.9, limits=(0.1, 1.))
        
        if simple_sig:           

            gauss1_f = Gaussian("mySingleGauss1F", mass, mean, sigma)
            gauss2_f = Gaussian("mySingleGauss2F", mass, mean, sigma2_f)
            signal_f = AddPDF("myGaussSumF",
                              pdfs=[gauss1_f, gauss2_f],
                              coefficients=[fraction_f])
            total_shape_names['signal'] = "GaussSum"

        else:
            tau_f = Parameter("tauF", value=0.0005, limits=(-0.005, 0.01))
            crystal_ball1_f = CrystalBall("myCBF", mass, mean, sigma, alpha, n)
            crystal_ball2_f = CrystalBall("myCB2F", mass, mean, sigma2_f, alpha, n)
            signal_f = AddPDF("myCBsumF",
                              pdfs=[crystal_ball1_f, crystal_ball2_f],
                              coefficients=[fraction_f])
            total_shape_names['signal'] = "CBsum"
        
        if simple_bkg:
            c0_f = Parameter("c0F", 1, (0,10000))
            c1_f = Parameter("c1F", -1, (-100,0))
            background_f = Polynomial("myPoly",mass,[c0_f, c1_f])
        else:
            background_f = Exponential("myExpoF", mass, tau_f)


        # in case of simultaneous fit efficiency is the common parameter
        efficiency_sig = Parameter("efficiency_sig", value = 0.9, limits=(0., 1.))
        efficiency_bkg = Parameter("efficiency_bkg", value = 0.9, limits=(0., 1.))
        signal_yield_pass = Formula("f_pass", "@0*@1", 
                                    parameters=[signal_yield, efficiency_sig])
        signal_yield_fail = Formula("f_fail", "@0*(1.-@1)",
                                    parameters=[signal_yield, efficiency_sig])
        background_yield_pass = Formula("f_pass_bkg", "@0*@1",
                                        parameters=[background_yield,
                                                    efficiency_bkg])
        background_yield_fail = Formula("f_fail_bkg", "@0*(1.-@1)",
                                        parameters=[background_yield,
                                                    efficiency_bkg])

        total_shape_pass = AddPDF('total_shape_pass',
                                  pdfs=[signal, background],
                                  coefficients=[signal_yield_pass,
                                                background_yield_pass])
        total_shape_fail = AddPDF('total_shape_fail',
                                  pdfs=[signal_f, background_f],
                                  coefficients=[signal_yield_fail,
                                                background_yield_fail])
        total_shape = Model("total_shape_sim",
                            pdfs=[total_shape_pass, total_shape_fail])
    return total_shape, total_shape_names

def do_fit(data: Union[np.array, List[np.array]],
           model: Model,
           parameters_to_reset: dict = {},
           parameters_to_fix: list = [],
           fix_values: bool = False,
           plot_name: str = "default.pdf",
           binned_fit: bool = False,
           plotTag: namedtuple = pltTag(None,None,None,None)):
    """Perform fit on data using model PDF.
    If parameters_to_reset is not an empty dictionary, parameters are reset to their standard values.
    parameters_to_fix: list of parameters that have to be fixed or not depending on fix_values value.
    plot_name: filename to be used to store the fit plot output.
    nbins: the number of bins to be used in the output plot."""

    # fix or unfix parameters in parameters_to_fix list
    if len(parameters_to_fix) != 0:
        for par in parameters_to_fix:
            if par in model.related:
                model[par].fixed = fix_values

    # reset parameters in parameters_to_reset to the values defined in the dictionary
    if parameters_to_reset != {}:
        for par, val in parameters_to_reset.items():
            if par in model.related:
                model[par].value = val

    # perform the unbinned likelihood fit
    options = {}
    if binned_fit:
        options['bins'] = fit_config["nbins"]
        options['print_level'] = 2
    fit = LikelihoodFit(data=data, model=model, **options)

    do_minos = False if binned_fit else True
    minuit = fit.minimise(hesse=True, minos=do_minos)

    # plot fit result
    plot = Plot(minuit=minuit,
                model=model,
                data=data,
                limit=model.observable.limits,
                n_bins=fit_config['nbins'])

    plot.plot(binnedFit = binned_fit,
                tag = plotTag.tag,
                mode = plotTag.mode,
                method = plotTag.method,
                muon_charge=plotTag.muon_charge,
             )

    plot.flush(plot_name)
    return minuit

def evaluate_efficiency(data: Union[list, np.array],
                        model: Model,
                        sim_fit: bool,
                        binned_fit: bool,
                        name: str, 
                        outputFiles: monitorFiles,
                        is_data: bool,
                        plotTag,
                        bin: Union[int, list]) -> list:
    """evaluate the efficiency fitting the invariant mass shape.
    data: list of 3 np.arrays
        data[0] = matched candidates
        data[1] = unmatched candidates
    model: mass model PDF
    sim_fit: True if it is a simultaneous fit
    status_file: fitting status output file
    warning_file: warnings output file
    is_data: True if fit performed on data, False in case of MC
    bin: integer number defining in which variable range the fit is performed
       no bins: [-1,-1]
       1D: [<bin>,-1]
       2D: [<binx>,<biny>]
    """
    
    __binned_fit = binned_fit

    # data[0] = matched candidates
    # data[1] = unmatched candidates

    # prepare string for binning
    string_bin = ""

    if bin[0] != -1:
        if bin[1] != -1:
            string_bin = " bin : (%d,%d)\t " % (bin[0], bin[1])
        else:
            string_bin = " bin : %d\t " % (bin[0])

    # values used to reset model values
    list_params = ["mean", "sigma", "sigma2", "frac", "alpha", "n"]
    dict_params = {"sigma": 30.,
                   "sigma2": 70.,
                   "sigma2F": 70.,
                   "frac": 0.5,
                   "fracF": 0.5,
                   "alpha": 1.,
                   "n": 1.,
                   "tau": -.0005,
                   "tauF": -0.005}
    # thresholds for warnings
    # Set the limit for chi2/NDOF. If above this value, warnings are printed and saved into warningFile
    
    result = [-1, -1, -1]
    matched_candidates = len(data[0])
    unmatched_candidates = len(data[1])

    OC.get_info_text(f"Number of matched candidates   :: {matched_candidates}")
    OC.get_info_text(f"Number of unmatched candidates :: {unmatched_candidates}")

    # if the dataset is large enough switch to binned_fit
    if matched_candidates > 1e4 and unmatched_candidates > 1e4:
        OC.get_warning_text("Number of matched and unmatched candidates above threshold")
        OC.get_warning_text("Fit forced to be a binned maximum likelihood fit")
        __binned_fit = True

    if matched_candidates == 0 or unmatched_candidates == 0:
        OC.get_error_text("Dataset Empty! Efficiency set to default value -1")
        return result #TODO test and fix this

    if not sim_fit:
        fraction = 0.2 if is_data else 0.9
        model["signal_yield"].value = fraction * matched_candidates
        model["background_yield"].value = (1 - fraction) * matched_candidates

        # define plot name
        plot_name = name + "_matched.pdf" 
        plotTag = plotTag._replace(tag=["Matched"])

        # perform fit for matched events and save plot
        minuit_m = do_fit(data=[data[0]],
                          model=model,
                          parameters_to_reset=dict_params,
                          parameters_to_fix=list_params,
                          fix_values=False,
                          plot_name=plot_name,
                          binned_fit=__binned_fit,
                          plotTag = plotTag)

        
        outputFiles.writeFitStatus("Matched", string_bin, minuit_m.fmin.is_valid)

        __contents, __bins = np.histogram(data[0], fit_config['nbins'], model.observable.limits)
        chi2 = chi_square(model.pdfs[0], __bins, __contents)
        outputFiles.writeCHI2("Matched",chi2)

        outputFiles.writeFitResult(minuit=minuit_m,
                                   string_name="MATCHED",
                                   binx=bin[0],
                                   biny=bin[1])

        matched_sig_yield = model["signal_yield"].value #not used anywhere, but keep it in case
        matched_bkg_yield = model["background_yield"].value

        # fit to all candidates
        all_candidates = matched_candidates + unmatched_candidates

        model["signal_yield"].value = 0.3 * all_candidates
        model["background_yield"].value = 0.7 * all_candidates

        # define plot name
        plot_name = name + "_full.pdf"
        plotTag = plotTag._replace(tag=["Full"])

        # perform fit for all candidates and save plot
        minuit_t = do_fit(data=[np.concatenate([data[0], data[1]])],
                          model=model,
                          parameters_to_reset={},
                          parameters_to_fix=list_params,
                          fix_values=True,
                          plot_name=plot_name,
                          binned_fit=__binned_fit,
                          plotTag = plotTag)

        status = minuit_t.fmin.is_valid
        outputFiles.writeFitStatus("All",string_bin,status)

        __contents, __bins = np.histogram(np.concatenate([data[0], data[1]]),
                                          fit_config['nbins'],
                                          model.observable.limits)
        chi2 = chi_square(model.pdfs[0], __bins, __contents)
        outputFiles.writeCHI2("All",chi2) 

        outputFiles.writeFitResult(minuit=minuit_t,
                                   string_name="All",
                                   binx=bin[0],
                                   biny=bin[1])

        # get the fit results
        all_sig_yield = model["signal_yield"].value #not used anywhere, but keep it in case
        all_bkg_yield = model["background_yield"].value

        efficiency, hi_error, lo_error = calc_eff_and_error(all_candidates, matched_candidates, all_bkg_yield,
                                                            matched_bkg_yield)
        result = [efficiency, hi_error, -lo_error]


    else:  # Simultaneous fit

        # expected fraction
        fraction = 0.5 if is_data else 0.99

        # set initial values for the signal yield
        model["signal_yield"].value = fraction * (matched_candidates + unmatched_candidates)
        model["background_yield"].value = (1. - fraction) * (matched_candidates + unmatched_candidates)

        plot_name = name + "_sim.pdf"        
        plotTag = plotTag._replace(tag=["Matched","Failed"])
        # reset efficiencies
        #        dict_params["efficiency_bkg"] = 0.5
        #        dict_params["efficiency_sig"] = 0.5
        # perform fit
        minuit = do_fit(data=[data[0], data[1]],
                        model=model,
                        parameters_to_reset=dict_params,
                        parameters_to_fix=list_params,
                        fix_values=False,
                        plot_name=plot_name,
                        binned_fit=__binned_fit,
                        plotTag = plotTag)

        status =  minuit.fmin.is_valid
        outputFiles.writeFitStatus("Simultanous Fit",string_bin,status)

        OC.get_info_text(f"Efficiency: {str(model['efficiency_sig'])}")

        __contents, __bins = np.histogram(data[0], fit_config['nbins'], model.observable.limits)
        chi2 = chi_square(model.pdfs[0], __bins, __contents)
        outputFiles.writeCHI2("Matched",chi2)

        __contents, __bins = np.histogram(data[1], fit_config['nbins'], model.observable.limits)
        chi2 = chi_square(model.pdfs[1], __bins, __contents)
        outputFiles.writeCHI2("Failed",chi2)

        outputFiles.writeFitResult(minuit=minuit,
                                   string_name="SIMULTANEOUS",
                                   binx=bin[0],
                                   biny=bin[1])

        efficiency = model["efficiency_sig"].value
        if hasattr(model['efficiency_sig'], "asym_error"):
            hi_error = model["efficiency_sig"].asym_error[0]
            lo_error = -model["efficiency_sig"].asym_error[1]
        else:
            hi_error = model["efficiency_sig"].error
            lo_error = -model["efficiency_sig"].error

        result = [efficiency, hi_error, lo_error]        

    outputFiles.writeEfficiency(result)
    return result


def getFitCuts(obs,method,polarity):
    cuts = "(Mother_ETA > 0) "
    # add cut on mass ranges
    cuts += "& (" + obs['name'] + " > " + str(obs['bounds'][method][0])
    cuts += ") & (" + obs['name'] + "<= " + str(obs['bounds'][method][1]) + ")"
    if (polarity == "MagDown"):
        cuts += " & (Polarity == -1)"
    elif (polarity == "MagUp"):
        cuts += " & (Polarity == 1)"
    return cuts


def Fit(mode, method, opts):
    # start printout
    multiple_lines_title(f'Welcome to TrackEff Fitter',
                         f'Running on {mode}',
                         f'Using {method} method')

    # set the mode
    is_data = True if mode == "Data" else False

    # check if binned fit is required
    binned_fit = opts.binned_fit
    
    # turn the requested muon charges into a list
    if int(opts.muon_charge) == 0:
        charges = [-1, 1]
    else:
        charges = [int(opts.muon_charge)]
        
    charge_label = {
        -1: "Minus",
         1: "Plus"
    }

    # create results path if not existing
    FileUtils.CreateFolder(Paths.getFitOutputFolder(opts,mode), opts.verbose)      

    #Open the status and warning files    
    outFiles = monitorFiles(opts,mode,method)
    
    # create the default variable/bin border dictionary
    bin_dict = create_bin_dictionary(opts.binning, opts.auto_binning, opts.fine_binning, opts.verbose)
    OC.get_debug_text(f"bin_dict {bin_dict}")
    bin_dict = one_bin_bin_dict(bin_dict,opts.just_one_bin)
    OC.get_debug_text(f"bin_dict {bin_dict}")

    list_of_vars = create_var_list(opts.variables, opts.variables_2D, bin_dict, opts.verbose)
    OC.get_debug_text(f"list_of_vars {list_of_vars}")

    # set observable
    obs = fit_config["obs"]
    OC.get_info_text("Observable :: %s" % obs['name'])
    OC.get_info_text("Boundaries considered :: (%f,%f) MeV/c^2" % (obs['bounds'][method][0], obs['bounds'][method][1]))

    branches = make_list_branches(list_of_vars)
    branches.append(obs["name"])
    branches.append("Mother_ETA")
    branches.append("Polarity")
    branches.append("matched")
    if opts.use_weights:
        branches.append("weight")

    if opts.verbose:
        OC.get_info_text("Branches used: %s"%", ".join(branches))

    ########################END OF ARGUMENT PARSING#########################

    OC.get_info_text(using("Reading in the prepared files"))
    # read input trees
    input_filename = Paths.getPrepareOutput(mode=mode,
                                            method=method,
                                            year=opts.year,
                                            official=opts.official,
                                            wg_production=opts.WGProduct,
                                            sim_ver=opts.sim_ver,
                                            match_crit=opts.match_crit,
                                            verbose=opts.verbose)

    # only for expert: if we want to inject a file with custom name and location
    if opts.input_filename != "":
        input_filename = opts.input_filename

    # Load the prepared trees
    try:
        OC.get_info_text(f"Using {input_filename} as an input for fitting.")
        # consider only the TTree for the specific method and transform into lazyarrays
        cache = uproot.cache.LRUCache("1GB")
        trees_orig = { charge : uproot.dask(
                          f"{adapt_path(input_filename)}:TrackEffTree{charge_label[charge]}{method}", 
                          filter_name=branches, cache = cache
                          ) for charge in charges }
    except Exception as ex:
        OC.get_error_text(f"Input file {input_filename} is not opened! Abort.")
        OC.get_error_text(ex)
        exit()

    # define additional selection on prepared ntuples
    cuts = getFitCuts(obs,method,opts.polarity)
    if opts.verbose:
        OC.get_info_text("Applying additional cut on trees: " + cuts)

    # Final dictionary with a dataframe per charge
    trees = {}
    actual_entries = sum([np.shape(trees_orig[charge].compute())[0] for charge in charges])

    OC.get_info_text("Sample Size before generic cuts: " + str(actual_entries))
    OC.get_info_text("Maximum number of entries requested: " + str(opts.max_entries))
    if opts.verbose:
        entries_to_be_processed = opts.max_entries if (
                opts.max_entries < actual_entries and opts.max_entries > 0) else actual_entries
        OC.get_info_text("Copying " + str((entries_to_be_processed) / 1.e6) + " million entries per charge.")
        OC.get_info_text(f"Old tuples contain {actual_entries} entries.")

    for charge in charges:
        if (int(opts.max_entries) > 0 and int(opts.max_entries) < actual_entries):
            selection = array_select(trees_orig[charge].compute(),cuts) 
            trees[charge] = (trees_orig[charge].compute()[selection])[:int(opts.max_entries / 2)]
        else:
            selection = array_select(trees_orig[charge].compute(),cuts) 
            trees[charge] = trees_orig[charge].compute()[selection]

    del trees_orig
    gc.collect()
    sample_size = sum([np.shape(trees[charge])[0] for charge in charges])
    OC.get_info_text(f"Sample size after generic cuts: {sample_size}")

    if opts.verbose:
        for charge in charges:
            OC.get_info_text(f"{charge_label[charge]} sample size: {len(trees[charge])}")

    OC.get_info_text(f"Using variables {list_of_vars}")
    
    # create mass pdf
    model, model_names = define_model(method=method,
                                     year=opts.year,
                                     simple_sig=opts.simple_sig,
                                     simple_bkg=opts.simple_bkg,
                                     sample_size=sample_size,
                                     sim_fit=opts.sim_fit
                                     )

    OC.get_info_text(f"Model created. Simultaneous fit set to {opts.sim_fit}")
    if binned_fit:
        OC.get_warning_text("Binned likelihood fit!")

    #Set the tag used for the control plots
    plotTag = pltTag("",mode,method,opts.muon_charge)

    # use make_datesets to produce 3 unbinned datasets: matched, unmatched and total.
    datasets = [None] * 3

    
    ###################################################
    #
    #                  Perform fits
    #
    ###################################################

    # Fit on the whole sample first
    # Start printout
    outFiles.writeIntro(opts)
   
    if not opts.no_full_fit:
        # charges = [-1, 1] or [opts.muon_charge]
        m = {charge: matching(trees[ charge ]) for charge in charges}

        datasets[0] = np.concatenate([trees[ charge ][ m[ charge ]][obs['name']].to_numpy() for charge in charges])
        datasets[1] = np.concatenate([trees[ charge ][~m[ charge ]][obs['name']].to_numpy() for charge in charges])

        if opts.verbose: OC.get_info_text("Full sample created!")
        
        outFiles.writeBinInfo("")

        result_full_sample = evaluate_efficiency(data=datasets,
                                                model=model,
                                                sim_fit=opts.sim_fit,
                                                binned_fit=binned_fit,
                                                name=Paths.getFitPlotName(mode,method,opts),
                                                outputFiles=outFiles,
                                                is_data=is_data,
                                                bin=[-1, -1],
                                                plotTag=plotTag)
        # create outputfile and dump result
        out_file = open(Paths.getFitOutput(opts,mode,method), "wb")
        pkl.dump(result_full_sample, out_file)
        info = fitInfo(year=opts.year,mode=mode,method=method,
                       sim_fit=opts.sim_fit,binned_fit=binned_fit,
                       polarity=opts.polarity,
                       match_crit=opts.match_crit,   
                       subfolder=opts.subfolder,
                       fit_model = model_names['signal']+"_"+model_names['background']
                       )
        info.dumpFitInfo(out_file)
        out_file.close()
        gc.collect()
        
    ######################## 
    # Binned part
    ######################## NEW

    results_binned = {}
    # loop over the variable used to bin the datasample

    for var in list_of_vars:

        # check if the split is performed on 1 or 2 variables
        lvar = var.split("-")
        mode_2D = False
        if len(lvar) == 2: mode_2D = True

        total_number_of_bins = 0
        results_binned[var] = {}
        if mode_2D:
            OC.get_info_text(f'2D efficiencies requested in {var}!')
            OC.get_info_text(f"Starting fits in {lvar[0]} and {lvar[1]} using " + str(
                len(bin_dict[lvar[0]]) - 1) + "x" + str(len(bin_dict[lvar[1]]) - 1) + " bins.")
            bounds = [None] * 2
            bounds[0] = np.asarray(bin_dict[lvar[0]])
            bounds[1] = np.asarray(bin_dict[lvar[1]])
            total_number_of_bins = (len(bounds[0]) - 1) * (len(bounds[1]) - 1)
            OC.get_debug_text(f"bounds[0] {bounds[0]}")
            OC.get_debug_text(f"bounds[1] {bounds[1]}")


        else:
            OC.get_info_text(f"Starting fits in {var} using {len(bin_dict[var]) - 1} bins.")
            bounds = np.asarray(bin_dict[var])
            total_number_of_bins = len(bounds) - 1

        out_file = open(Paths.getFitOutput(opts, mode, method), "wb")

        info = fitInfo(year=opts.year,
                       mode=mode,
                       method=method,
                       sim_ver=opts.sim_ver,
                       sim_fit=opts.sim_fit,
                       binned_fit=binned_fit,
                       polarity=opts.polarity,
                       match_crit=opts.match_crit,   
                       subfolder=opts.subfolder,                       
                       fit_model = model_names['signal']+"_"+model_names['background']
                       ) 
        info.dumpFitInfo(out_file)        

        outFiles.writeBinInfo(var)

        # loop over all bins (it gets bins from the dictionary of datasets)
        # bin is a global bin in 2D case
        for ibin in range(total_number_of_bins):
            OC.get_info_text(standard_text['bold_line'])
            OC.get_info_text(f"Started fit in {var}")
            OC.get_info_text(f"bin {ibin}")

            results_binned[var][ibin] = {}
            # prepare for 2D efficiencies
            bin_x = ibin
            bin_y = -1
            if mode_2D:
                #Get the boundaries from the bin edges
                bin_y, bin_x = binEdges.getXYBin(ibin,[bin_dict[lvar[1]],bin_dict[lvar[0]]])
                #I could use binEdges.getXYBin and then reverse, but I need bin_x and bin_y anyhow
                results_binned[var][ibin]["boundaries"] = ([bin_dict[lvar[0]][bin_x],
                                                            bin_dict[lvar[0]][bin_x + 1]],
                                                           [bin_dict[lvar[1]][bin_y],
                                                            bin_dict[lvar[1]][bin_y + 1]])
                OC.get_debug_text("results_binned[var][ibin]['boundaries']")                  
                OC.get_debug_text(results_binned[var][ibin]['boundaries'])    
            else:
                results_binned[var][ibin]["boundaries"] = ([bin_dict[var][ibin],
                                                            bin_dict[var][ibin + 1]])

            # make an array of datasets
            # store the datasets with matched [0] and unmatched [1] candidates:
            datasets[0], datasets[1] = subsample(df=trees,
                                                 bin_dict=bin_dict,
                                                 variables=lvar,
                                                 obs=obs['name'],
                                                 index=[bin_x, bin_y],
                                                 charge=opts.muon_charge)  # charge == 0 means both charges

            if mode_2D:
                OC.get_info_text(f"Subdataset for ({bin_x},{bin_y}) in {var} cretated.\n")
            else:
                OC.get_info_text(f"Subdataset for {ibin} in {var} cretated.\n")
            
            # WRT previous version: naming changed bins starts from 0 not from 1            
            if (mode_2D):
                outputPlotName = Paths.getFitPlotName(mode,method,opts,lvar,[bin_x,bin_y])
            else:
                outputPlotName = Paths.getFitPlotName(mode,method,opts,var,ibin)

            results_binned[var][ibin]["efficiency"] = evaluate_efficiency(data=datasets,
                                                                          model=model,
                                                                          sim_fit=opts.sim_fit,
                                                                          binned_fit=binned_fit,
                                                                          name=outputPlotName,
                                                                          outputFiles=outFiles,
                                                                          is_data=is_data,
                                                                          bin=[bin_x,
                                                                          bin_y],
                                                                          plotTag=plotTag)


            gc.collect()

            if mode_2D:
                OC.get_info_text(f"Fit in {var} bin ({bin_x},{bin_y}) done!\n")
            else:
                OC.get_info_text(f"Fit in {var} bin {ibin} done!\n")


        OC.get_ok_text(f"Fits done in {var}")
        pkl.dump(results_binned, out_file)
        out_file.close()
        gc.collect()

    #Close the warning and status files
    outFiles.close()

    gc.collect()
    # Delete the registry of Parameters used for the fit
    ParameterDatabase.parameter_registry = {}
    OC.get_ok_text(f"All fits for {method} method are done!")

import numpy as np
from numba import njit, vectorize, float64, int64
from math import exp, gamma, sqrt
from scipy.stats import norm, binomtest
from typing import Tuple

from Utilities import OutputColoring as OC

@vectorize([float64(int64, float64)], fastmath=True)
def poisson(k, mu):
    """Poisson probability"""
    #Not using the one from scipy, cause njit won't compile non-common scipy modules
    if mu < 0 or k < 0:
        return 0
    return mu ** k * exp(-mu) / gamma(k + 1)

@vectorize([float64(int64, float64)], fastmath = True)
def poisson_cdf(k,mu):
    "Poisson comulative"
    #Not using the one from scipy, cause njit won't compile non-common scipy modules
    prob = 0.0
    if k < 0:
        return 0
    for x in np.arange(0,k+1):
        prob += poisson(x,mu)
    return prob


@njit(fastmath=True, parallel=True)
def check_belt(mu: float64, nobs: int64, nbkg: float64, cl: float64) -> bool:
    """Feldman-Cousin Belt construction with Poissonian probability"""
    size = 50*(1+int64(nobs/25.)) + 1
    prob = np.zeros(size)
    ratio = np.zeros(size)
    for n in range(0, size):

        prob[n] = poisson(int64(n), float64(mu + nbkg))
        mu_best = max(0, n - nbkg)
        prob_best = poisson(int64(n), float64(mu_best + nbkg))

        if prob_best > 0:
            ratio[n] = float64(-1. * prob[n] / prob_best)
        else:
            ratio[n] = 0

    index_sorted = np.argsort(ratio)
    index_min = index_sorted[0]
    index_max = index_sorted[0]

    cumulative = np.float64(0.0)

    for i in range(len(index_sorted)):

        if index_sorted[i] < index_min:
            index_min = index_sorted[i]

        if index_sorted[i] > index_max:
            index_max = index_sorted[i]
        cumulative += prob[index_sorted[i]]

        if cumulative >= cl:
            break

    return True if (nobs <= index_max) and (nobs >= index_min) else False


@njit
def fc_limits(nobs, nbkg, cl):
    """Feldman-Cousins boundaries evaluation"""
    mu_lim = []
    switcher = False
    # small hack to speed up the evaluation
    # start close to the boundary
    # mu = 0
    mu = nobs - nobs**cl
    if nobs > 50:
        return nobs - np.sqrt(nobs), nobs + np.sqrt(nobs)
    while True:
        mu += 0.005
        good_choice = check_belt(mu, nobs, nbkg, cl)

        if good_choice != switcher:
            switcher = not switcher
            mu_lim.append(mu)
            # hack to speed up, move close to the next boundary
            mu += nobs
            if len(mu_lim) == 2:
                break

    return mu_lim[0], mu_lim[1]


def sigma_to_cl(sigma: float) -> float:
    """From significance in sigma to confidence level"""
    return 2*norm.cdf(sigma)-1


def cl_to_sigma(cl: float) -> float:
    """From confidence level to significance in sigma"""
    return norm.ppf(0.5*(1+cl))


@njit(fastmath=True)
def frequentist_uncertainty(nobs, cl):
    """Fequentist recipe for uncertainty evaluation (à la RooPlot)"""
    cond1 = False
    cond2 = False
    l = 0.0
    llo = 0.0
    lhi = 0.0
    if nobs == 0:
        cond1 = True
    while not (cond1 and cond2):
        l += 1e-5
        prob1 = poisson_cdf(nobs - 1, l)

        prob2 = prob1 + poisson(nobs, l)
        if prob1 > (1 + cl)/2:
            cond1 = True
            llo = l
        if prob2 < (1 - cl)/2:
            cond2 = True
            lhi = l

    return llo - nobs, lhi - nobs


def efficiency(events: int, passed: int, confidence_level: float) -> Tuple[float64, float64, float64]:
    """calculate efficiency and uncertainties from signal without background"""
    if (passed == 0): 
        OC.get_error_text(f"Empty bin! returning -1{OC.PM}-1!")
        OC.get_error_text("This should not be happening, please conact the authors.")
        return -1.0, -1.0, -1.0
    test = binomtest(passed, events)
    eff = test.proportion_estimate
    lower, higher = test.proportion_ci(confidence_level = confidence_level, method = 'wilsoncc')

    return eff, higher - eff, eff - lower


def calc_eff_and_error_no_bkg(n_tot, n_match, sigma: float = 1.0, debug: bool = True):
    """Evaluate efficiency and uncertainty using total number of events and matched events
    Args:
        n_tot: total number of events
        n_match: number of matched events
        sigma: statistical significance in n of sigma
        debug: print out efficiency and uncertainties
    """
    cl = sigma_to_cl(sigma)
    eff, hi_err, lo_err = efficiency(int(n_tot), int(n_match), cl)
    if debug: OC.get_info_text_blue("Total efficiency = %4.3f + %4.3f - %4.3f" % (eff, hi_err, lo_err))
    return eff, hi_err, lo_err


def calc_eff_and_error(n_tot, n_match, n_bkg, n_match_bkg, sigma = 1.0, debug = True):
    #number of sigmas calculated from CL
    cl = sigma_to_cl(sigma)

    #Physical bound
    n_match_bkg = min(n_bkg, n_match_bkg)

    #background efficiency
    eff_bkg = n_match_bkg/n_bkg

    #calculate signal yield and match
    n_sig = n_tot - n_bkg
    n_sig = max(n_sig, 1.)

    n_match_sig = n_match - n_match_bkg
    n_match_sig = max(0., n_match_sig)

    #Physical bound
    n_match_sig = min(n_sig, n_match_sig)

    #calculate efficiency
    eff_matched_sig =  np.divide(n_match_sig,n_sig, where = n_sig>0)

    ###calculate binomial uncertainties###
    err_matched_sig_up   = sigma*sqrt( eff_matched_sig*(1.0-eff_matched_sig) / n_sig )
    err_matched_sig_down = err_matched_sig_up

    #use FC for low statistics and close to boundaries
    if n_sig < 2000 and (5.*err_matched_sig_up+eff_matched_sig >= 1. or -5.*err_matched_sig_down+eff_matched_sig <= 0.):
        eff_matched_sig, err_matched_sig_up, err_matched_sig_down = efficiency(int(n_sig), int(n_match_sig), cl)

    if debug:
        OC.get_info_text_blue("Binomial efficiency error: +" + str(err_matched_sig_up) + " -" + str(err_matched_sig_down))

    #background uncertainty
    var_n_bkg_down = sigma*sqrt(n_bkg)
    var_n_bkg_up   = var_n_bkg_down

    #low background statistics
    if n_bkg < 25:
        lo, hi = fc_limits(n_bkg, 0, sigma_to_cl(1))
        var_n_bkg_down = n_bkg - lo
        var_n_bkg_up   = hi - n_bkg

    if var_n_bkg_down >= n_sig: var_n_bkg_down = n_sig - 1.

    #calculate effect of background fluctuations
    eff_bkg_up   = (n_match_sig + eff_bkg*var_n_bkg_up) / (n_sig + var_n_bkg_up)
    eff_bkg_down = (n_match_sig - eff_bkg*var_n_bkg_down) / (n_sig - var_n_bkg_down)

    eff_bkg_down = max(0, eff_bkg_down)
    err_n_bkg_up   = eff_matched_sig - eff_bkg_up if eff_bkg_up < eff_bkg_down else eff_matched_sig - eff_bkg_down
    err_n_bkg_down = -eff_matched_sig - eff_bkg_down if eff_bkg_up > eff_bkg_down else -eff_matched_sig + eff_bkg_up

    if debug:
        OC.get_info_text_blue("Poissonian background error: +" + str(err_n_bkg_up) + " -" + str(err_n_bkg_down))


    ###Binomial error for matched bkg###
    err_matched_bkg_up   = 0. if n_bkg == 0. else sigma*sqrt(1.0-eff_bkg)/n_bkg
    err_matched_bkg_down = err_matched_bkg_up

    #use FC for low statistics and close to boundaries
    if n_bkg < 2000 and (3./sigma*err_matched_bkg_up+eff_bkg >= 1. or -3./sigma*err_matched_bkg_down+eff_matched_sig <= 0.):
        #binomial error on matched bkg from FC
        eff_matched_bkg, err_matched_bkg_up, err_matched_bkg_down = efficiency(int(n_bkg), int(n_match_bkg), cl)

    #convert to frame of signal eff.
    err_matched_bkg_down = err_matched_bkg_down*n_bkg/n_sig
    err_matched_bkg_up   = err_matched_bkg_up*n_bkg/n_sig

    if debug:
        OC.get_info_text_blue("Binomial background error: +" + str(err_matched_bkg_up) + " -" + str(err_matched_bkg_down))


    ###Add up all errors quadratically
    err_tot_down = sqrt(err_matched_sig_down**2 + err_matched_bkg_down**2 + err_n_bkg_down**2)
    err_tot_up   = sqrt(err_matched_sig_up**2 + err_matched_bkg_up**2 + err_n_bkg_up**2)

    err_tot_down = min(err_tot_down, eff_matched_sig)
    err_tot_up = min(1 - eff_matched_sig, err_tot_up)

    if debug:
        OC.get_info_text_blue("Total efficiency =" + str(eff_matched_sig) + " +" + str(err_tot_up) + " -" + str(err_tot_down))

    return eff_matched_sig, err_tot_up, err_tot_down

#=============================================================
#
#           Final Method
#
#=============================================================

def getWeightFromErr(err): #See Renata's dissertation p.155
    if (err == 0):
        OC.get_warning_text("Your uncertainty is zero. Returning weight as zero.")
        return 0.0
    else: return 1.0/(err**2) 

def getFinalMethodErr(err_L, err_C):  #See Renata's dissertation p.155

    if (err_L==0 and err_C==0): 
        OC.get_warning_text("Both uncertanities for Long and Combined methods are zero. Returning zero.")
        return 0.0
    if (np.sign(err_L) == -np.sign(err_C)):
        OC.get_debug_text(f"err_L: {err_L}, err_C: {err_C}")
        #raise Exception(OC.get_error_text("Somehow your upper and lower uncertanties got mixed up. Abort."))
            
    w_L = getWeightFromErr(err_L)
    w_C = getWeightFromErr(err_C)
        
    return np.sign(err_L)*(1.0/np.sqrt(w_L+w_C))

def getFinalMethodVal(eff_L, err_L_D, err_L_U, eff_C, err_C_D, err_C_U): #See Renata's dissertation p.155
    w_L = getWeightFromErr(0.5*(-err_L_D+err_L_U))
    w_C = getWeightFromErr(0.5*(-err_C_D+err_C_U))
    if (w_L+w_C ==0): return 0.0    
    else: return (eff_L*w_L + eff_C*w_C)/(w_L+w_C)


def getFinalMethod(eff_list_Long, eff_err_list_Long, eff_list_Comb, eff_err_list_Comb):        
    eff_list_Final = list(map(getFinalMethodVal,eff_list_Long,eff_err_list_Long[0],eff_err_list_Long[1],eff_list_Comb,eff_err_list_Comb[0],eff_err_list_Comb[1]))
        
    eff_err_list_Final = []
    eff_err_list_Final.append(list(map(getFinalMethodErr,eff_err_list_Long[0],eff_err_list_Comb[0])
                                )
                            )
    eff_err_list_Final.append(list(map(getFinalMethodErr,eff_err_list_Long[1],eff_err_list_Comb[1])
                                )
                            )
    #Protect from errors > 1
    eff_list_Final = np.array(eff_list_Final)
    eff_err_list_Final = np.array(eff_err_list_Final)
    mask = eff_list_Final + eff_err_list_Final[1] > 1.0
    eff_err_list_Final[1,np.where(mask)] = 1-eff_list_Final[np.where(mask)]

    return eff_list_Final, eff_err_list_Final 

#=============================================================
#
#           Combined Method
#
#=============================================================

#Combined method
def getCombEffErr(val_a,err_a,val_b,err_b):

    if (err_a==0 and err_b==0): 
        OC.get_warning_text("Both uncertanities for T and Velo methods are zero. Returning 0.")
        return 0.0
    if (np.sign(err_a) == -np.sign(err_b)):
        OC.get_debug_text(f"err_a: {err_a}, err_b: {err_b}")
        raise Exception(OC.get_error_text("Somehow your upper and lower uncertanties got mixed up. Abort."))
        
    return np.sign(err_a)*np.sqrt( (val_b*err_a)**2 + (val_a*err_b)**2 )

def getCombEffVal(val_a,val_b):
    return val_a*val_b

def getCombinedMethod(eff_list_Velo, eff_err_list_Velo, eff_list_T, eff_err_list_T):

    eff_list_Combined = list(map(getCombEffVal,eff_list_Velo,eff_list_T))

    eff_err_list_Combined = [] 
    eff_err_list_Combined.append(list(map(getCombEffErr,eff_list_Velo, eff_err_list_Velo[0], 
                                        eff_list_T, eff_err_list_T[0])
                                    )
                                )
    eff_err_list_Combined.append(list(map(getCombEffErr,eff_list_Velo, eff_err_list_Velo[1], 
                                        eff_list_T, eff_err_list_T[1])
                                    )
                                )

    return eff_list_Combined, eff_err_list_Combined


#=============================================================
#
#           Ratio
#
#=============================================================

def getRatioVal(Data,MC): # Can be anything, but for simplicity we use data and MC as identifiers
                    # as that is the one we are most likely to use
    '''
        Takes the list efficiencies in 'data' and 'mc',
        returns the list of the 'data'/'mc' ratio
    '''
    if (MC == 0):
        OC.get_warning_text("The denominator (typically MC) is zero! Returning 0")
        return 0.0
    else: return Data/MC                       
    

def getRatioErr(Data,Data_err,MC, MC_err):    
    '''
        Takes the list of errors for 'data' and 'mc' points,
        returns the list of the 'data'/'mc' errors
    '''
    if (MC == 0):
        OC.get_warning_text("The denominator (typically MC) is zero! Returning 0")
        return 0
    return np.sqrt( (Data_err/MC)**2 + (Data*MC_err/(MC**2))**2 )
        
    
def getRatio(Data,Data_err,MC, MC_err):    

    ratio_val = list(map(getRatioVal,Data,MC))

    ratio_err = []
    ratio_err.append(list(map(getRatioErr,Data,Data_err[0],MC, MC_err[0])))
    ratio_err.append(list(map(getRatioErr,Data,Data_err[1],MC, MC_err[1])))        

    return ratio_val, ratio_err

#=============================================================
#
#           Realtive difference
#
#=============================================================


def getRelativeDiffVal(first, second):
    if ((first+second)==0): return 0.0
    return (second-first)/(first+second)

def getRelativeDiffErr(first_val, first_err, second_val, second_err):
    if ((first_val+second_val)==0): return 0.0
    d_first  = (2.0*second_val)/np.power(first_val+second_val,2)
    d_second = (2.0*first_val)/np.power(first_val+second_val,2)
    return np.sqrt(np.power(d_first*first_err,2)+np.power(d_second*second_err,2))

def getRelativeDiff(first_val, first_err, second_val, second_err):
    diff_val = list(map(getRelativeDiffVal,first_val, second_val))
    diff_err = []
    diff_err.append(list(map(getRelativeDiffErr, first_val, first_err[0], second_val, second_err[0])))
    diff_err.append(list(map(getRelativeDiffErr, first_val, first_err[1], second_val, second_err[1])))
    return diff_val, diff_err

def get_total_uncertainty(err_list): 
    '''
        Returns the total uncertainty: sqrt(E1*E1+...+EN*EN)
    '''
    return np.sqrt(np.sum(np.power(err_list,2)))
import numpy as np
from numba import njit, vectorize, float64, int64
from math import exp, gamma
from scipy.stats import norm
from scipy import stats as st
from typing import Tuple, Union, List
from ..TCFit.model import Model
from ..TCFit.basic_pdfs import PDF

__values = np.linspace(0, 100, 10000)
sigma = st.norm.cdf(-1)


def eval_unc(nobs: int) -> Tuple[float, float]:
    if nobs < 51:
        prob = st.poisson.sf(nobs - 1, __values)
        lower = nobs - __values[np.argmin(np.abs(prob - sigma))]

        prob = st.poisson.cdf(nobs, __values)
        upper = __values[np.argmin(np.abs(prob - sigma))] - nobs

    else:
        lower = nobs ** 0.5
        upper = lower
    return lower, upper


@vectorize([float64(int64, float64)], fastmath=True)
def poisson(k, mu):
    """Poisson probability"""
    if mu < 0 or k < 0:
        return 0
    return mu ** k * exp(-mu) / gamma(k + 1)


@vectorize([float64(int64, float64)], fastmath=True)
def poisson_cdf(k, mu):
    "Poisson comulative"
    prob = 0.0
    if k < 0:
        return 0
    for x in np.arange(0, k + 1):
        prob += poisson(x, mu)
    return prob


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
        if prob1 > (1 + cl) / 2:
            cond1 = True
            llo = l
        if prob2 < (1 - cl) / 2:
            cond2 = True
            lhi = l

    return [nobs - llo, lhi - nobs]


def sigma_to_cl(sigma: float) -> float:
    """From significance in sigma to confidence level"""
    return 2 * norm.cdf(sigma) - 1


def cl_to_sigma(cl: float) -> float:
    """From confidence level to significance in sigma"""
    return norm.ppf(0.5 * (1 + cl))


@njit(fastmath=True, parallel=True)
def check_belt(mu: float64, nobs: int64, nbkg: float64, cl: float64) -> bool:
    """Feldman-Cousin Belt construction with Poissonian probability"""
    size = 50 * (1 + int64(nobs / 25.)) + 1
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
    mu = nobs - nobs ** cl
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


def chi_square(pdf: Union[PDF, Model], boundaries: List[float], values: List[float], errors: List[List[float]] = None):
    """Chi2/NDOF evaluation
    Parameters
    ----------
        pdf : PDF or Model
        boundaries : binning boundaries of the observable, List of floats
        values : List of bin contents
        errors : List of lower and upper uncertainties


    Return
    ------
        chi_2 evaluation
    """
    __chi_square = 0
    __expected = pdf.integrate(boundaries)
    __errors = errors

    if errors is None:
        __errors = [eval_unc(value) for value in values]
    for index, value in enumerate(values):
        if value > __expected[index]:
            pull = (value - __expected[index]) / __errors[index][0]
        else:
            pull = (value - __expected[index]) / __errors[index][1]

        __chi_square += pull * pull

    __n_params = len(pdf.parameters)
    return __chi_square / (len(values) - __n_params)

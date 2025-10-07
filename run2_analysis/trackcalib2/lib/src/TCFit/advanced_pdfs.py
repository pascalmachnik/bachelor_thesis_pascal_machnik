import numba
import numpy as np
import math

from scipy.special import wofz
import scipy.integrate

from ..TCFit.pdf import PDF, RANGE_TYPE
from ..TCFit.error_messages import ErrorMessages
from ..TCFit.formula import Formula
from ..TCFit.parameter import Parameter
from ..TCFit.dataset import check_data_compatibility

import warnings

from typing import Union, Tuple, List

##########################################################################################################
#                                      GLOBAL CONSTANTS                                                  #
##########################################################################################################

PI = np.pi
SQRT2 = np.sqrt(2)

#########################################################################################################
#                           DOUBLESIDEDCRYSTALBALL IMPLEMENTATION                                       #
#########################################################################################################

@numba.njit(fastmath = {'reassoc', 'arcp', 'contract', 'afn' })
def double_sided_crystalball_integral(mu: np.float64, 
                                      sigma_L: np.float64, 
                                      alpha_L: np.float64, 
                                      n_L: np.float64, 
                                      sigma_R: np.float64, 
                                      alpha_R: np.float64, 
                                      n_R: np.float64, 
                                      limits: Union[list, tuple, np.array]) -> np.float64:

    use_log_L = False
    if np.abs(n_L - 1.0) < 1e-5:
        use_log_L = True

    use_log_R = False        
    if np.abs(n_R - 1.0) < 1e-5:
        use_log_R = True

    abs_sigma_L = np.abs(sigma_L)
    abs_alpha_L = np.abs(alpha_L)
    abs_sigma_R = np.abs(sigma_R)
    abs_alpha_R = np.abs(alpha_R)

    value = np.float64(0.0)

    z_min_L = (limits[0] - mu) / abs_sigma_L
    z_max_L = (limits[1] - mu) / abs_sigma_L
    z_min_R = (limits[0] - mu) / abs_sigma_R
    z_max_R = (limits[1] - mu) / abs_sigma_R

    if z_min_L < - abs_alpha_L:
        if z_max_L < - abs_alpha_L:
            a_L = np.power((n_L / abs_alpha_L), n_L) * np.exp(-0.5 * abs_alpha_L * abs_alpha_L)
            b_L = n_L / abs_alpha_L - abs_alpha_L
            
            if use_log_L:
                value += a_L * abs_sigma_L * (np.log(b_L - z_min_L) - np.log(b_L - z_max_L))
            else:
                value += a_L * abs_sigma_L / (1.0 - n_L) * (
                        1.0 / np.power(b_L - z_min_L, n_L - 1.0) - 1.0 / np.power(b_L - z_max_L, n_L - 1.0))
                    
        elif z_max_L <= 0:
            a_L = np.power((n_L / abs_alpha_L), n_L) * np.exp(-0.5 * abs_alpha_L * abs_alpha_L)
            b_L = n_L / abs_alpha_L - abs_alpha_L

            first_term = np.float64(0.0)
            
            if use_log_L:
                first_term += a_L * abs_sigma_L * (np.log(b_L - z_min_L) - np.log(b_L + abs_alpha_L))
            else:
                first_term += a_L * abs_sigma_L / (1.0 - n_L) * (
                            1.0 / np.power(b_L - z_min_L, n_L - 1.0) - 1.0 / np.power(b_L + abs_alpha_L, n_L - 1.0))

            second_term = abs_sigma_L * np.sqrt(PI / 2) * (
                        math.erf(z_max_L / SQRT2) - math.erf(-abs_alpha_L / SQRT2))

            value += first_term + second_term
        
        elif z_max_R <= abs_alpha_R:
            a_L = np.power((n_L / abs_alpha_L), n_L) * np.exp(-0.5 * abs_alpha_L * abs_alpha_L)
            b_L = n_L / abs_alpha_L - abs_alpha_L

            first_term = np.float64(0.0)
            
            if use_log_L:
                first_term += a_L * abs_sigma_L * (np.log(b_L - z_min_L) - np.log(b_L + abs_alpha_L))
            else:
                first_term += a_L * abs_sigma_L / (1.0 - n_L) * (
                            1.0 / np.power(b_L - z_min_L, n_L - 1.0) - 1.0 / np.power(b_L + abs_alpha_L, n_L - 1.0))

            second_term = abs_sigma_L * np.sqrt(PI / 2) * (
                        0 - math.erf(-abs_alpha_L / SQRT2)) # math.erf(0 / SQRT2) = 0

            third_term = abs_sigma_R * np.sqrt(PI / 2) * (
                        math.erf(z_max_R / SQRT2) - 0) #  - math.erf(0 / SQRT2) = 0
            
            value += first_term + second_term + third_term
            
        else:
            a_L = np.power((n_L / abs_alpha_L), n_L) * np.exp(-0.5 * abs_alpha_L * abs_alpha_L)
            b_L = n_L / abs_alpha_L - abs_alpha_L
            
            a_R = np.power((n_R / abs_alpha_R), n_R) * np.exp(-0.5 * abs_alpha_R * abs_alpha_R)
            b_R = n_R / abs_alpha_R - abs_alpha_R

            first_term = np.float64(0.0)
            
            if use_log_L:
                first_term += a_L * abs_sigma_L * (np.log(b_L - z_min_L) - np.log(b_L + abs_alpha_L))
            else:
                first_term += a_L * abs_sigma_L / (1.0 - n_L) * (
                            1.0 / np.power(b_L - z_min_L, n_L - 1.0) - 1.0 / np.power(b_L + abs_alpha_L, n_L - 1.0))

            second_term = abs_sigma_L * np.sqrt(PI / 2) * (
                        0 - math.erf(-abs_alpha_L / SQRT2)) # math.erf(0 / SQRT2) = 0

            third_term = abs_sigma_R * np.sqrt(PI / 2) * (
                        math.erf(abs_alpha_R / SQRT2) - 0) #  - math.erf(0 / SQRT2) = 0

            if use_log_R:
                fourth_term = a_R * abs_sigma_R * (np.log(b_R + z_max_R) - np.log(b_R + abs_alpha_R))
            else:
                fourth_term = a_R * abs_sigma_R / (1.0 - n_R) * (
                            1.0 / np.power(b_R + z_max_R, n_R - 1.0) - 1.0 / np.power(b_R + abs_alpha_R, n_R - 1.0))
            
            value += first_term + second_term + third_term + fourth_term
            
    elif z_min_L <= 0:
        if z_max_L <= 0:
            value += abs_sigma_L * np.sqrt(PI / 2) * (
                        math.erf(z_max_L / SQRT2) - math.erf(z_min_L / SQRT2))
        
        elif z_max_R <= abs_alpha_R:

            first_term = np.float64(0.0)

            first_term = abs_sigma_L * np.sqrt(PI / 2) * (
                        0 - math.erf(z_min_L / SQRT2)) # math.erf(0 / SQRT2) = 0

            second_term = abs_sigma_R * np.sqrt(PI / 2) * (
                        math.erf(z_max_R / SQRT2) - 0) #  - math.erf(0 / SQRT2) = 0
            
            value += first_term + second_term
            
        else:
            a_R = np.power((n_R / abs_alpha_R), n_R) * np.exp(-0.5 * abs_alpha_R * abs_alpha_R)
            b_R = n_R / abs_alpha_R - abs_alpha_R

            first_term = np.float64(0.0)
            
            first_term = abs_sigma_L * np.sqrt(PI / 2) * (
                        0 - math.erf(z_min_L / SQRT2)) # math.erf(0 / SQRT2) = 0

            second_term = abs_sigma_R * np.sqrt(PI / 2) * (
                        math.erf(abs_alpha_R / SQRT2) - 0) #  - math.erf(0 / SQRT2) = 0

            if use_log_R:
                third_term = a_R * abs_sigma_R * (np.log(b_R + z_max_R) - np.log(b_R + abs_alpha_R))
            else:
                third_term = a_R * abs_sigma_R / (1.0 - n_R) * (
                            1.0 / np.power(b_R + z_max_R, n_R - 1.0) - 1.0 / np.power(b_R + abs_alpha_R, n_R - 1.0))
            
            value += first_term + second_term + third_term
            
    elif z_min_R <= abs_alpha_R:
        if z_max_R <= abs_alpha_R:
            value += abs_sigma_R * np.sqrt(PI / 2) * (
                        math.erf(z_max_R / SQRT2) - math.erf(z_min_L / SQRT2))
            
        else:
            a_R = np.power((n_R / abs_alpha_R), n_R) * np.exp(-0.5 * abs_alpha_R * abs_alpha_R)
            b_R = n_R / abs_alpha_R - abs_alpha_R

            first_term = np.float64(0.0)

            first_term = abs_sigma_R * np.sqrt(PI / 2) * (
                        math.erf(abs_alpha_R / SQRT2) - math.erf(z_min_L / SQRT2))

            if use_log_R:
                second_term = a_R * abs_sigma_R * (np.log(b_R + z_max_R) - np.log(b_R + abs_alpha_R))
            else:
                second_term = a_R * abs_sigma_R / (1.0 - n_R) * (
                            1.0 / np.power(b_R + z_max_R, n_R - 1.0) - 1.0 / np.power(b_R + abs_alpha_R, n_R - 1.0))
            
            value += first_term + second_term
            
    else:
        a_R = np.power((n_R / abs_alpha_R), n_R) * np.exp(-0.5 * abs_alpha_R * abs_alpha_R)
        b_R = n_R / abs_alpha_R - abs_alpha_R
        
        if use_log_R:
            value += a_R * abs_sigma_R * (np.log(b_R + z_max_R) - np.log(b_R + z_min_R))
        else:
            value += a_R * abs_sigma_R / (1.0 - n_R) * (
                            1.0 / np.power(b_R + z_max_R, n_R - 1.0) - 1.0 / np.power(b_R + z_min_R, n_R - 1.0))
            
    if value > 1e-323: return value
    else:              return np.float64(1e-38)
    
#@numba.njit(fastmath = {'reassoc', 'arcp', 'contract', 'afn' }, parallel = True)
def double_sided_crystalball_integral_vec(mu: np.float64, 
                                          sigma_L: np.float64, 
                                          alpha_L: np.float64, 
                                          n_L: np.float64, 
                                          sigma_R: np.float64, 
                                          alpha_R: np.float64, 
                                          n_R: np.float64, 
                                          ranges: Union[list, tuple, np.array]) -> np.float64:

    size = ranges.shape[0]
    if size < 2:
        raise ValueError("Size of limits < 2")

    integrals = np.zeros(shape = size-1, dtype = np.float64)
            
    for i in numba.prange(size - 1):
        integrals[i] = double_sided_crystalball_integral(mu, sigma_L, alpha_L, n_L, sigma_R, alpha_R, n_R, [ranges[i], ranges[i+1]])
    return integrals

def fxn():
    warnings.warn("runtime", RuntimeWarning)

#@numba.njit(fastmath = {'reassoc', 'arcp', 'contract', 'afn' }) 
def double_sided_crystalball_unnormalized(x: Union[np.float64, np.ndarray],
                             mu: Union[np.float64, np.ndarray],
                             sigma_l: np.float64, 
                             alpha_l: np.float64, 
                             n_l: np.float64, 
                             sigma_r: np.float64, 
                             alpha_r: np.float64, 
                             n_r: np.float64) -> Union[np.float64, np.ndarray]:
                
    """Implementation of the double-sided crystalball function accelerated by numba

    Parameters
    ----------
        x       : int, float, np.ndarray
            `x` is evaluated on the double-sided crystalball function
        mu      : int, float
            Determines the central location of the double-sided crystalball function
        sigma_l   : int, float
            Determines the left shape of the double-sided crystalball function
        alpha_l : int, float
            Controls the left core tail
        n_l     : int
            Controls the shape of the left tail
        sigma_r   : int, float
            Determines the right shape of the double-sided crystalball function
        alpha_r : int, float
            Controls the right core tail
        n_r     : int
            Controls the shape of the right tail

    Returns
    -------
        int, float, np.ndarray
            The return value has the same type and shape as `x`
    """
    
    #Sanity check for NaNs
    if (np.isnan(mu) or np.isnan(sigma_l) or np.isnan(alpha_l) or np.isnan(n_l)
                     or np.isnan(sigma_r) or np.isnan(alpha_r) or np.isnan(n_r) ):
        print("Some value is NaN")
        return x - x + 1e-323
    #Check the parameter values
    if ((np.abs(alpha_l) <= 1e-12) or (np.abs(alpha_r) <= 1e-12) or 
        (sigma_l == 0.0) or (sigma_r == 0.0)):
        return x - x + 1e-323 #Returns an array with 1e-323 everywhere
    
    z_l = (x - mu) / (sigma_l)
    a_l = np.power((n_l / np.abs(alpha_l)), n_l) * np.exp(-0.5 * alpha_l * alpha_l)
    b_l = (n_l / np.abs(alpha_l)) - np.abs(alpha_l)
    z_r = (x - mu) / (sigma_r)
    a_r = np.power((n_r / np.abs(alpha_r)), n_r) * np.exp(-0.5 * alpha_r * alpha_r)
    b_r = (n_r / np.abs(alpha_r)) - np.abs(alpha_r)

    cond_l = np.less(z_l, -np.abs(alpha_l))
    cond_cl = np.greater_equal(z_l, -np.abs(alpha_l)) & np.less_equal(z_l, 0)
    cond_cr = np.greater(z_r, 0) & np.less_equal(z_r, np.abs(alpha_r))
    cond_r = np.greater(z_r, np.abs(alpha_r))
    
    value = np.select([cond_l, 
                       cond_cl,
                       cond_cr, 
                       cond_r],
                      [a_l * np.power(b_l - z_l, -n_l), 
                       np.exp(-0.5 * np.power(z_l, 2)), 
                       np.exp(-0.5 * np.power(z_r, 2)),
                       a_r * np.power(b_r + z_r, -n_r)])
                
    value = np.where(np.isnan(value), np.full(value.shape, 1e-323), value)
    
    return value

# The fastmath optins are needed in order to be able to check nan and inf.
# If just set to fastmath=True, numba does not return inf/nan as infs and nans,
# leading to problems when doing sanity checks on the parameters.
# The parameters can be nan as iminuit can return nan.

#@numba.njit(fastmath = {'reassoc', 'arcp', 'contract', 'afn' }) 
def double_sided_crystalball(x: Union[np.float64, np.ndarray],
                             mu: Union[np.float64, np.ndarray],
                             sigma_l: np.float64, 
                             alpha_l: np.float64, 
                             n_l: np.float64, 
                             sigma_r: np.float64, 
                             alpha_r: np.float64, 
                             n_r: np.float64, 
                             limits: Union[list, tuple, np.array]) -> Union[np.float64, np.ndarray]:
    """Implementation of the double-sided crystalball function accelerated by numba

    Parameters
    ----------
        x       : int, float, np.ndarray
            `x` is evaluated on the double-sided crystalball function
        mu      : int, float
            Determines the central location of the double-sided crystalball function
        sigma_l   : int, float
            Determines the left shape of the double-sided crystalball function
        alpha_l : int, float
            Controls the left core tail
        n_l     : int
            Controls the shape of the left tail
        sigma_r   : int, float
            Determines the right shape of the double-sided crystalball function
        alpha_r : int, float
            Controls the right core tail
        n_r     : int
            Controls the shape of the right tail

    Returns
    -------
        int, float, np.ndarray
            The return value has the same type and shape as `x`
    """
    normalisation  = double_sided_crystalball_integral(mu, sigma_l, alpha_l, n_l, sigma_r, alpha_r, n_r, limits)
    DSCB_value = double_sided_crystalball_unnormalized(x, mu, sigma_l, alpha_l, n_l, sigma_r, alpha_r, n_r)
    
    if normalisation > 1e-323:
        value = DSCB_value / normalisation
    else:
        value = DSCB_value - DSCB_value + 1e323

    return value


class DoubleSidedCrystalBall(PDF):
    """Double-sided crystalball probability density function

    Parameters
    ----------
        name       : str
            The name of the double-sided crystalball pdf
        observable : parameter.Parameter
            The physical observable
        mu         : parameter.Parameter, utils.formula.Formula
            Positional shift of the double-sided crystalball
        sigma      : parameter.Parameter, utils.formula.Formula
            Shape of the double-sided crystalball
        alpha_l    : parameter.Parameter, utils.formula.Formula
            Controls the left core tail
        n_l        : parameter.Parameter, utils.formula.Formula
            Controls the shape of the left tail
        alpha_r    : parameter.Parameter, utils.formula.Formula
            Controls the right core tail
        n_r        : parameter.Parameter, utils.formula.Formula
            Controls the shape of the right tail

    Attributes
    ----------
        parameters : list[parameter.Parameter]
            parameters contains `mu`, `sigma`, `alpha_l`, `n_l`, `alpha_r`, `n_r`, in this order
    """

    def __init__(self, name: str, observable: Parameter, mu: Union[Parameter, Formula],
                 sigma_l: Union[Parameter, Formula], alpha_l: Union[Parameter, Formula], n_l: Union[Parameter, Formula],
                 sigma_r: Union[Parameter, Formula], alpha_r: Union[Parameter, Formula], n_r: Union[Parameter, Formula]):
        super().__init__(name, observable)
        self.mu = mu
        self.sigma_l = sigma_l
        self.alpha_l = alpha_l
        self.n_l = n_l
        self.sigma_r = sigma_r
        self.alpha_r = alpha_r
        self.n_r = n_r
        self.set_parameters([mu, sigma_l, alpha_l, n_l, sigma_r, alpha_r, n_r])

    def evaluate(self, x: Union[np.float64, np.ndarray]) -> Union[np.float64, np.ndarray]:
        """`double_sided_crystalball` is evaluated on `x` with the initialised parameters. Before evaluation, `x`
        is constrained to only those values which are inside of the observable limits
        
        Parameters
        ----------
            x : int, float, np.ndarray
                `x` has to be a numeric value or array-like so that the `double_sided_crystalball` can be evaluated

        Returns
        -------
            int, float, np.ndarray
                The return value has the same type as `x`
        """

        return double_sided_crystalball(x, self.mu.value, 
                                        self.sigma_l.value, self.alpha_l.value, self.n_l.value,
                                        self.sigma_r.value, self.alpha_r.value, self.n_r.value,
                                        self.observable.limits)

    def integrate(self, ranges: Union[list, tuple, np.ndarray]) -> Union[Tuple[float, float], Tuple[List[float], List[float]]]:

        normalisation = double_sided_crystalball_integral(self.mu.value,
                                                          self.sigma_l.value,
                                                          self.alpha_l.value,
                                                          self.n_l.value, 
                                                          self.sigma_r.value,
                                                          self.alpha_r.value,
                                                          self.n_r.value, 
                                                          self.observable.limits)

        return double_sided_crystalball_integral_vec(self.mu.value,
                                                     self.sigma_l.value,
                                                     self.alpha_l.value,
                                                     self.n_l.value, 
                                                     self.sigma_r.value,
                                                     self.alpha_r.value,
                                                     self.n_r.value, 
                                                     ranges)/normalisation
        
    @staticmethod
    def func(x: Union[np.float64, np.ndarray], mu: np.float64,
             sigma_l: np.float64, alpha_l: np.float64, n_l: np.float64,
             sigma_r: np.float64, alpha_r: np.float64, n_r: np.float64, 
             **kwargs) -> Union[np.float64, np.ndarray]:
        """Static Version of the `evaluate` method
        
        Parameters
        ----------
            x          : int, float, np.ndarray
                `x` has to be a numeric value or array-like so that the `double_sided_crystalball` can be evaluated
            mu         : int, float
                Positional shift of the double-sided crystall ball
            sigma      : int, float
                Shape of the double-sided crystal ball
            alpha_l    : int, float
                Controls the left core tail
            n_l        : int
                Controls the shape of the left tail
            alpha_r    : int, float
                Controls the right core tail
            n_r        : int
                Controls the shape of the right tail
            **kwargs:
                limits : array-like
                    `limits` has to be of size `2` and is used in order to correctly normalise the
                    double-sided crystalball

        Returns
        -------
            int, float, np.ndarray
                The return value has the same type and shape as `x`
        """
        if "limits" in kwargs:
            if not isinstance(kwargs["limits"], (list, tuple, np.ndarray)) or len(kwargs["limits"]) != 2:
                raise TypeError(ErrorMessages.get_error_message(4))
            x = check_data_compatibility(kwargs["limits"], x)
        else:
            kwargs["limits"] = (-1e10 * sigma_l - mu, +1e10 * sigma_r + mu)

        return double_sided_crystalball(x, mu, sigma_l, alpha_l, n_l, sigma_r, alpha_r, n_r, kwargs["limits"])


### Missing to be added: argus, johnson

#########################################################################################################
#                           ARGUS IMPLEMENTATION                                                        #
#########################################################################################################
@numba.jit(nopython=True, fastmath=True)
def argus_integral(m0: float, c: float, p: float, low: Union[int, float], high: Union[int, float]):
    """Integral of the argus function between low and high values, used to properly normalise the argus pdf.
    The analytical function is used in case p == 0.5.

    Parameters
    ----------
        name : str
            The name of the crystall ball function
        observable : parameter.Parameter
            Physical observable
        m0: Resonance mass
        c: Slope parameter
        p: Power
        low   : int, float
            `low` is used to determine the lower boundary
        high  : int, float
            `high` is used to determine the higher boundary

    Returns
    ----------
        float
            The computed area under the crystalball function is returned
    """
    if p == 0.5:
        f1 = (1. - np.pow(low / m0, 2))
        f2 = (1. - np.pow(high / m0, 2))
        if c < 0.:
            a_low  = -0.5 * m0 * m0 * (np.exp(c * f1) * np.sqrt(f1) / c + 0.5 / np.pow(-c, 1.5) * np.sqrt(PI) * math.erf(np.sqrt(-c * f1)))
            a_high = -0.5 * m0 * m0 * (np.exp(c * f2) * np.sqrt(f2) / c + 0.5 / np.pow(-c, 1.5) * np.sqrt(PI) * math.erf(np.sqrt(-c * f2)))
        elif c == 0.:
            a_low  = -m0 * m0 / 3. * f1 * np.sqrt(f1)
            a_high = -m0 * m0 / 3. * f1 * np.sqrt(f2)
        else:
            a_low  = 0.5 * m0 * m0 * np.exp(c * f1) / (c * np.sqrt(c)) * (0.5 * np.sqrt(PI) * (wofz(np.sqrt(c * f1))).imag - np.sqrt(c * f1))
            a_high = 0.5 * m0 * m0 * np.exp(c * f2) / (c * np.sqrt(c)) * (0.5 * np.sqrt(PI) * (wofz(np.sqrt(c * f2))).imag - np.sqrt(c * f2))

        func = a_high - a_low

    else:
        func, _ = scipy.integrate(argus, low, high)

    return func


@numba.jit(nopython=True, fastmath=True)
def argus(x: Union[int, float, np.ndarray], m0: float, c: float, p: float, low: Union[int, float], high: Union[int, float]):
    """Implementation of the argus function between low and high values, used to properly normalise the argus pdf.
    The analytical function is used in case p == 0.5.

    Parameters
    ----------
        name : str
            The name of the crystall ball function
        observable : parameter.Parameter
            Physical observable
        m0: Resonance mass
        c: Slope parameter
        p: Power (by default set to 0.5)
        low   : int, float
            `low` is used to determine the lower boundary
        high  : int, float
            `high` is used to determine the higher boundary

    Returns
    ----------
        float, np.ndarray
            The return value has the same type as `x`
    """
    t = x / m0;
    if t >= 1:
        return 0;

    u = 1 - t * t;
    return x * np.pow(u, p) * np.exp(c * u)/argus_integral(m0, c, p, low, high)


class Argus(PDF):
    """Argus probability density function

    Parameters
    ----------
        name : str
            The name of the crystall ball function
        observable : parameter.Parameter
            Physical observable
        m0: Resonance mass
        c: Slope parameter
        p: Power (by default set to 0.5)

    Attributes
    ----------
        parameters: list[parameter.Parameter or utils.formula.Formula]
            `parameters` is a list which contains `m0`, `c` and `p`. This list is used in order
            to keep tracks and simply access to the parameters of the Argus
        """
    def __init__(self, name: str, observable: Parameter,
                 m0: Union[Parameter,Formula],
                 c: Union[Parameter,Formula],
                 **kwargs):
        super().__init__(name, observable)
        self.m0 = m0
        self.c = c
        if "p" in kwargs:
            self.p = kwargs["p"]
        else:
            self.p = Parameter(f"{name}_p", 0.5, (0.5, 0.5)).fixed

        self.set_parameters([self.m0, self.c, self.p])


    def evaluate(self, x: Union[int, float, np.ndarray]) -> Union[int, float, np.ndarray]:
        """`argus` is evaluated on `x` with the initialised parameters. Before evaluation, `x`
        is constrained to only those values which are inside of the observable limit

        Parameters
        ----------
            x : int, float, np.ndarray
                `x` has to be a numeric value or array-like so that the `Argus` can be evaluated

        Returns
        -------
            int, float, np.ndarray
                The return value has the same type as `x`
        """

        return argus(x,
                     self.m0.value,
                     self.c.value,
                     self.p.value,
                     self.observable.limits[0],
                     self.observable.limits[1])

    @staticmethod
    def func(x: Union[int, float, np.ndarray], m0: float, c: float, p: float, **kwargs) -> Union[int, float, np.ndarray]:
        """Static Version of the `evaluate` method

        Parameters
        ----------
            x          : int, float, np.ndarray
                `x` has to be a numeric value or array-like so that the `crystalball` can be evaluated
            m0: float,
                Resonance mass
            c: float
                Slope parameter
            p: float
                Power
            **kwargs:
                limits : array-like
                    `limits` has to be of size `2` and is used in order to correctly normalise the
                    argus

        Returns
        -------
            float, np.ndarray
                The return value has the same type as `x`
        """
        if "limits" in kwargs:
            if not isinstance(kwargs["limits"], (list, tuple, np.ndarray)) or len(kwargs["limits"]) != 2:
                raise TypeError(ErrorMessages.get_error_message(4))
            x = check_data_compatibility(kwargs["limits"], x)
        else:
            kwargs["limits"] = (-1e10 - m0, +1e10 + m0)

        return argus(x, m0, c, p, kwargs["limits"][0], kwargs["limits"][1])



#########################################################################################################
#                                      JOHNSON IMPLEMENTATION                                           #
#########################################################################################################

def johnson_integral(
    mu: np.float64, 
    lambd: np.float64, 
    gamma: np.float64, 
    delta: np.float64, 
    limits: Union[list, tuple, np.array]
) -> np.float64:

    #value, _ = scipy.integrate.quad(johnson_unnormalized, limits[0], limits[1], args = (mu, lambd, gamma, delta))
    value, _ = scipy.integrate.quad(lambda x: johnson_unnormalized(x, mu, lambd, gamma, delta), limits[0], limits[1])
    
    if value > 1e-323: return value
    else:              return np.float64(1e-38)


#@numba.njit(fastmath = {'reassoc', 'arcp', 'contract', 'afn' }, parallel = True)
def johnson_integral_vec(
    mu: np.float64, 
    lambd: np.float64, 
    gamma: np.float64, 
    delta: np.float64, 
    ranges: Union[list, tuple, np.ndarray]) -> np.ndarray:

    size = ranges.shape[0]
    if size < 2:
        raise ValueError("Size of limits < 2")

    integrals = np.zeros(shape = size-1, dtype = np.float64)
        
    for i in numba.prange(size - 1):
        integrals[i] = johnson_integral(mu, lambd, gamma, delta, [ranges[i], ranges[i+1]])
    
    return integrals

# The fastmath optins are needed in order to be able to check nan and inf.
# If just set to fastmath=True, numba does not return inf/nan as infs and nans,
# leading to problems when doing sanity checks on the parameters.
# The parameters can be nan as iminuit can return nan.

@numba.njit(fastmath = {'reassoc', 'arcp', 'contract', 'afn' }) 
def johnson_unnormalized(
    x: Union[np.float64, np.ndarray],
    mu: np.float64, 
    lambd: np.float64, 
    gamma: np.float64, 
    delta: np.float64, 
    ) -> Union[np.float64, np.ndarray]:
                
    """Implementation of the Johnson function without normalization

    Parameters
    ----------
        x : float
            `x` is evaluated on the Johnson function
        mu    : float
            Determines the central location of the Johnson function
        lambd : float
            Determines the shape of the Johnson function
        gamma : float
            Controls the central shape
        delta     : float
            Controls the shape of the tails

    Returns
    -------
        int, float, np.ndarray
            The return value has the same type as `x`
    """

    #Sanity check for NaNs
    if (np.isnan(mu) or np.isnan(lambd) or np.isnan(gamma) or np.isnan(delta)):
        return x - x + 1e-323
    #Check the parameter values
    if (lambd == 0.0) or (mu < 0):
        return x - x + 1e-323 #Returns an array with 1e-323 everywhere

    z = (x - mu) / lambd
    value = np.exp(-0.5 * np.power(gamma + delta * np.arcsinh(z), 2))
    
    return value

#@numba.njit(fastmath = {'reassoc', 'arcp', 'contract', 'afn' }) 
def johnson(
    x: Union[np.float64, np.ndarray],
    mu: np.float64, 
    lambd: np.float64, 
    gamma: np.float64, 
    delta: np.float64, 
    limits: Union[list, tuple, np.array]
    ) -> Union[np.float64, np.ndarray]:
                
    """Implementation of the Johnson function

    Parameters
    ----------
        x : float
            `x` is evaluated on the Johnson function
        mu    : float
        lambd : float
        gamma : float
        delta : float
        limits : array-like                 
            `limits` has to be of size `2` and is used for normalisation of the Johnson function

    Returns
    -------
        int, float, np.ndarray
            The return value has the same type as `x`
    """

    normalisation  = johnson_integral(mu, lambd, gamma, delta, limits)
    johnson_value = johnson_unnormalized(x, mu, lambd, gamma, delta)
    
    if normalisation > 1e-323:
        value = johnson_value / normalisation
    else:
        value = johnson_value - johnson_value + 1e323

    return value


class Johnson(PDF):
    """Johnson probability density function

    Parameters
    ----------
        name : str
            The name of the crystall ball function
        observable : parameter.Parameter
            Physical observable
        mu         : parameter.Parameter, utils.formula.Formula
        lambd      : parameter.Parameter, utils.formula.Formula
        gamma      : parameter.Parameter, utils.formula.Formula
        delta      : parameter.Parameter, utils.formula.Formula
            
    Attributes
    ----------
        parameters : list[parameter.Parameter or utils.formula.Formula]
            `parameters` is a list which contains `mu`, `lambd`, `gamma` and `delta`. This list is used in order 
            to keep tracks and simply access to the parameters of the Johnson function
    """

    def __init__(
        self, name: str, 
        observable: Parameter,
        mu: Union[Parameter, Formula], 
        lambd: Union[Parameter, Formula],
        gamma: Union[Parameter, Formula], 
        delta: Union[Parameter, Formula]
        ):
        super().__init__(name, observable)
        self.mu = mu
        self.lambd = lambd
        self.gamma = gamma
        self.delta = delta
        self.set_parameters([self.mu, self.lambd, self.gamma, self.delta])

    def evaluate(self, x: Union[np.float64, np.ndarray]) -> Union[np.float64, np.ndarray]:
        """`Johnson` is evaluated on `x` with the initialised parameters. Before evaluation, `x`
        is constrained to only those values which are inside of the observable limit
        
        Parameters
        ----------
            x : int, float, np.ndarray
                `x` has to be a numeric value or array-like so that the `johnson` function can be evaluated

        Returns
        -------
            int, float, np.ndarray
                The return value has the same type as `x`
        """

        return johnson(
            x,
            self.mu.value,
            self.lambd.value,
            self.gamma.value,
            self.delta.value,
            self.observable.limits
        )

    def integrate(self, ranges: Union[list, tuple, np.ndarray]) -> Union[Tuple[float, float], Tuple[List[float], List[float]]]:

        normalisation = johnson_integral(
            self.mu.value,
            self.lambd.value,
            self.gamma.value,
            self.delta.value, 
            self.observable.limits
            )
        return johnson_integral_vec(
            self.mu.value,
            self.lambd.value,
            self.gamma.value,
            self.delta.value, 
            ranges
            )/normalisation

    @staticmethod
    def func(
        x: Union[np.float64, np.ndarray], 
        mu: np.float64, 
        lambd: np.float64, 
        gamma: np.float64, 
        delta: np.float64, 
        **kwargs
        ) -> Union[np.float64, np.ndarray]:
        """Static Version of the `evaluate` method
        
        Parameters
        ----------
            x          : int, float, np.ndarray
                `x` has to be a numeric value or array-like so that the `johnson` can be evaluated
            mu         : parameter.Parameter, utils.formula.Formula
            lambd      : parameter.Parameter, utils.formula.Formula
            gamma      : parameter.Parameter, utils.formula.Formula
            delta      : parameter.Parameter, utils.formula.Formula
            **kwargs:
                limits : array-like
                    `limits` has to be of size `2` and is used in order to correctly normalise the Johnson function

        Returns
        -------
            int, float, np.ndarray
                The return value has the same type as `x`
        """
        if "limits" in kwargs:
            if not isinstance(kwargs["limits"], (list, tuple, np.ndarray)) or len(kwargs["limits"]) != 2:
                raise TypeError(ErrorMessages.get_error_message(4))
            x = check_data_compatibility(kwargs["limits"], x)
        else:
            kwargs["limits"] = (-1e10 * lambd - mu, +1e10 * lambd + mu)


        return johnson(x, mu, lambd, gamma, delta, kwargs["limits"])
    
#########################################################################################################
#                           CHEBYSHEV IMPLEMENTATION                                                   #
#########################################################################################################

#@numba.njit(fastmath = {'reassoc', 'arcp', 'contract', 'afn' })
def chebyshev_single_polynomial(x: Union[np.float64, np.ndarray], degree: int) -> Union[np.float64, np.ndarray]:
    # first few Chebyshev polynomials of the first kind
    if degree == 1:
        return x
    elif degree == 2:
        return 2 * np.power(x,2) - 1
    elif degree == 3:
        return 4 * np.power(x,3) - 3 * x
    elif degree == 4:
        return 8 * np.power(x,4) - 8 * np.power(x,2) + 1
    elif degree == 5:
        return 16 * np.power(x,5) - 20 * np.power(x,3) + 5 * x
    
#@numba.njit(fastmath = {'reassoc', 'arcp', 'contract', 'afn' })
def chebyshev_integral(coefficients: np.ndarray, 
                        limit: Union[list, tuple, np.ndarray],
                        main_limit: Union[list, tuple, np.ndarray]) -> np.float64:
    # Integral is calculated iteratively following: https://en.wikipedia.org/wiki/Chebyshev_polynomials

    # as first coeffecient is implicitly 1, the number of coefficients is the degree - 1
    degree_m1 = coefficients.shape[0]

    halfrange = (main_limit[1] - main_limit[0]) / 2.
    mid = (main_limit[1] + main_limit[0]) / 2.
    
    # scale integral range to be between -1 and 1
    a = (limit[0] - mid) / halfrange
    b = (limit[1] - mid) / halfrange

    summed = b - a # zeroth order integral
    
    if degree_m1 > 0:
        c = coefficients[0]
        summed += (0.5 * (b + a) * (b - a)) * c # first order integral
        
        # recursively calculate higher order integrals
        if degree_m1 > 1:
            # calculate lower edge iteration
            a_curr = a
            a_prev = 1
            
            # T_(n+1)(x) = 2x*T_n(x) - T_(n-1)(x)
            new_val = 2 * a * a_curr - a_prev
            a_prev  = a_curr
            a_curr  = new_val
            
            # calculate upper edge iteration
            b_curr = b
            b_prev = 1
            
            # T_(n+1)(x) = 2x*T_n(x) - T_(n-1)(x)
            new_val = 2 * b * b_curr - b_prev
            b_prev = b_curr
            b_curr = new_val
            
            # n_m1 = 1 is second order integral
            # n_m1 = 2 is third order integral
            # etc.
            for n_m1 in range(1, degree_m1):
                # calculate T_(n-1) integral
                c = coefficients[n_m1]
                
                # term including T_(n-1)(x)
                term_nm1 = (b_prev - a_prev) / n_m1
                
                # T_(n+1)(x) = 2x*T_n(x) - T_(n-1)(x)
                new_val = 2 * a * a_curr - a_prev
                a_prev = a_curr
                a_curr = new_val
                
                # T_(n+1)(x) = 2x*T_n(x) - T_(n-1)(x)
                new_val = 2 * b * b_curr - b_prev
                b_prev = b_curr
                b_curr = new_val
                
                # term including T_(n+1)(x)
                term_np1 = (b_curr - a_curr) / (n_m1 + 2)
                
                # full term: integral of T_n(x) in [-1, 1] = 0.5 * (T_(n+1)/(n+1) - T_(n-1)/(n-1))
                int_Tn = 0.5 * (term_np1 - term_nm1)
                summed += int_Tn * c

    return halfrange*summed


#@numba.njit(fastmath = {'reassoc', 'arcp', 'contract', 'afn' }, parallel = True)
def chebyshev_integral_vec(coefficients: np.ndarray, 
                           limits: Union[list, tuple, np.ndarray],
                           main_limit: Union[list, tuple, np.ndarray]) -> np.ndarray:

    size = limits.shape[0]
    if size < 2:
        raise ValueError("Size of limits < 2")

    integrals = np.zeros(shape = size-1, dtype = np.float64)
        
    for i in numba.prange(size - 1):
        integrals[i] = chebyshev_integral(coefficients, [limits[i], limits[i+1]], main_limit)
    
    return integrals
    

#@numba.njit(fastmath = {'reassoc', 'arcp', 'contract', 'afn' })
def chebyshev(x: Union[np.float64, np.ndarray], 
               coefficients: np.ndarray, 
               limit: Union[list, tuple, np.ndarray]) -> Union[np.float64, np.ndarray]:
    # need to scale x to use Chebyshev polynomials outside of [-1,1]
    s = (2*x - (limit[0] + limit[1])) / (limit[1] - limit[0])

    degree = coefficients.shape[0] + 1
    value = np.float64(1)
    for n in np.arange(1,degree):
        value += chebyshev_single_polynomial(s, n) * coefficients[n-1]

    return value/chebyshev_integral(coefficients, limit, limit)


class Chebyshev(PDF):
    """Chebyshev probability density function of the first kind

    Parameters
    ----------
        name       : str
            The name of the Chebyshev pdf
        observable : parameter.Parameter
            The physical observable
        coefficients : List[parameter.Parameter]
            Arbitrary set of coefficient for f(s) = 1 + coeff_0*T_1(s) + coeff_1*T_2(s) + ...
            s = (2x - (x_min + x_max)) / (x_max - x_min) as Chebyshev polynomials are defined on [-1,1]
            See: https://en.wikipedia.org/wiki/Chebyshev_polynomials

    Attributes
    ----------
        parameters : list[parameter.Parameter]
            parameters contains an arbitrary ordered set of coefficients
    """

    def __init__(self, name: str, observable: Parameter, coefficients: Union[Parameter,List[Parameter]]):
        super().__init__(name, observable)
        # limit number of coefficients to 5
        if len(coefficients) > 5:
            raise ValueError("Chebyshev polynomials of degree > 5 are not supported")
        
        if isinstance(coefficients, Parameter):
            self.coefficients = [coefficients]
        else:
            self.coefficients = coefficients
        self.set_parameters(coefficients)

    def evaluate(self, x: Union[np.float64, np.ndarray]) -> Union[np.float64, np.ndarray]:
        """`chebyshev` is evaluated on `x` with the initialised parameters. Before evaluation, `x`
        is constrained to only those values which are inside of the observable limits

        Parameters
        ----------
            x : np.float64, np.ndarray
                `x` has to be a numeric value or array-like so that the `polynomial` can be evaluated

        Returns
        -------
            np.float64, np.ndarray
                The return value has the same type as `x`
        """
        return chebyshev(x, np.array([coefficient.value for coefficient in self.coefficients], dtype=np.float64), self.observable.limits)

    def integrate(self, ranges: Union[list, tuple, np.ndarray]) -> Union[np.float64, np.ndarray]:
        """
        Integral of the polynomial function between two limit values.

        Parameters
        ----------
            ranges: list, tuple, or np.array
                    `limits` are the boundaries where the integral is computed

        Returns
        -------
            np.array
                The return values have the size of limits -1
        """
        coefficients = np.array([coefficient.value for coefficient in self.coefficients], dtype=np.float64)
        normalisation = chebyshev_integral(coefficients, self.observable.limits, self.observable.limits)
        return chebyshev_integral_vec(coefficients, ranges, self.observable.limits)/normalisation


    @staticmethod
    def func(x: Union[np.float64, np.ndarray], coefficients: Union[np.float64, List[np.float64]], **kwargs) -> Union[
        np.float64, np.ndarray]:
        """Static Version of the `evaluate` method

        Parameters
        ----------
            x           : int, float, np.ndarray
                `x` has to be a numeric value or array-like so that the `double_sided_crystalball` can be evaluated
            coefficients: float, List[float]
                list of coefficient for f(x) = coeff_0 + coeff_1*T_1(x) + coeff_2*T_2(x) + ...
            **kwargs:
                limits  : array-like
                    `limits` has to be of size `2` and is used in order to correctly normalise the
                    chebyshev

        Returns
        -------
            int, float, np.ndarray
                The return value has the same type and shape as `x`
        """
        
        if "limits" in kwargs:
            if not isinstance(kwargs["limits"], (list, tuple, np.ndarray)) or len(kwargs["limits"]) != 2:
                raise TypeError(ErrorMessages.get_error_message(4))
            x = check_data_compatibility(kwargs["limits"], x)
        else:
            kwargs["limits"] = (-1e10, +1e10)

        if np.isscalar(coefficients):
            return chebyshev(x, np.array([coefficients], dtype=np.float64), kwargs["limits"], kwargs["limits"])
        else:
            return chebyshev(x, coefficients, kwargs["limits"], kwargs["limits"])
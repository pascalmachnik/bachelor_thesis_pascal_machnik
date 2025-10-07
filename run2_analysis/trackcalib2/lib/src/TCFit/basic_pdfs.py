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

from typing import Union, Tuple, List

""" TODO:
  - Add generators - sample from distribution and data points should be in the boundary (2nd order problem)
"""

##########################################################################################################
#                                      GLOBAL CONSTANTS                                                  #
##########################################################################################################

PI = np.pi
SQRT2 = np.sqrt(2)

#########################################################################################################
#                                      GAUSSIAN IMPLEMENTATION                                          #
#########################################################################################################

#@numba.njit(fastmath = {'reassoc', 'arcp', 'contract', 'afn' })
def gaussian_integral(mu: np.float64, sigma: np.float64, limits: Union[list, tuple, np.ndarray]) -> np.float64:
    _sigma_sqrt2 = sigma * SQRT2

    x_min = (limits[0] - mu) / _sigma_sqrt2
    x_max = (limits[1] - mu) / _sigma_sqrt2
    scale = (np.sqrt(np.pi) * _sigma_sqrt2)

    integral = np.float64(0.0)
    if x_min * x_max < 0:
        integral = scale * 0.5 * (2 - math.erfc(-x_min) - math.erfc(x_max))
    elif x_max < 0:
        integral = scale * 0.5 * (math.erfc(np.abs(x_max)) - math.erfc(np.abs(x_min)))
    else:
        integral = scale * 0.5 * (math.erfc(np.abs(x_min)) - math.erfc(np.abs(x_max)))

    return integral


#@numba.njit(fastmath = {'reassoc', 'arcp', 'contract', 'afn' }, parallel = True)
def gaussian_integral_vec(mu: np.float64, sigma: np.float64, ranges: Union[list, tuple, np.ndarray]) -> np.ndarray:

    size = ranges.shape[0]
    if size < 2:
        raise ValueError("Size of limits < 2")

    integrals = np.zeros(shape = size-1, dtype = np.float64)
        
    for i in numba.prange(size - 1):
        integrals[i] = gaussian_integral(mu, sigma, [ranges[i], ranges[i+1]])
    
    return integrals


#@numba.njit(fastmath = {'reassoc', 'arcp', 'contract', 'afn' })
def gaussian(x: Union[np.float64, np.ndarray],
             mu: np.float64,
             sigma: np.float64,
             limits: Union[list, tuple, np.ndarray]) -> Union[np.float64, np.ndarray]:
    """Implementation of the gaussian pdf with numba acceleration
    
    Parameters
    ----------
        x     : int, float, np.ndarray
            The `gaussian` pdf is evaluated on x
        mu    : int, float
            The central position of the `gaussian`
        sigma : int, float
            Determines the shape of the `gaussian`
        limits : list, tuple, np.ndarray
            `limit` is used to determine the normalisation of the gaussian

    Returns
    -------
        int, float, np.ndarray
            The return value has the same type as `x`
    """
    z = (x - mu) / sigma
    normalisation = gaussian_integral(mu, sigma, limits)
    return np.exp(-0.5 * z ** 2) / normalisation


class Gaussian(PDF):
    """ Gaussian probability density function

    Parameters
    ----------
        name       : str
            The `name` of the Gaussian pdf
        observable : parameter.Parameter
            Physical observable
        mu         : parameter.Parameter, utils.formula.Formula
            Positional shift of the Gaussian
        sigma      : parameter.Parameter, utils.formula.Formula
            Shape of the Gaussian

    Attributes
    ----------
        parameters : list[parameter.Parameter or utils.formula.Formula]
            `mu` and `sigma` are stored in this parameter list. These are needed such that 
            later during fitting, the parameters get updated and are accessable through the 
            pdf generically
    """

    def __init__(self, name: str, observable: Parameter,
                 mu: Union[Parameter, Formula],
                 sigma: Union[Parameter, Formula]):
        super().__init__(name = name, observable = observable)
        self.mu = mu
        self.sigma = sigma
        self.set_parameters([self.mu, self.sigma])

    def evaluate(self, x: Union[np.float64, np.ndarray]) -> Union[np.float64, np.ndarray]:
        """`gaussian` is evaluated on `x` with the initialised parameters. Before evaluation, `x`
        is constrained to only those values which are inside of the observable limit
        
        Parameters
        ----------
            x : int, float, np.ndarray
                `x` has to be a numeric value or array-like so that the `gaussian` can be evaluated

        Returns
        -------
            int, float, np.ndarray
                The return value has the same type as `x`
        """
        return gaussian(x, self.mu.value, self.sigma.value, self.observable.limits)

    def integrate(self, ranges: Union[list, tuple, np.ndarray]) -> Union[np.float64, np.ndarray]:
        """
        Integral of the gaussian function between two limit values.

        Parameters
        ----------
            ranges: list, tuple, or np.array
                    `limits` are the boundaries where the integral is computed

        Returns
        -------
            np.array
                The return values have the size of limits -1
        """
        normalisation = gaussian_integral(self.mu.value, self.sigma.value, self.observable.limits)
        return gaussian_integral_vec(self.mu.value, self.sigma.value, ranges)/normalisation
        
    @staticmethod
    def func(x: Union[np.float64, np.ndarray], mu: np.float64, sigma: np.float64, **kwargs) -> Union[np.float64, np.ndarray]:
        """Static Version of the `evaluate` method
        
        Parameters
        ----------
            x  : int, float, np.ndarray
                `x` has to be a numeric value or array-like so that the `gaussian` can be evaluated
            mu : int, float
                `mu` determines the central location of the `gaussian`
            sigma : int, float
                `sigma` determines the shape of the `gaussian`
            **kwargs:
                limits : list, tuple, np.ndarray
                    `limits` has to be of size `2` and sets the boundaries for the values `x`

        Returns
        -------
            int, float, np.ndarray
                The return value has the same type as `x`
        
        Raises
        ------
            TypeError
                Raises if limit does not fulfill the correct type, to be array-like with size `2`
        """
        if not "limits" in kwargs:
            raise KeyError(ErrorMessages.get_error_message(12, "limits", "kwargs"))
        if not isinstance(kwargs["limits"], (list, tuple, np.ndarray)) or len(kwargs["limits"]) != 2:
            raise TypeError(ErrorMessages.get_error_message(4))
        x = check_data_compatibility(kwargs["limits"], x)

        return gaussian(x, mu, sigma, kwargs["limits"])


#########################################################################################################
#                                      EXPONENTIAL IMPLEMENTATION                                       #
#########################################################################################################
#@numba.njit(fastmath = {'reassoc', 'arcp', 'contract', 'afn' })
def exponential_integral(slope: np.float64, limits: Union[list, tuple, np.ndarray]) -> np.float64:

    term_1 = np.float64(0.0)
    term_2 = np.float64(0.0)
    if limits[0] > 0:
        term_1 = np.exp(-slope * limits[0])

    if limits[1] > 0:
        term_2 = np.exp(-slope * limits[1])

    integral = np.float64(term_1 - term_2 + 1e-323)
    return integral


#@numba.njit(fastmath = {'reassoc', 'arcp', 'contract', 'afn' }, parallel = True)
def exponential_integral_vec(slope: np.float64, ranges: Union[list, tuple, np.ndarray]) -> np.ndarray:

    size = ranges.shape[0]
    if size < 2:
        raise ValueError("Size of limits < 2")

    integrals = np.zeros(shape = size-1, dtype = np.float64)
        
    for i in numba.prange(size - 1):
        integrals[i] = exponential_integral(slope, [ranges[i], ranges[i+1]])
    
    return integrals

#@numba.njit(fastmath = {'reassoc', 'arcp', 'contract', 'afn' })
def exponential(x: Union[np.float64, np.ndarray], slope: np.float64, limits: Union[list, tuple, np.ndarray]) -> Union[
    np.float64, np.ndarray]:
    """Exponential pdf implementation accelerated by numba

    Parameters
    ----------
        x     : float, np.ndarray
            `x` is evaluated on the exponential pdf
        slope : float
            `slope` determines the shape of the exponential pdf
        limits : array-like
            `limit` has to be of size `2` since it is used for the normalisation
            of the exponential pdf

    Returns
    -------
        float, np.ndarray
            The return value has the same type as `x` and is the computed value for the exponential pdf
    """
    integral = exponential_integral(slope, limits)

    value = np.where(x > 0.0, slope * np.exp(-slope * x) / integral, 0.0)

    return value

class Exponential(PDF):
    """Exponential probability density function

    Parameters
    ----------
        name       : str
            The name of the exponential pdf
        observable : parameter.Parameter
            The physical observable
        slope      : parameter.Parameter, utils.formula.Formula
    
    Attributes
    ----------
        parameters : list[parameter.Parameter]
            `parameters` is a list which contains the `slope` since it is the only parameter
            of the exponential pdf
    """

    def __init__(self, name: str, observable: Parameter, slope: Union[
        Parameter, Formula]):
        super().__init__(name, observable)
        self.slope = slope
        self.set_parameters([self.slope])

    def evaluate(self, x: Union[np.float64, np.ndarray]) -> Union[int, float, np.ndarray]:
        """`exponential` is evaluated on `x` with the initialised parameters. Before evaluation, `x`
        is constrained to only those values which are inside of the observable limit
        
        Parameters
        ----------
            x : int, float, np.ndarray
                `x` has to be a numeric value or array-like so that the `gaussian` can be evaluated

        Returns
        -------
            int, float, np.ndarray
                The return value has the same type as `x`
        """
        return exponential(x, self.slope.value, self.observable.limits)

    def integrate(self, ranges: Union[list, tuple, np.ndarray]) -> Union[np.float64, np.ndarray]:
        """
        Integral of the gaussian function between two limit values.

        Parameters
        ----------
            ranges: list, tuple, or np.array
                    `limits` are the boundaries where the integral is computed

        Returns
        -------
            np.array
                The return values have the size of limits -1
        """
        normalisation = exponential_integral(self.slope.value, self.observable.limits)
        return exponential_integral_vec(self.slope.value, ranges)/normalisation

    @staticmethod
    def func(x: Union[np.float64, np.ndarray], slope: np.float64, limits: np.float64,
             **kwargs) -> Union[np.float64, np.ndarray]:
        """Static Version of the `evaluate` method
        
        Parameters
        ----------
            x     : int, float, np.ndarray
                `x` has to be a numeric value or array-like so that the `exponential` can be evaluated
            slope : int, float
                `slope` determines the shape of the `exponential` pdf
            limits : array-like
                `limits` is used to check the data compatibility, limits is an interval and has to be of size
                `2`
            **kwargs, no options

        Returns
        -------
            int, float, np.ndarray
                The return value has the same type as `x`
        """
        return exponential(x, slope, limits)


#########################################################################################################
#                                      CRYSTALBALL IMPLEMENTATION                                       #
#########################################################################################################

@numba.njit(fastmath = {'reassoc', 'arcp', 'contract', 'afn' })
def crystalball_integral(mu: np.float64, 
                         sigma: np.float64, 
                         alpha: np.float64, 
                         n: np.float64, 
                         limits: Union[list, tuple, np.array]) -> np.float64:

    abs_sigma = np.abs(sigma)
    abs_alpha = np.abs(alpha)

    value = np.float64(0.0)

    use_log = False
    if np.abs(n - 1.0) < 1e-5:
        use_log = True

    z_min = (limits[0] - mu) / abs_sigma
    z_max = (limits[1] - mu) / abs_sigma

    if alpha < 0:
        tmp = z_min
        z_min = -z_max
        z_max = -tmp

    if z_min >= - abs_alpha:
        value += abs_sigma * np.sqrt(PI / 2) * (math.erf(z_max / SQRT2) - math.erf(z_min / SQRT2))
    elif z_max < - abs_alpha:
        a = np.power((n / abs_alpha), n) * np.exp(-0.5 * abs_alpha * abs_alpha)
        b = n / abs_alpha - abs_alpha
        if (use_log):
            value += a * abs_sigma * (np.log(b - z_min) - np.log(b - z_max))
        else:
            value += a * abs_sigma / (1.0 - n) * (
                    1.0 / np.power(b - z_min, n - 1.0) - 1.0 / np.power(b - z_max, n - 1.0))
    else:
        a = np.power((n / abs_alpha), n) * np.exp(-0.5 * abs_alpha * abs_alpha)
        b = n / abs_alpha - abs_alpha

        first_term = np.float64(0.0)
        if use_log:
            first_term += a * abs_sigma * (np.log(b - z_min) - np.log(n / abs_alpha))
        else:
            first_term += a * abs_sigma / (1.0 - n) * (
                    1.0 / np.power(b - z_min, n - 1.0) - 1.0 / np.power(n / abs_alpha, n - 1.0))

        second_term = abs_sigma * np.sqrt(PI / 2) * (
                    math.erf(z_max / SQRT2) - math.erf(-abs_alpha / SQRT2))

        value += first_term + second_term

    #python3.9 trackcalib.py fit -year "2012_50ns_strip" -max_entries 50000 -v -sim_fit -method "Long" -mode "MC" -variables_2D "P-ETA" -variables ""


    if value > 1e-323: return value
    else:              return np.float64(1e-38)


@numba.njit(fastmath = {'reassoc', 'arcp', 'contract', 'afn' }, parallel = True)
def crystalball_integral_vec(mu: np.float64, 
                             sigma: np.float64, 
                             alpha: np.float64, 
                             n: np.float64,
                             ranges: Union[list, tuple, np.ndarray]) -> np.ndarray:

    size = ranges.shape[0]
    if size < 2:
        raise ValueError("Size of limits < 2")

    integrals = np.zeros(shape = size-1, dtype = np.float64)
        
    for i in numba.prange(size - 1):
        integrals[i] = crystalball_integral(mu, sigma, alpha, n, [ranges[i], ranges[i+1]])
    
    return integrals

# The fastmath optins are needed in order to be able to check nan and inf.
# If just set to fastmath=True, numba does not return inf/nan as infs and nans,
# leading to problems when doing sanity checks on the parameters.
# The parameters can be nan as iminuit can return nan.

@numba.njit(fastmath = {'reassoc', 'arcp', 'contract', 'afn' }) 
def crystalball(x: Union[np.float64, np.ndarray],
                mu: np.float64, 
                sigma: np.float64, 
                alpha: np.float64, 
                n: np.float64, 
                limits: Union[list, tuple, np.array]) -> Union[np.float64, np.ndarray]:
                
    """Implementation of the crystal ball function

    Parameters
    ----------
        x : float
            `x` is evaluated on the crystal ball function
        mu    : float
            Determines the central location of the crystalball function
        sigma : float
            Determines the shape of the crystalball function
        alpha : float
            Controls the central shape
        n     : float
            Controls the shape of the tails
        low   : float
            `low` is used to determine the minimum z value
        high  : float
            `high` is used to determine the maximum z value

    Returns
    -------
        int, float, np.ndarray
            The return value has the same type as `x`
    """
    
    sign_alpha = np.sign(alpha)
    abs_alpha = np.abs(alpha)

    #Sanity check for NaNs
    if (np.isnan(mu) or np.isnan(sigma) or np.isnan(n) or np.isnan(alpha)):
        return x - x + 1e-323
    #Check the parameter values
    if ((abs_alpha <= 1e-12) or (sigma == 0.0) or (n < 1.0)):
        return x - x + 1e-323 #Returns an array with 1e-323 everywhere

    z = (x - mu) / sigma * sign_alpha

    a = np.power((n / abs_alpha), n) * np.exp(-0.5 * abs_alpha*abs_alpha)

    b = (n / abs_alpha) - abs_alpha

    value = np.where(z > -abs_alpha, np.exp(-0.5 * z * z), np.divide(a, np.power(b - z*sign_alpha, n)))

    integral  = crystalball_integral(mu, sigma, alpha, n, limits)
    if integral > 1e-323:
        value = value / integral
    else:
        value = value - value + 1e323

    return value


class CrystalBall(PDF):
    """CrystalBall probability density function

    Parameters
    ----------
        name : str
            The name of the crystall ball function
        observable : parameter.Parameter
            Physical observable
        mu         : parameter.Parameter, utils.formula.Formula
            Positional shift of the crystall ball
        sigma      : parameter.Parameter, utils.formula.Formula
            Shape of the crystal ball
        alpha      : parameter.Parameter, utils.formula.Formula
            Tails of the crystal ball
        n          : parameter.Parameter, utils.formula.Formula
            Tails of the crystall ball

    Attributes
    ----------
        parameters : list[parameter.Parameter or utils.formula.Formula]
            `parameters` is a list which contains `mu`, `sigma`, `alpha` and `n`. This list is used in order 
            to keep tracks and simply access to the parameters of the crystal ball
    """

    def __init__(self, name: str, observable: Parameter,
                 mu: Union[Parameter, Formula], sigma: Union[
                Parameter, Formula],
                 alpha: Union[Parameter, Formula], n: Union[
                Parameter, Formula]):
        super().__init__(name, observable)
        self.mu = mu
        self.sigma = sigma
        self.alpha = alpha
        self.n = n
        self.set_parameters([self.mu, self.sigma, self.alpha, self.n])

    def evaluate(self, x: Union[np.float64, np.ndarray]) -> Union[np.float64, np.ndarray]:
        """`crystalball` is evaluated on `x` with the initialised parameters. Before evaluation, `x`
        is constrained to only those values which are inside of the observable limit
        
        Parameters
        ----------
            x : int, float, np.ndarray
                `x` has to be a numeric value or array-like so that the `crystalball` can be evaluated

        Returns
        -------
            int, float, np.ndarray
                The return value has the same type as `x`
        """

        return crystalball(x,
                           self.mu.value,
                           self.sigma.value,
                           self.alpha.value,
                           self.n.value,
                           self.observable.limits)

    def integrate(self, ranges: Union[list, tuple, np.ndarray]) -> Union[Tuple[float, float], Tuple[List[float], List[float]]]:

        normalisation = crystalball_integral(self.mu.value,
                                             self.sigma.value,
                                             self.alpha.value,
                                             self.n.value, 
                                             self.observable.limits)
        return crystalball_integral_vec(self.mu.value,
                                        self.sigma.value,
                                        self.alpha.value,
                                        self.n.value, 
                                        ranges)/normalisation

    @staticmethod
    def func(x: Union[np.float64, np.ndarray], mu: np.float64, sigma: np.float64,
             alpha: np.float64, n: np.float64, **kwargs) -> Union[np.float64, np.ndarray]:
        """Static Version of the `evaluate` method
        
        Parameters
        ----------
            x          : int, float, np.ndarray
                `x` has to be a numeric value or array-like so that the `crystalball` can be evaluated
            mu         : int, float
                Positional shift of the crystall ball
            sigma      : int, float
                Shape of the crystal ball
            alpha      : int, float
                Tails of the crystal ball
            n          : int
                Tails of the crystall ball
            **kwargs:
                limits : array-like
                    `limits` has to be of size `2` and is used in order to correctly normalise the
                    crystalball

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
            kwargs["limits"] = (-1e10 * sigma - mu, +1e10 * sigma + mu)


        return crystalball(x, mu, sigma, alpha, n, kwargs["limits"])


#########################################################################################################
#                           POLYNOMIAL IMPLEMENTATION                                                   #
#########################################################################################################
@numba.njit(fastmath = {'reassoc', 'arcp', 'contract', 'afn' })
def polynomial_integral(coefficients: np.ndarray, 
                        limits: Union[list, tuple, np.ndarray]) -> np.float64:

    degree = coefficients.shape[0]
    value = np.float64(0.0)
    for n in range(0,degree):
        value += (np.power(limits[1], n+1) - np.power(limits[0], n+1)) * coefficients[n] /(n+1)

    return value


@numba.njit(fastmath = {'reassoc', 'arcp', 'contract', 'afn' }, parallel = True)
def polynomial_integral_vec(coefficients: np.ndarray, 
                            ranges: Union[list, tuple, np.ndarray]) -> np.ndarray:

    size = ranges.shape[0]
    if size < 2:
        raise ValueError("Size of limits < 2")

    integrals = np.zeros(shape = size-1, dtype = np.float64)
        
    for i in numba.prange(size - 1):
        integrals[i] = polynomial_integral(coefficients, [ranges[i], ranges[i+1]])
    
    return integrals
    

@numba.njit(fastmath = {'reassoc', 'arcp', 'contract', 'afn' })
def polynomial(x: Union[np.float64, np.ndarray], 
               coefficients: np.ndarray, 
               limits: Union[list, tuple, np.ndarray]) -> Union[np.float64, np.ndarray]:

    degree = coefficients.shape[0]
    value = (x - x) + np.float64(coefficients[0])
    for n in np.arange(1,degree):
        value += np.power(x,n) * coefficients[n]

    return value/polynomial_integral(coefficients, limits)


class Polynomial(PDF):
    """Polynomial probability density function

    Parameters
    ----------
        name       : str
            The name of the double-sided crystalball pdf
        observable : parameter.Parameter
            The physical observable
        coefficients : List[parameter.Parameter]
            Arbitrary set of coefficient for f(x) = coeff_0 + coeff_1*x + coeff_2*x*x + ...


    Attributes
    ----------
        parameters : list[parameter.Parameter]
            parameters contains an arbitrary ordered set of coefficients
    """

    def __init__(self, name: str, observable: Parameter, coefficients: Union[Parameter,List[Parameter]]):
        super().__init__(name, observable)
        if isinstance(coefficients, Parameter):
            self.coefficients = [coefficients]
        else:
            self.coefficients = coefficients
        self.set_parameters(coefficients)

    def evaluate(self, x: Union[np.float64, np.ndarray]) -> Union[np.float64, np.ndarray]:
        """`double_sided_crystalball` is evaluated on `x` with the initialised parameters. Before evaluation, `x`
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
        return polynomial(x, np.array([coefficient.value for coefficient in self.coefficients], dtype=np.float64), self.observable.limits)

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
        normalisation = polynomial_integral(coefficients, self.observable.limits)
        return polynomial_integral_vec(coefficients, ranges)/normalisation


    @staticmethod
    def func(x: Union[np.float64, np.ndarray], coefficients: Union[np.float64, List[np.float64]], **kwargs) -> Union[
        np.float64, np.ndarray]:
        """Static Version of the `evaluate` method

        Parameters
        ----------
            x           : int, float, np.ndarray
                `x` has to be a numeric value or array-like so that the `double_sided_crystalball` can be evaluated
            coefficients: float, List[float]
                list of coefficient for f(x) = coeff_0 + coeff_1*x + coeff_2*x*x + ...
            **kwargs:
                limits  : array-like
                    `limits` has to be of size `2` and is used in order to correctly normalise the
                    double-sided crystalball

        Returns
        -------
            int, float, np.ndarray
                The return value has the same type and shape as `x`
        """

        if np.isscalar(coefficients):
            return polynomial(x, np.array([coefficients], dtype=np.float64))
        else:
            return polynomial(x, coefficients)


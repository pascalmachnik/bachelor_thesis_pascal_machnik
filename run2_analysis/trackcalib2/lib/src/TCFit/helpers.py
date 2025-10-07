from typing import Union, Tuple

import numba as nb
import numpy as np
from scipy import stats as st

xspace = np.linspace(0, 100, 10000)
one_sigma = st.norm.cdf(-1)


@nb.njit(parallel = True)
def x_log_x_y(x: Union[np.ndarray, np.float64], y: Union[np.array, np.float64]) -> Union[np.ndarray, np.float64]:
    return np.dot(x, np.log(x)+np.log(y))


@nb.njit(parallel = True)
def x_log_x_over_y(x: Union[np.ndarray, np.float64], y: Union[np.array, np.float64]) -> Union[np.ndarray, np.float64]:
    return np.dot(x, np.log(x)-np.log(y))


@nb.njit(parallel = True)
def xmin(prob):
    return xspace[np.argmin(np.abs(prob - one_sigma))]


@nb.njit(parallel = True)
def xmax(prob):
    return xspace[np.argmin(np.abs(prob-one_sigma))]


def freq_error(k: Union[np.float64, float, int]) -> Tuple[np.float64]:
    if k < 20:
        prob = st.poisson.sf(k-1, xspace)
        minimum = xmin(prob)-k

        prob = st.poisson.cdf(k,xspace)
        maximum = xmax(prob)-k
    else:
        maximum = np.sqrt(k)
        minimum = -maximum
    return (minimum, maximum)
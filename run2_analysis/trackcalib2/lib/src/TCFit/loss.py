import abc
from abc import abstractmethod
from typing import Tuple, List, Sequence, Union

import numpy as np
from numba import jit
from scipy.special import loggamma
from ..TCFit.model import Model
from ..TCFit.basic_pdfs import PDF
from ..TCFit.error_messages import ErrorMessages


def check_dimensions(bins: Sequence[Union[List[float], float]], histo: Sequence[Union[List[float], float]]) -> bool:
    """Check dimension of bins contents and binned dataset"""
    shape_b = np.shape(bins)
    shape_h = np.shape(histo)
    shape_dim = np.shape(shape_h)
    if shape_b[0] == shape_dim[0]:
        for index, sub_bins in enumerate(bins):
            if len(sub_bins) - 1 != shape_h[index]:
                return False

        return True
    return False


class IFLoss(metaclass=abc.ABCMeta):

    @abstractmethod
    def __call__(self, *args, **kwargs):
        pass

    def __set_parameters(self, *parameters):
        self.parameters = []
        for parameter in parameters:
            self.parameters.append(parameter)

        self.names = [parameter.name for parameter in self.parameters]
        self.values = [parameter.value for parameter in self.parameters]
        self.limits = [parameter.limits for parameter in self.parameters]
        self.errors = [parameter.error for parameter in self.parameters]
        self.fixed = [parameter.fixed for parameter in self.parameters]

    @abstractmethod
    def __update_parameters(self, *args, **kwargs):
        pass

class UnbinnedLikelihood(IFLoss):

    def __init__(self, data, model, extended=False):
        self.data, self.model = self.__set(data, model)
        self.extended = extended
        self.errordef = 0.5
        self.__set_parameters(model.parameters)

    def __call__(self, *values):
        #Update the parameters

        negative_log_likelihood = 0.0
        sum_of_normalisations = 0.0

        likelihoods = self.model.evaluate(self.data)
        for likelihood in likelihoods:
            if any(likelihood <= 0.0):
                return 1e100
            negative_log_likelihood += np.log(likelihood).sum()

        if (self.extended):
            for normalisation in self.model.normalisations:
                sum_of_normalisations += normalisation.value

        return 2. * (sum_of_normalisations - negative_log_likelihood + 1e-323)

    @staticmethod
    def __set(data, model):
        if isinstance(model, Model):
            if len(data) == len(Model.pdfs):
                return data, model
            else:
                raise ValueError("Data not compatible with Model used")

        else:
            return data, model


class BinnedLikelihood(IFLoss):

    def __init__(self, data, bin_container, model, extended=False):
        self.data, self.bin_container, self.model = self.__set(data, bin_container, model)
        self.extended = extended
        self.errordef = 0.5

    def __call__(self, *values):

        negative_log_likelihood = 0.0

        expected_container = self.model.integrate(bin_container=self.bin_container)
        for expected, frequency in zip(expected_container, self.data):
            if np.any(np.array(expected) <= 0):
                return 1e323
            negative_log_likelihood += np.sum(frequency * np.log(expected) - expected - loggamma(frequency + 1))
            if np.any(negative_log_likelihood) == 0:
                return 1e323

        return -1.0 * (negative_log_likelihood + 1e-323)

    @staticmethod
    def __set(data, bin_container, model):
        # TODO: check model dimensions
        if isinstance(model, PDF):
            # TODO: check dimensions of the PDF (1D, 2D or 3D pdf)
            pass
        elif isinstance(model, Model):
            # TODO: check if categories are present
            pass
        else:
            raise TypeError(ErrorMessages.get_error_message(11, "model", "PDF or Model"))

        if check_dimensions(bin_container, data):
            return data, bin_container, model
        else:
            raise TypeError(ErrorMessages.get_error_message(7, "bins container", "binned data"))


class binned_chi2(IFLoss):
    pass


class regression(IFLoss):
    pass

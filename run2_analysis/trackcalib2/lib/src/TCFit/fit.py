import abc
from abc import abstractmethod
from typing import Tuple, List, Sequence, Union

import iminuit
import numpy as np

from ..TCFit.error_messages import ErrorMessages
from ..TCFit.model import Model
from ..TCFit.monitoring import MonitorLikeLihood
from ..TCFit.parameter import Parameter
from ..TCFit.basic_pdfs import PDF
from ..TCFit.loss import IFLoss
from ..TCFit.dataset import check_data_compatibility
from scipy.special import loggamma

##########################################################################################################
#                              Interface LikelihoodFit Implementation                                    #
##########################################################################################################


class IFFit(metaclass=abc.ABCMeta):
    """This interface method defines which APIs should be used
    if the user want to fit the given input data with the user defined pdf models

    Parameters
    ----------
        data   : array-like
            `data` contains the data points which will be used in iminuit in order to fit
            the pdf models to it
        model  : pdfs.model.Model
            `model` contains the pdfs, one pdf for each category

    Optionals
    ---------
        bins           : int, List[float], np.arrays, List[List[arrays]
            number of bins in case of a int otherwise interpreted as a list of boundaries
        binned_dataset : bool
            set True in case the dataset passed is a binned dataset

    Methods (abstract)
    ------------------
        minimize() iminuit.Minuit
            The definition of the minimisations take place here, also the Minuit
            object is defined which holds the necessary configurations for the fitting
            process
    """

    @abstractmethod
    def minimise(self) -> iminuit.Minuit:
        """When `data` and `model` are ready to be fitted, the user can start
        the minimization process with this method, iminuit is used in order to 
        do the fits

        Returns
        -------
            iminuit.Minuit
                The `Minuit` object is returned which contains all the relevant information
                when it is done with the fitting, e.g. the minimized parameters, the calculated
                errors, etc.
        """
        pass


##########################################################################################################
#                                  LikelihoodFit Implementation                                          #
##########################################################################################################

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


class LikelihoodFit(IFFit):
    """The model is a representation of the fit model which is used 
    in order to perform a fit on the data, this class uses the negative
    log-likelihood in order to estimate the optimal parameters
    
    Parameters
    ----------
        data   : array-like
            `data` is the set of data points on which the model will be fitted on. The size of the `data`
            array has to be the same as the `model` pdfs array since each data is used for each pdf
        model  : pdfs.model.Model
            `model` contains the pdfs which are used as fit models in order to perform the fit on the `data` 

    Attributes
    ----------
        parameters : list[utils.parameter.Parameter]
            The model holds parameters, these are fetched because these contain all 
            the relevant information for iminuit in order to start the fitting process
        minuit     : iminuit.Minuit
            The object Minuit holds all the information about the fitting procedure as e.g.
            the estimated errors, parameter values, ...
        names      : list[str]
            These are the parameter names
        values     : list[int, float]
            These are the parameter values
        limits     : list[tuple]
            These are the parameter limits on the values
        errors     : list[int, float]
            These are the parameter errors on the values
        fixed      : list[bool]
            These determine if the parameters are hold fixed during the fitting procedure

    Methods
    -------
        minimize(hesse : bool, minos : bool, print_level : int) iminuit.Minuit
            If `hesse` is set, then the uncertainties will be computed based on a covariance
            matrix for all fitted parameters in order to also respect the correlation between
            fitting parameters
            If `minos` is set, then confidence intervals for the fitting parameters will 
            be approximately computed, no parameter correlation is taken into account here
            One can set `print_level` to ...
    """

    def __init__(self, data: Sequence[Union[float,List[float],np.ndarray]], model: Union[Model,PDF], monitor_on: bool = False, **kwargs):
        self.__check_compatibility_between_data_and_model(data, model)

        self.model = model
        self.parameters = self.__set_parameters(model)
        if "print_level" in kwargs and isinstance(kwargs["print_level"], int):
            self.print_level = kwargs["print_level"]
        else:
            self.print_level = 1
        self.__monitor = MonitorLikeLihood(parameters=self.parameters, monitor_on=monitor_on)
        self.binned = False
        self.binned_dataset = False
        if "binned_dataset" in kwargs:
            self.binned_dataset = kwargs["binned_dataset"]

        self.data = []
        if self.binned_dataset == False:
            if isinstance(self.model, Model):
                for _data, pdf in zip(data, self.model.pdfs):
                    self.data.append(np.array(check_data_compatibility(pdf.observable.limits, _data), dtype=np.float64))
            if isinstance(self.model, PDF):
                self.data = np.array(check_data_compatibility(self.model.observable.limits, data), dtype = np.float64)

        else:
            self.data = np.array(data, dtype=np.float64)
        self.minuit = self.__prepare_settings()
        if "bins" in kwargs:
            __bins = kwargs['bins']
            if isinstance(__bins, int) and self.binned_dataset is True:
                raise ValueError("If dataset is binned, set of boundaries has to be provided")

            self.binned = True
            self.frequencies, self.bin_container = self.__set_bins(kwargs["bins"])


    def minimise(self, hesse: bool = False, minos: bool = False, strategy: int = 1) -> iminuit.Minuit:
        """With starting this method, the fitting procedure will be triggered and the model will be
        fitted on the data

        Parameters
        ----------
            hesse       : bool
                If `hesse` is `True`, then the uncertainties will be computed based on the second-derivative
                matrix with respect to all fitted parameters. The error matrix is obtained by inverting the
                second-derivative. If matrix is not positive defined a warning message is printed.
            minos       : bool
                If `minos` is `True`, then confidence intervals for the fitting parameters will 
                be accurately computed. It takes into account parameter correlations.
        
        Returns
        -------
            iminuit.Minuit
                The object Minuit holds all the information about the fitting procedure as e.g.
                the estimated errors, parameter values, ...
        """
        self.minuit.strategy = strategy

        print(self.minuit.migrad())
        valid = self.minuit.valid

        if not valid:
            self.minuit.migrad(ncall=4000, iterate=4)

        if hesse:
            print(self.minuit.hesse())
            if self.minuit.valid:
                self.minuit.migrad()
                print(self.minuit.hesse())

        if minos and self.minuit._fmin.is_valid:
            print(self.minuit.minos())

        # final update on parameters
        for parameter, value, error in zip(self.parameters, self.minuit.values, self.minuit.errors):
            self.__monitor(name=parameter.name, value=value)
            parameter.value = value
            parameter.error = error

        if minos:
            for parameter, value in zip(self.parameters, self.minuit.merrors.values()):
                parameter.asym_errors = (value.lower, value.upper)
        self.__monitor.display_parameter_evolution(epochs=self.minuit.nfcn + 2)

        return self.minuit

    def __set_parameters(self, model: Model) -> List[Parameter]:
        """The parameters are extracted from the model

        Parameters
        ----------
            model : pdfs.model.Model
                The model which were given as an input from the user. The model holds
                all the parameters (remember that normalisations are in the list of parameters)

        Returns
        -------
            list[utils.parameter.Parameter]
                The filtered parameters are returned
        """
        parameters = []
        for parameter in model.parameters:
            parameters.append(parameter)

        return parameters

    def __prepare_settings(self) -> iminuit.Minuit:
        """The settings consist of all the parameter information, the values, the
        names, the limits, ..., these are extracted and with them the Minuit object
        will be instantiated which will be needed to start the fitting procedure

        Returns
        -------
            iminuit.Minuit
                The Minuit object will be used to start the fitting procedure, it already
                holds all the relevant information for fitting, as a loss function, the negative
                log-likelihood will be used
        """
        self.names = [parameter.name for parameter in self.parameters]
        self.values = [parameter.value for parameter in self.parameters]
        self.limits = [parameter.limits for parameter in self.parameters]
        self.errors = [parameter.error for parameter in self.parameters]
        self.fixed = [parameter.fixed for parameter in self.parameters]

        minuit = iminuit.Minuit(self.__evaluate_likelihood, self.values, name=self.names)
        minuit.limits = self.limits
        minuit.errors = self.errors
        minuit.fixed = self.fixed
        minuit.errordef = 0.5
        minuit.print_level = self.print_level

        return minuit

    def __evaluate_likelihood(self, *values):
        """The negative log-likelihood will be computed
        for the given model and data, the result will indicate
        iminuit how to change the parameters

        Parameters
        ----------
            *values : *list[int, float]
                The parameter values are given as input for the computation
                of the negative log-likelihood
        """
        # First step: Update the parameters
        for parameter, value in zip(self.parameters, values[0]):
            self.__monitor(name=parameter.name, value=value)
            parameter.value = value

        # Second step: Calculate the total negative log-likelihood
        negative_log_likelihood = np.float64(0.0)
        sum_of_normalisations = np.float64(0.0)
        result = np.float64(0)

        if not self.binned:
            likelihoods = self.model.evaluate(self.data)
            for likelihood in likelihoods:
                if any(likelihood <= 0.0):
                    return 1e100
                negative_log_likelihood += np.log(likelihood).sum()

            if self.model.normalisations is not None:
                for norm in self.model.normalisations:
                    sum_of_normalisations += norm.value

            result = np.float64(2.0 * (sum_of_normalisations - negative_log_likelihood + 1e-38))

        else:

            expected_container = self.model.integrate(bin_container=self.bin_container)
            for expected, frequency in zip(expected_container, self.frequencies):

                if np.any(np.array(expected) <= 0):
                    return 1e100
                negative_log_likelihood += np.sum(frequency * np.log(expected) - expected - loggamma(frequency + 1))

            result = np.float64(-1.0 * (negative_log_likelihood + 1e-38))


        self.__monitor(name="Norm", value=sum_of_normalisations)
        self.__monitor(name="NLL", value=negative_log_likelihood)
        self.__monitor(name="Loss", value=result)

        return result

    def __check_compatibility_between_data_and_model(self, datasets: Sequence[np.ndarray], model: Model):
        """For the subsequent section we call `data`, `datasets` since `data` is a 2-dimensional array, 
        providing multiple data for each pdf in model because ach pdf in model belongs to a distinct category.
        For that reason we need to check the dimensional requirements, the length of `datasets` has to match the length 
        of the pdfs inside of `model`
        
        Parameters
        ----------
            datasets : array-like[np.ndarray]
                `datasets` contains multiple data for each pdf in model

            model    : pdfs.model.Model
                Each model contains n pdfs which belong to n distinct categories, each of these
                pdfs are evaluated on their appropriate data

        Raises
        ------
            TypeError
                Is raised when the length of `datasets` is not equal to the length of the pdfs inside of `model`
        """
        if len(datasets) != len(model.pdfs):
            raise TypeError(ErrorMessages.get_error_message(7, "data", "model pdfs"))

    def __set_bins(self, bins: Union[int, Sequence[list]]) -> Tuple[list, list]:
        """The user has the possibility to give a sequence of sets of bins,
        these bins will be used in order to determine the frequencies, if the user
        decides to give us no bins, then we will automatically choose by default a binning
        of 10 with equal widths
        
        Parameters
        ----------
            bin_container : array-like[list] or int
                The `bin_container` contains a list of bins, each list is for each 
                dataset in `data`
                If `bin_container` is an int, it is considered as the number of bins

        Returns
        -------
            array-like[list], array-like[list[tuple]]
                A `frequency_container` and a `bin_container` is returned, both have the same
                length as `data` and contain the frequencies for each bin and the binning respectively

        Raises
        ------
            TypeError
                Is raised when the length of the bin_container does not match the length of `data`, each dataset
                needs a defined binning
        """
        if isinstance(bins, int):
            bins_container = []
            frequency_container = []
            if not self.binned_dataset:
                for dataset in self.data:
                    _limits = self.model.observable.limits
                    _bins = np.linspace(_limits[0], _limits[1], bins + 1)
                    frequencies, _ = np.histogram(dataset, bins=_bins)
                    bins_container.append(_bins)
                    frequency_container.append(frequencies)

                return frequency_container, bins_container

            else:
                raise TypeError(ErrorMessages.get_error_message(16))

        elif isinstance(bins, list):
            if self.binned_dataset:
                if not check_dimensions(bins, self.data):
                    raise TypeError(ErrorMessages.get_error_message(7, "bins container", "binned data"))
                else:
                    return self.data, bins
            else:
                frequency_container = []
                bins_container = []
                if len(bins) == len(self.data):
                    for index, sub_bins in enumerate(bins):
                        _limits = self.model.observable.limits
                        if isinstance(sub_bins, int):
                            _bins = np.linspace(_limits[0], _limits[1], sub_bins + 1)
                        elif isinstance(sub_bins, list):
                            if bins[0] >= _limits[0] and bins[-1] < _limits[1]:
                                _bins = sub_bins
                            else:
                                raise TypeError(ErrorMessages.get_error_message(7, "bins container", "binned data"))
                        else:
                            raise TypeError(ErrorMessages.get_error_message(7,"bins"))
                        _frequencies, _ = np.histogram(self.data[index], bins=_bins)

                        bins_container.append(_bins)
                        frequency_container.append(_frequencies)

                return frequency_container, bins_container

        else:
            raise TypeError(ErrorMessages.get_error_message(14, "bins"))




class Fit(IFFit):


    def __init__(self, loss: IFLoss, **kwargs):
        """Class to bound loss classes to the minimisation algorithm

         Parameters
         ----------
             loss   : TCFit.IFLoss
                 `loss` function that is minimised. It already include the information of the model used and the dataset

         Attributes
         ----------
             loss          : TCFit.Loss
                 Loss function
             minimiser     : iminuit.Minuit
                 The object Minuit holds all the information about the fitting procedure as e.g.
                 the estimated errors, parameter values, ...

         Methods
         -------
             minimise(hesse : bool, minos : bool, print_level : int) iminuit.Minuit
                 If `hesse` is set, then the uncertainties will be computed based on a covariance
                 matrix for all fitted parameters in order to also respect the correlation between
                 fitting parameters
                 If `minos` is set, then confidence intervals for the fitting parameters will
                 be approximately computed, no parameter correlation is taken into account here
                 One can set `print_level` to ...
         """

        self.loss = loss
        self.__minimiser = self.__set_minimiser()

    def __set_minimiser(self):

        minuit = iminuit.Minuit(self.__evaluate_likelihood, self.loss.values, self.loss.names)
        minuit.limits = self.loss.limits
        minuit.errors = self.loss.errors
        minuit.fixed = self.loss.fixed
        minuit.errordef = self.loss.errordef
        minuit.print_level = 1

        return minuit

    def __evaluate_likelihood(self, *values):
        """The negative log-likelihood will be computed
        for the given model and data, the result will indicate
        iminuit how to change the parameters

        Parameters
        ----------
            *values : *list[int, float]
                The parameter values are given as input for the computation
                of the negative log-likelihood
        """
        # First step: Update the parameters
        result = self.loss(values)

#        self.__monitor(name="Loss", value=result)

        return result

    def minimise(self, hesse: bool = False, minos: bool = False, print_level: int = 1, **kwargs):
        """With starting this method, the fitting procedure will be triggered and the model will be
        fitted on the data

        Parameters
        ----------
            hesse       : bool
                If `hesse` is `True`, then the uncertainties will be computed based on the second-derivative
                matrix with respect to all fitted parameters. The error matrix is obtained by inverting the
                second-derivative. If matrix is not positive defined a warning message is printed.
            minos       : bool
                If `minos` is `True`, then confidence intervals for the fitting parameters will
                be accurately computed. It takes into account parameter correlations.
            print_level : int
                Minuit printing level

        Returns
        -------
            self
                The object Minuit holds all the information about the fitting procedure as e.g.
                the estimated errors, parameter values, ...
        """
        print(self.__minimiser.migrad())
        if hesse:
            print(self.__minimiser.hesse())
        if minos:
            print(self.__minimiser.minos())

        # final update on parameters
        for parameter, value, error in zip(self.loss.parameters, self.__minimiser.values, self.__minimiser.errors):
            parameter.value = value
            parameter.error = error

        if minos:
            for parameter, value in zip(self.loss.parameters, self.__minimiser.merrors.values()):
                parameter.asym_errors = (value.lower, value.upper)

        return self

    def simplex(self):
        self.__minimiser.simplex()
        for parameter, value, error in zip(self.loss.parameters, self.__minimiser.values, self.__minimiser.errors):
            parameter.value = value
            parameter.error = error

        return self

    def migrad(self):
        self.__minimiser.migrad()
        for parameter, value, error in zip(self.loss.parameters, self.__minimiser.values, self.__minimiser.errors):
            parameter.value = value
            parameter.error = error

        return self

    def hesse(self):
        self.__minimiser.hesse()
        for parameter, value, error in zip(self.loss.parameters, self.__minimiser.values, self.__minimiser.errors):
            parameter.value = value
            parameter.error = error

        return self

    def minos(self, **kwargs):
        self.__minimiser.minos(**kwargs)
        for parameter,  value, error, asym_error in zip(self.loss.parameters,
                                                        self.__minimiser.values,
                                                        self.__minimiser.errors,
                                                        self.__minimiser.merrors.values()):
            parameter.value = value
            parameter.error = error
            parameter.asym_errors = (asym_error.lower, asym_error.upper)

        return self

    def __str__(self):
        return self.__minimiser.__str__()

    def __repr__(self):
        return self.__minimiser.__repr__()
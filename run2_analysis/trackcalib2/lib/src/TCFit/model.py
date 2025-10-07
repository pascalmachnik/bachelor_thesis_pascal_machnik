import abc
from abc import abstractmethod
from typing import Union, List, Tuple, Sequence

import numpy as np

from ..TCFit.basic_pdfs import PDF
from ..TCFit.error_messages import ErrorMessages
from ..TCFit.exceptions import ValidationError
from ..TCFit.extensions import ExtendedPDF, RapidExtendedPDFGenerator
from ..TCFit.parameter import Parameter
from ..TCFit.formula import Formula
from ..TCFit.standards import OrderedSet, String

RANGE_TYPE = "ndarray, list, tuple or array"

def remove_duplicate_attributes(attributes: List[Parameter]) -> List[
    Parameter]:
    """Attributes here are understood as either parameters or
    normalisations, duplicated attributes are removed by this function

    Parameters
    ----------
        attributes : list[utils.parameter.Parameter]
            `attributes` can be either a list of parameters or normalisations

    Returns
    -------
        list[utils.parameter.Parameter]
            The returned list contains attributes which are all unique
    """
    names_encountered = []
    filtered_attributes = []
    for attribute in attributes:
        if attribute.name not in names_encountered:
            filtered_attributes.append(attribute)
            names_encountered.append(attribute.name)

    return filtered_attributes


##########################################################################################################
#                                      Add PDF Implementation                                            #
##########################################################################################################

class AddPDF(PDF):
    """A new pdf will be constructed by summing up several pdfs
    There are 3 cases which determine how this class works
    1st case:
    The input has to be extended pdfs with no fractions given
    2nd case:
    The input has to be unextended pdfs with normalisations given.
    3rd case:
    The input has to be unextended pdfs with fractions given. These fractions
    will be multiplied with the corresponding pdfs
    
    Parameters
    ----------
        name         : str
            Name of the AddPDF instance
        pdfs         : list[standard.IFPDF]
            The pdf container which lists the pdfs which should be summed up
        fractions    : list[utils.parameter.Parameter], optional
            If the user does not provide any fractions, the pdfs need to be extended pdfs,
            otherwise fractions need to be provided

    Attributes
    ----------
        pdfs_type      : str
            Is either `ExtendedPDF` if all input pdfs are extended pdfs or `UnextendedPDF` if
            all input pdfs are unextended pdfs
        parameters     : list[utils.parameter.Parameter]
            The pdfs hold parameters, these are extracted and duplicates are removed
        normalisations : list[utils.parameter.Parameter]
            The pdfs can hold normalisations if these are e.g. `ExtendedPDF`, these are extracted and duplicates are removed
        
    Methods
    -------
        evaluate(x : int, float, array-like)
            The summed pdf is evaluated at `x`, `x` is either an int, float or array-like object
        integrate(range : array-like)
            The extended pdf is integrated in the region range
    """
    pdfs: OrderedSet(same_type=True, clones=False)

    def __init__(self, name: str, pdfs: List[PDF], coefficients: List[Union[Parameter, Formula]] = None):

        if coefficients is not None:
            if len(coefficients) == len(pdfs):
                pdfs = RapidExtendedPDFGenerator(pdfs=pdfs, normalisations=coefficients).generate_extended_pdfs()
            elif len(coefficients) != len(pdfs) - 1:
                raise ValidationError(
                    ErrorMessages.get_error_message(7, "fractions", "pdfs - 1"))

        self.pdfs_type = self.check_pdfs(pdfs=pdfs)
        self.check_observables(pdfs=pdfs)
        super().__init__(name=name, observable=pdfs[0].observable)

        self.coefficients = []

        self.pdfs = pdfs
        self.related = {pdf.name: pdf for pdf in pdfs}
        self.set_parameters()
        self.set_coefficients(coefficients)
        
    def get_coefficient(self, name: str) -> Parameter:
        """Returns the coefficient with the given name
        
        Parameters
        ----------
            name : str
                The name of the coefficient
        
        Returns
        -------
            utils.parameter.Parameter
                The coefficient with the given name
        """
        for coefficient in self.coefficients:
            if coefficient.name == name:
                return coefficient
        raise ValueError(f"Could not find coefficient with name {name}")
    
    @property
    def get_coefficient_names(self) -> List[str]:
        """Returns the names of the coefficients
        
        Returns
        -------
            list[str]
                The names of the coefficients
        """
        return [coefficient.name for coefficient in self.coefficients]

    def evaluate(self, x: Union[int, float, np.ndarray]) -> Union[int, float, np.ndarray]:
        """The summed pdf is evaluated at `x`
        
        Parameters
        ----------
            x : int, float, array-like
                `x` represents the values which are evaluated by the pdf
        
        Returns
        ----------
            int, float, array-like
                The data type is the same as `x`, when `x` is array-like, the extended pdf will be evaluated for all elements in `x`
        """
        result = np.float64(0.0)
        if self.pdfs_type == "ExtendedPDF":
            for pdf in self.pdfs:
                result += pdf.evaluate(x)

        elif self.pdfs_type == "UnextendedPDF":
            #apparently a bit slow
            self.__rescale_fractions()
            residual = np.float64(1.0)
            # faster
            for coefficient, pdf in zip(self.coefficients, self.pdfs):
                residual -= coefficient.value
                result += coefficient.value * pdf.evaluate(x)
            #last element has to be evaluated
            result += residual * self.pdfs[-1].evaluate(x)

        return result

    def integrate(self, ranges: Union[list, tuple, np.ndarray]) -> Union[np.float64, np.ndarray]:
        
        if not isinstance(ranges, (list, tuple, np.ndarray)) or len(ranges) < 2:
            raise ValueError(ErrorMessages.get_error_message(15, RANGE_TYPE, "1"))

        integral = np.zeros(len(ranges)-1, dtype=np.float64)
        if self.pdfs_type == "ExtendedPDF":
            for pdf in self.pdfs:
                integral += pdf.integrate(ranges)

        elif self.pdfs_type == "UnextendedPDF":
            #apparently a bit slow
            self.__rescale_fractions()
            residual = np.float64(1.0)
            # faster
            for coefficient, pdf in zip(self.coefficients, self.pdfs):
                residual -= coefficient.value
                integral += coefficient.value * pdf.integrate(ranges)
            #last element has to be evaluated
            integral += residual * self.pdfs[-1].integrate(ranges)

        return integral

    def set_coefficients(self, coefficients: List[Parameter]):
        """This method distinguished between two cases:
        1st case: User inputs unextended pdfs and fractions, then these coefficients are checked if they are valid,
                  if necessary, also they will be rescaled if their sum exceeds 1. 
        2nd case: User inputs extended pdf and no coefficients, so normalisations are stored

        Parameters
        ----------
            coefficients : array-like[utils.parameter.Parameter]
                When the pdfs are unextended, fractions are needed

        """

        if self.pdfs_type == "ExtendedPDF":
            for pdf in self.pdfs:
                self.coefficients.append(pdf.normalisation)

        elif self.pdfs_type == "UnextendedPDF":
            sum_of_fractions = np.sum([fraction.value for fraction in coefficients])
            self.coefficients = coefficients
            self.parameters.extend(coefficients)
            if sum_of_fractions > 1.0:
                self.__rescale_fractions()

        else:
            raise ValueError("The length of the parameters and the number of pdfs have to be compatible!\
                             In case of extended PDFs")

    def __rescale_fractions(self) -> None:
        """If the sum of fractions exceeds 1, then the fractions need to be re?scaled
        because the sum of fractions has to be always 1
        """
        sum_of_fractions = np.sum([fraction.value for fraction in self.coefficients])

        if (1 - sum_of_fractions) < 0:
            for fraction in self.coefficients:
                fraction.value = fraction.value / sum_of_fractions
                if fraction.error is not None:
                    fraction.error = fraction.error / sum_of_fractions
                if fraction.asym_errors is not None:
                    fraction.asym_errors = (fraction.asym_errors[0] / sum_of_fractions,
                                            fraction.asym_errors[1] / sum_of_fractions)

    def set_parameters(self, *args):
        """Every parametrised pdf has parameters which in turn will be also the parameters for the AddPDF instance.
        It automatically takes also the normalisations from the pdf.parameters
        """
        parameters = []
        for pdf in self.pdfs:
            for parameter in pdf.parameters:
                parameters.append(parameter)
                self.related[parameter.name] = parameter

        self.parameters = remove_duplicate_attributes(parameters)

    @staticmethod
    def check_pdfs(pdfs: List[PDF]) -> str:
        """Pdfs have to be the same, either all are extended or all are unextended

        Parameters
        ----------
            pdfs : list[pdfs.basic_pdfs.IFPDF]
                A list of pdfs which should be checked

        Returns
        -------
            str
                Either `ExtendedPDF` in the case all pdfs are extended or `UnextendedPDF` if all are unextended

        Raises
        ------
            ValidationError
                If pdfs are not of the same type in the described sense, then this error will be raised
        """
        first_pdf_type = ("ExtendedPDF" if isinstance(pdfs[0], ExtendedPDF) else "UnextendedPDF")
        is_valid = False
        if first_pdf_type == "ExtendedPDF":
            is_valid = np.all([isinstance(pdf, ExtendedPDF) for pdf in pdfs])
        elif first_pdf_type == "UnextendedPDF":
            is_valid = not np.any([isinstance(pdf, ExtendedPDF) for pdf in pdfs])

        if not is_valid:
            raise ValidationError(
                ErrorMessages.get_error_message(6, "Pdfs", "extended or unextended pdfs"))

        return first_pdf_type

    @staticmethod
    def check_observables(pdfs: Union[list, tuple, np.ndarray]) -> None:
        """Checks whether the observables of the pdfs match or not

        Parameters
        ----------
            pdfs : array-like
                Those are the pdfs which will be summed up in order to form the new pdf

        Raises
        ------
            ValidationError
                If the observable is not the same for all pdfs, then it will raise an error
                since AddPDF cannot have an undefinable observable
        """
        first_observable = pdfs[0].observable
        for index in range(1, len(pdfs)):
            if first_observable != pdfs[index].observable:
                raise ValidationError(
                    ErrorMessages.get_error_message(5, "Observables"))


##########################################################################################################
#                                      Interface Model Implementation                                    #
##########################################################################################################


class IFModel(metaclass=abc.ABCMeta):
    """Model is a superclass which keeps track of the pdfs given to it,
    it is the central object in the fitting procedure, since the fitting
    framework will fetch the parameters and evaluations of these pdfs
    from the model.

    Methods (abstract)
    ------------------
        evaluate(data : int, float, array-like) int, float
            The evaluation of the model is essentially the evaluation of the underlying pdfs with their
            respective normalisations and parameters
    """

    @abstractmethod
    def evaluate(self, data: Sequence[Union[float, list, tuple, np.ndarray]]) -> Union[int, float,
                                                                                np.ndarray, List[np.ndarray]]:
        """Model holds multiple pdfs, each pdf belonging to a distinct category, 
        each pdf in a category is evaluated on its `data`. The result is a list
        of pdf evaluations for a category, the categories are then summed up
        
        Parameters
        ----------
            data : int, float, np.ndarray
                `data` is evaluated on each appropriate category, so the first
                data element is evaluated on the first category, thus on the first pdf

        Returns
        -------
            int, float
                The subcategories are summed up and this sum is returned
        """
        pass

    @abstractmethod
    def integrate(self, bin_container: Sequence[Union[float,List[float]]]) -> Sequence[List[float]]:
        """For doing binned fits, one has to integrate pdfs for given bins, each `data` has multiple datasets,
        the evaluation of the binned log-likelihood has to take place within each category. The result
        is the total log-likelihood for binned `data`

        Parameters
        ----------
            bin_container       : array-like[list[float]]
                `bin_container` has the same length as `data` since for each dataset in `data` a set of bins
                have to be provided, each bin list in the container is such a set of bins

        Returns
        -------
            float
                The total log-likelihood is returned for binned data, this is done by computing
                the log-likelihood for each dataset and then summing up over categories
        """
        pass


##########################################################################################################
#                                         Model Implementation                                           #
##########################################################################################################

class Model(IFModel):
    """Model is a superclass which keeps track of the pdfs given to it,
    it is the central object in the fitting procedure, since the fitting
    framework will fetch the parameters and evaluations of these pdfs
    from the model.

    Parameters
    ----------
        name : str
            The name of the model, e.g. total_shape
        pdfs : basic_pdfs.PDF or list[basic_pdfs.PDF]
            The underlying pdfs which will be managed by the model

    Attributes
    ----------
        normalisations  : list[utils.parameter.Parameter]
            Pdfs may hold normalisations, these are extracted from the respective pdfs and also
            their duplicates are removed
        parameters      : list[utils.parameter.Parameter]
            Pdfs hold parameters, these are extracted from the respective pdfs and also
            their duplicates are removed

    Methods
    -------
        evaluate(data : array-like) int, float
            The evaluation of the model is essentially the evaluation of the underlying pdfs with their
            respective normalisations and parameters
        integrate(bin_container: array-like[list[float]], frequency_container: array-like[List[float]]) float
            The integration method is needed for the binned fits since they need to compute the integral over the
            pdf in well defined ranges which are given by the bins, they are weighted then with the respective frequencies
            per bin
    """
    name = String(min_size=2, max_size=20, lowercase=True)

    def __init__(self, name: str, pdfs: Union[PDF,List[PDF]]):
        self.name = name
        self.pdfs, self.observable = self.__set_pdfs(pdfs)
        self.related = {pdf.name: pdf for pdf in pdfs}
        self.normalisations = self.__set_normalisations(pdfs)
        self.parameters = self.__set_parameters(pdfs)

    def evaluate(self, data: Sequence[Union[float, list, tuple, np.ndarray]]) -> Union[int, float,
                                                                                       np.ndarray, List[np.ndarray]]:
        """Model holds multiple pdfs, each pdf belonging to a distinct category, 
        each pdf in a category is evaluated on its `data`. The result is a list
        of pdf evaluations for a category, the categories are then summed up
        
        Parameters
        ----------
            data : int, float, np.ndarray
                `data` is evaluated on each appropriate category, so the first
                data element is evaluated on the first category, thus on the first pdf

        Returns
        -------
            int, float, ndarray
                The subcategories are summed up and this sum is returned
        """

        result_per_category = []
        if len(self.pdfs) == 1 and isinstance(data[0], float):
            data = np.array([data], dtype = np.float64)
        for pdf, dataset in zip(self.pdfs, data):
            dataset = np.array(dataset, dtype = np.float64)
            result_per_category.append(pdf.evaluate(dataset))

        return result_per_category

    def integrate(self, bin_container: Sequence[Union[float,List[float]]]) -> Sequence[List[float]]:
        """For doing binned fits, one has to integrate pdfs for given bins, each `data` has multiple datasets,
        the evaluation of the binned log-likelihood has to take place within each category. The result
        is the total log-likelihood for binned `data`

        Parameters
        ----------
            bin_container       : array-like[list[float]]
                `bin_container` has the same length as `data` since for each dataset in `data` a set of bins
                have to be provided, each bin list in the container is such a set of bins

        Returns
        -------
            float
                The total log-likelihood is returned for binned data, this is done by computing
                the log-likelihood for each dataset and then summing up over categories
        """
        likelihood_per_category = []
        if len(self.pdfs) == 1 and isinstance(bin_container[0], float):
            bin_container = [bin_container]
        for pdf, bins in zip(self.pdfs, bin_container):
            result_per_bin = pdf.integrate(bins)
            likelihood_per_category.append(result_per_bin)

        likelihood_per_category = np.array(likelihood_per_category, dtype = np.float64)
        return likelihood_per_category

    def __set_pdfs(self, pdfs: Union[PDF, List[PDF]]) -> Tuple[List[PDF], Parameter]:

        __pdfs = []
        if isinstance(pdfs, PDF):
            __pdfs.append(pdfs)
        else:
            self.__check_observables(pdfs)
            return pdfs, pdfs[0].observable

        return __pdfs, __pdfs[0].observable

    def __set_normalisations(self, pdfs: List[PDF]) -> List[Parameter]:
        """The normalisations have to be extracted from the pdfs if they are available, 
        the list will be empty if no underlying pdf holds a normalisation. Also, if normalisations
        are there multiple times, these are removed

        Parameters
        ----------
            pdfs : list[pdfs.basic_pdfs.IFPDF]
                These are the underlying pdfs which may hold normalisations 
        
        Returns
        -------
            list[utils.parameter.Parameter]
                A list of normalisations is returned where all normalisations are unique in it
        """
        normalisations = []
        for pdf in pdfs:
            if isinstance(pdf, AddPDF) and pdf.pdfs_type == "ExtendedPDF":
                for normalisation in pdf.coefficients:
                    normalisations.append(normalisation)
                    self.related[normalisation.name] = normalisation
            elif isinstance(pdf, ExtendedPDF):
                normalisations.append(pdf.normalisation)
                self.related[pdf.normalisation.name] = pdf.normalisation

        return remove_duplicate_attributes(normalisations)

    def __set_parameters(self, pdfs: List[PDF]) -> List[Parameter]:
        """The parameters have to be extracted from the pdfs, if multiple parameters
        are the same, the duplicates will be removed

        Parameters
        ----------
            pdfs : list[pdfs.basic_pdfs.IFPDF]
                These are the underlying pdfs which hold the parameters 
        
        Returns
        -------
            list[utils.parameter.Parameter]
                A list of parameters is returned where all parameters are unique in it
        """
        parameters = []
        for pdf in pdfs:
            for parameter in pdf.parameters:
                parameters.append(parameter)
                self.related[parameter.name] = parameter

        return remove_duplicate_attributes(parameters)

    @staticmethod
    def __check_observables(pdfs: Union[list, tuple, np.ndarray]) -> None:
        """Checks whether the observables of the pdfs match or not
        
        Parameters
        ----------
            pdfs : array-like
                Those are the pdfs which will be summed up in order to form the new pdf

        Raises
        ------
            ValidationError
                If the observable is not the same for all pdfs, then it will raise an error
                since AddPDF cannot have an undefinable observable
        """
        first_observable = pdfs[0].observable
        for index in range(1, len(pdfs)):
            if first_observable != pdfs[index].observable:
                raise ValidationError(
                    ErrorMessages.get_error_message(5, "Observables"))

    def __getitem__(self, item):
        """Model object is subscriptable"""
        return self.related[item]



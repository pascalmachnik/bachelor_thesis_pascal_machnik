import abc
from abc import abstractmethod
import scipy.integrate
from typing import Union, List
import numpy as np

from ..TCFit.parameter import Parameter
from ..TCFit.formula import Formula
from ..TCFit.standards import String
from ..TCFit.error_messages import ErrorMessages

##########################################################################################################
#                                      GLOBAL CONSTANTS                                                  #
##########################################################################################################

RANGE_TYPE = "ndarray, list, tuple or array"

##########################################################################################################
#                                      GLOBAL FUNCTIONS                                                  #
##########################################################################################################


def integrate(func, lower_bound: list, upper_bound: list):
    result = []
    for low, high in zip(lower_bound, upper_bound):
        __result, _ = scipy.integrate.quad(func,low,high)
        result.append(__result)

    return result


#########################################################################################################
#                                      INTERFACE PDF IMPLEMENTATION                                     #
#########################################################################################################


class IFPDF(metaclass=abc.ABCMeta):
    """Interface class for all subsequent pdf classes. This interface should
    guarantee that all pdfs fulfill the same requirements in regard to APIs
    
    Parameters
    ----------
        name       : str
            The `name` of the pdf. The minimum length of the string has to be `2` and not longer than
            `20`.
        observable : parameter.Parameter
            The `observable` of the pdf
    """
    name = String(min_size=2, max_size=20)

    def __init__(self, name: str, observable: Parameter):
        self.name = name
        self.observable = observable
        self.parameters = []
        self.related = {}

    @abstractmethod
    def evaluate(self, x: Union[int, float, np.ndarray]) -> Union[int, float, np.ndarray]:
        """The pdf will be evaluated at `x`

        Parameters
        ----------
            x : int, float, np.ndarray
                `x` has to be a numeric value or array-like so that the pdf can be evaluated
        
        Returns
        -------
            int, float, np.ndarray
                The return type is the same as the input type `x`
        """
        pass

    @abstractmethod
    def integrate(self, ranges: Union[list, tuple, np.ndarray]) -> float:
        """Integrating the pdf in the region `range`
        
        Parameters
        ----------
            ranges : array-like
                `range` must be array-like with size = 2 since the first element represents the lower integration boundary and 
                the second element represents the upper integration boundary
        
        Returns
        -------
            float
                The integration of a pdf yields its area under the curve, thus the area is returned
        """
        pass

    @staticmethod
    @abstractmethod
    def func(x: Union[int, float, np.ndarray], **kwargs) -> Union[int, float, np.ndarray]:
        """This method is essentially the same as evaluate but as a static version of it. This enables the user to compute
        the pdf without the need to instantiate an object of that pdf

        Parameters
        ----------
            x : numeric or array-like
                The values on which the pdf is evaluated

        Returns
        -------
            int, float, np.ndarray
                The return type is the same as the input type `x`
        """
        pass

    @abstractmethod
    def set_parameters(self, parameters: List[Union[Parameter, Formula]]) -> List[
        Parameter]:
        """The user can give as input a parameter or a formula, the pdfs need to store
        only parameters, thus the parameters from formula are extracted

        Parameters
        ----------
            parameters : list[parameter.Parameter or formula.Formula]
                A list of parameters and formulas, if formulas are present, the parameters of the formula
                have to be extracted
        """
        pass

class PDF(IFPDF):

    def __init__(self, name: str, observable: Parameter):
        super().__init__(name, observable)

    def evaluate(self, x: Union[int, float, np.ndarray]) -> Union[int, float, np.ndarray]:
        """ This method has to be implemented since the interface requires it """
        pass

    def integrate(self, ranges: Union[list, tuple, np.ndarray]) -> Union[np.float64, np.ndarray]:
        """The `pdf` is integrated on the `range`

        Parameters
        ----------
            ranges : array-like
                `range` has to be of size `2` since the first element represents the lower integration boundary
                and the upper element represents the upper integration boundary

        Returns
        -------
            float, float
                The integration of a pdf yields its area and its uncertainty under the curve

        Raises
        ------
            ValueError
                Is raised when the length of `range` is not greater than `1`
        """
                # if not isinstance(ranges, (list, tuple, np.ndarray)) or len(ranges) < 2:
                #     raise ValueError(ErrorMessages.get_error_message(15, RANGE_TYPE, "1"))

                # if len(ranges) == 2:
                #     return scipy.integrate.quad(self.evaluate, ranges[0], ranges[1])[0]

                # else:
                #     return integrate(self.evaluate, ranges[:-1], ranges[1:])
    pass

    @staticmethod
    def func(x: Union[int, float, np.ndarray], **kwargs) -> Union[int, float, np.ndarray]:
        """ This method has to be implemented since the interface requires it """
        pass

    def set_parameters(self, parameters: List[Union[Parameter, Formula]]):
        """The user can give as input a parameter or a formula, the pdfs need to store
        only parameters, thus the parameters from formula are extracted

        Parameters
        ----------
            parameters : list[parameter.Parameter or formula.Formula]
                A list of parameters and formulas, if formulas are present, the parameters of the formula
                have to be extracted
        """

        for parameter in parameters:
            if isinstance(parameter, Parameter):
                self.parameters.append(parameter)
                self.related[parameter.name] = parameter
            if isinstance(parameter, Formula):
                self.related[parameter.name] = parameter
                for formula_parameter in parameter.parameters:
                    self.parameters.append(formula_parameter)

    def __getitem__(self, item: str):
        return self.related[item]
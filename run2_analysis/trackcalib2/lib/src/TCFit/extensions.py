from typing import Union, List

import numpy as np
import scipy.integrate

from ..TCFit.basic_pdfs import PDF
from ..TCFit.error_messages import ErrorMessages
from ..TCFit.exceptions import ValidationError
from ..TCFit.formula import Formula
from ..TCFit.parameter import Parameter

RANGE_TYPE = "ndarray, list, tuple or array"
##########################################################################################################
#                                      Extended PDF Implementation                                       #
##########################################################################################################


class ExtendedPDF(PDF):
    """Extended Pdf is a pdf which is extended by a multiplicative normalisation
    
    Parameters
    ----------
        name          : str
            Name of the extended pdf
        pdf           : standard.IFPDF
            A probability density function
        normalisation : list[utils.parameter.Parameter]
            `normalisation` is a list of parameters in the case of formulas, otherwise
            only one parameter is stored in the list

    Attributes
    ----------
        parameters : list[utils.parameter.Parameter]
            The extended pdf inherits the parameters from the underlying pdf instance
        normalisation_formula : formula.Formula
            Is used for the evaluation of the pdf if the normalisation is derived from a formula

    Methods
    -------
        evaluate(x : int, float, array-like)
            The extended pdf is evaluated at `x`, `x` is either an int, float or array-like object
        integrate(range : array-like)
            The extended pdf is integrated in the region range
    """

    def __init__(self, name: str, pdf: PDF, normalisation: Union[Parameter, Formula]):
        super().__init__(name, pdf.observable)
        self.pdf = pdf
        self.normalisation = normalisation
        self.related = self.pdf.related
        self.related[normalisation.name] = normalisation
        self.set_parameters()

    def evaluate(self, x: Union[int, float, np.ndarray]):
        """The extended pdf is evaluated at `x`
        
        Parameters
        ----------
            x : int, float, array-like
                `x` represents the values which are evaluated by the pdf
        
        Returns
        ----------
            int, float, array-like
                The data type is the same as `x`, when `x` is array-like,
                the extended pdf will be evaluated for all elements in `x`
        """

        return self.normalisation.value * self.pdf.evaluate(x)

    def integrate(self, ranges: Union[list, tuple, np.ndarray]):
        """Integrating the extended pdf in the region `range`
        
        Parameters
        ----------
            range : array-like
                `range` must be array-like with size = 2 since the first element represents the lower integration boundary and 
                the second element represents the upper integration boundary
        
        Returns
        -------
            float
                The integration of a pdf yields its area under the curve, thus the area is returned

        Raises
        ------
            ValueError
                Is raised when the length of `range` is not equal to `2` or is not array-like
        """
        # if not isinstance(ranges, (list, tuple, np.ndarray)) or len(ranges) != 2:
        #     raise ValueError(
        #         ErrorMessages.get_error_message(1, "ndarray, list, tuple or array", "2"))

        # return scipy.integrate.quad(self.evaluate, ranges[0], ranges[1])
        if not isinstance(ranges, (list, tuple, np.ndarray)) or len(ranges) < 2:
            raise ValueError(ErrorMessages.get_error_message(15, RANGE_TYPE, "1"))

        return self.normalisation.value * self.pdf.integrate(ranges)

    @staticmethod
    def func(pdf: PDF, normalisation: Union[int, float], x: Union[int, float, np.ndarray], **kwargs) -> Union[
        int, float, np.ndarray]:
        """This method is essentially the same as evaluate but as a static version of it. This enables the user to compute
        the pdf without the need to instantiate an object of that pdf

        Parameters
        ----------
            pdf           : basic_pdfs.IFPDF
                The pdf on which the extended version should be computed on
            normalisation : int, float
                The normalisation which is multiplied with the pdf
            x             : numeric or array-like
                The values on which the pdf is evaluated

        Returns
        -------
            int, float, np.ndarray
                The return type is the same as the input type `x`
        """

        return normalisation * pdf.evaluate(x)

    def set_parameters(self, *args):
        """The user can give as input a parameter or a formula, the pdfs need to store
        only parameters, thus the parameters from formula are extracted

        Parameters
        ----------
            normalisation : utils.parameter.Parameter or formula.Formula
                Normalisation can be either given as a parameter or as a formula, formula consists of multiple
                parameters

        Returns
        -------
            list[utils.parameter.Parameter]
                If the normalisation is a formula, multiple parameters are stored in the list, otherwise
                only one parameter is stored
        """

        for parameter in self.pdf.parameters:
            self.parameters.append(parameter)
        if isinstance(self.normalisation, Parameter):
            self.parameters.append(self.normalisation)
        elif isinstance(self.normalisation, Formula):
            for parameter in self.normalisation.parameters:
                self.parameters.append(parameter)

    def __str__(self):
        return f"{self.normalisation.name} * {self.pdf.name}"

    def __getitem__(self, item: str):
        if item in self.related:
            return self.related[item]
        else:
            print("Item %s not present" % item)


##########################################################################################################
#                                 Rapid Extended PDF Generator Implementation                            #
##########################################################################################################


class RapidExtendedPDFGenerator:
    """Simplifies the instantiation of several extended pdfs

    Parameters
    ----------
        pdfs : list[basic_pdfs.IFPDF]
            These are the pdfs which should be extended
        normalisations : list[utils.parameter.Parameter, formula.Formula]
            The normalisations are used for the construction of extended pdfs since an extended pdf
            is equal to its `normalisation` times the `pdf`

    Raises
    ------
        ValidationError
            When the length of the unextended pdfs does not match the normalisations, then this error will be
            raised
    """

    def __init__(self, pdfs: List[PDF], normalisations: List[Union[
        Parameter, Formula]]):
        if len(pdfs) != len(normalisations):
            raise ValidationError(
                ErrorMessages.get_error_message(1, "PDFs", "normalisations"))
        self.pdfs = pdfs
        self.normalisations = normalisations

    def generate_extended_pdfs(self) -> List[ExtendedPDF]:
        """Extends the given pdfs 

        Returns
        -------
            list[ExtendedPDF]
                The extended pdf instances are returned
        """
        extended_pdfs = []
        for pdf, normalisation in zip(self.pdfs, self.normalisations):
            extended_pdfs.append(ExtendedPDF(name=f"extended_{pdf.name}", pdf=pdf, normalisation=normalisation))

        return extended_pdfs

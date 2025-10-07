import abc
import re
from abc import abstractmethod
from math import sin, cos, tan, sqrt, exp, log, log10
from typing import Sequence, Union, List, Tuple

import numpy as np
import numexpr as ne

from ..TCFit.error_messages import ErrorMessages
from ..TCFit.exceptions import ValidationError
from ..TCFit.parameter import Parameter
from ..TCFit.standards import String
from collections import OrderedDict

from uncertainties import ufloat

##########################################################################################################
#                                    Interface Parser Implementation                                     #
##########################################################################################################


class IFParser(metaclass=abc.ABCMeta):
    """Interface for any parser which needs to evaluate
    string expressions dynamically

    Methods (abstract)
    ------------------
        parse(array-like[int, float]) int, float
            Dynamically evaluates a mathematical string expression, it is also advisible to implement in
            the parse method the grammar rules, such that only well-defined
            string expressions are evaluated
    """

    @abstractmethod
    def parse(self, values: Sequence[Union[int, float]]) -> Union[int, float]:
        """Dynamically evaluates a mathematical string expression after checking if the grammar rules
        apply

        Parameters
        ----------
            values : array-like[int, float]
                `values` contains all the values which are substituted into the mathematical
                expression with the placeholders which are defined in the string expression

        Returns
        -------
            int, float
                When the evaluation of the expression was successful, the result is returned
        """
        pass


##########################################################################################################
#                                     FormulaParser Implementation                                       #
##########################################################################################################


class FormulaParser(IFParser):
    """In order to prevent security issues, expressions
    which are evaluated should be made as safe as possible,
    for this reason this class is introduced which checks what
    input is given to the class and filters out only allowed 
    expressions
    
    Parameters
    ----------
        formula : str
            `formula` is a mathematical string expression, care that formula
            is a very restriced expression which does not accept all subexpressions,
            if you want to see which subexpressions are allowed, please use the method
            `get_information`

    Attributes
    ----------
        safe_dict : dict
            Contains a list of functions which can be used in a formula, these methods
            mostly derive from the package `math`

    Methods
    -------
        parse(values : array-like[int, float]) int, float
            `parse` dynamically evaluates a valid `formula`, the placeholders
            which are marked by @ are substituted with the corresponding values in order

        get_information()
            Prints all the rules which one has to follow in order to build a valid string expression
    """
    __operations = ["+", "-", "*", "/"]
    __base_geometrics = ["sin", "cos", "tan"]
    __base_functions = ["sqrt", "abs", "exp", "log", "log10"]

    def __init__(self, formula: str):
        self.formula = self.__check_formula(formula)
        self.safe_dict = {"sin": sin, "cos": cos, "tan": tan, "sqrt": sqrt, "exp": exp, "log": log, "log10": log10}

    def parse(self, values: Sequence[Union[int, float]]) -> Union[int, float]:
        """Dynamically evaluates a mathematical string expression, no grammatical rules
        are checked

        Parameters
        ----------
            values : array-like[int, float]
                `values` contains all the values which are substituted into the mathematical
                expression with the placeholders which are defined in the string expression

        Returns
        -------
            int, float
                When the evaluation of the expression was successful, the result is returned

        Raises
        ------
            TypeError
                If the length of placeholders is not equal to the number of provided values, 
                the string expression cannot be evaluated and thus an error will be raised
        """
        matches = re.findall("@\d", self.formula)
        if len(matches) != len(values):
            raise TypeError(ErrorMessages.get_error_message(7, "matches", "values"))
        expression = re.sub("@\d", r"{\g<0>}", self.formula).replace("@", "").format(*values)
        code = compile(expression, "<string>", "eval")

        return eval(code, {"__builtins__": {"abs": abs, "nan": np.nan}}, self.safe_dict)

    def get_information(self) -> None:
        info = f"Welcome to the Formula Parser\n \
                 These are the following rules which one has to follow in order to \
                 build a valid mathematical string expression which is called formula:\n \
                 - Operands which one can use: {self.__operations}\n \
                 - Geometric functions which one can use: {self.__base_geometrics}\n \
                 - Base functions which one can use: {self.__base_functions}\n \
                 - Placeholder variables have to be named in the following format: @ followed by \
                   a number between 0 and 9, thus one can only use maximally 10 placeholder variables\n \
                 - You can only use ( and ) brackets, but make sure that all open brackets follow a closed bracket\n \
                 - When using a function, please write it as `func()`, e.g. if you want to use a sin, write \
                   sin(@number)"
        print(info)

    def __check_formula(self, formula: str) -> str:
        """Filters only those subexpressions from the formula
        which are following the rules for formula, the user can still 
        build invalid expressions but these checks will be not done, one 
        only checks for using safe subexpressions which are allowed

        Parameters
        ----------
            formula : str
                `formula` is the mathematical string expression which should be checked

        Raises
        ------
            ValidationError
                If the number of open brackets does not match the number of closed brackets,
                it will raise an error
        """
        formula = re.sub("\s+", "", formula)
        result = ""
        count_brackets = 0
        for position, letter in enumerate(formula):
            if position + 3 <= len(formula) and formula[position:position + 3] in self.__base_geometrics:
                result += f"{formula[position:position + 3]}"
                continue
            for index in range(3, 6):
                if position + index <= len(formula) and formula[position:position + index] in self.__base_functions:
                    result += f"{formula[position:position + index]}"
                    continue
            if letter == "@":
                result += "@"
                continue
            if letter >= "0" and letter <= "9":
                result += f"{letter}"
                continue
            if letter in self.__operations:
                result += f"{letter}"
                continue
            if letter == "(":
                result += f"{letter}"
                count_brackets += 1
                continue
            if letter == ")":
                result += f"{letter}"
                count_brackets -= 1
                continue
        if count_brackets != 0:
            raise ValidationError(ErrorMessages.get_error_message(10,
                                                                  "count_open_bracket",
                                                                  "count_closed_bracket",
                                                                  "for every open bracket there missed a closed bracket"))

        return result

    def get_compiled_formula(self, parameters_names: Tuple[str]) -> str:
        return re.sub(r'@\d+', r'{\g<0>}', self.formula).replace("@", "").format(*parameters_names)


##########################################################################################################
#                                      Formula Implementation                                            #
##########################################################################################################


class SlowFormula:
    """Defines a mathematical expression as a string which can be executed
    dynamically 

    Parameters
    ----------
        name       : str
            The name of the formula, name should be not longer than 20 characters but also
            not shorter than 2 characters, also it is case insensitive
        formula    : str
            A mathematical string expression, for more information on the rules for building
            a formula, please use the method `get_information`
        parameters : list[parameter.Parameter]
            List of physical parameters which will be substituted into the formula

    Attributes
    ----------
        value : int, float
            The result when evaluating the given formula

    Methods
    -------
        get_information()
            Prints all the information/rules on how to build a valid formula
    
    Raises
    ------
        TypeError
            If the parameters are not belonging to the parameter.Parameter class, then
            this error will be raised
    """
    name = String(min_size=2, max_size=20, lowercase=True)
    formula = String()

    def __init__(self, name: str, formula: str, parameters: list):
        self.name = name
        self.formula = formula
        if not np.all([isinstance(parameter, Parameter) for parameter in parameters]):
            raise TypeError(
                ErrorMessages.get_error_message(11, "parameter", "parameter.Parameter"))
        self.parameters = self.set_parameters(parameters)
        self.related = {item.name: item for item in self.parameters}

        self.__parser = FormulaParser(formula=formula)

    def get_information(self) -> None:
        """Prints the information about how to build a valid formula expression"""
        self.__parser.get_information()

    @staticmethod
    def set_parameters(parameters) -> List[Parameter]:
        """Store input parameters into self.parameters. In case of the parameter is not a parameter.Parameter class
        it is checked if the parameter has parameters attribute as for Formula"""
        _parameters = []
        for parameter in parameters:
            if isinstance(parameter, Parameter):
                _parameters.append(parameter)

            else:
                if hasattr(parameter, 'parameters'):
                    _parameters.extend(parameter)

        return _parameters

    @property
    def value(self):
        return self.__parser.parse(values=[parameter.value for parameter in self.parameters])

    def __str__(self):
        return self.formula

    def __getitem__(self, item: str):
        return self.related[item]


class IFFormula(metaclass=abc.ABCMeta):

    def __init__(self, name: str, formula: str, parameters: list):
        self.name = name
        self.formula = formula

        self.parameters = parameters


class Formula(IFFormula):
    """Defines a mathematical expression as a string which can be executed
    dynamically

    Parameters
    ----------
        name       : str
            The name of the formula, name should be not longer than 20 characters but also
            not shorter than 2 characters, also it is case insensitive
        formula    : str
            A mathematical string expression, for more information on the rules for building
            a formula, please use the method `get_information`
        parameters : list[parameter.Parameter]
            List of physical parameters which will be substituted into the formula

    Attributes
    ----------
        value : int, float
            The result when evaluating the given formula

    Methods
    -------
        get_information()
            Prints all the information/rules on how to build a valid formula

    Raises
    ------
        TypeError
            If the parameters are not belonging to the parameter.Parameter class, then
            this error will be raised
    """
#    name = String(min_size=2, max_size=20, lowercase=True)
#    formula = String()

    def __init__(self, name: str, formula: str, parameters: list):
        self.name = name

        self.parameters = self.set_parameters(parameters)
        if not np.all([isinstance(parameter, Parameter) or isinstance(parameter, IFFormula) for parameter in parameters]):
            raise TypeError(
                ErrorMessages.get_error_message(11, "parameter", "Parameter or Formula"))
        self.related = {item.name: item for item in self.parameters}
        self.dictionary = {item.name: item.value for item in self.parameters}

        self.formula = FormulaParser(formula=formula).get_compiled_formula(self.related.keys())

    def get_information(self) -> None:
        """Prints the information about how to build a valid formula expression"""
        self.__parser.get_information()

    @staticmethod
    def set_parameters(parameters) -> List[Parameter]:
        """Store input parameters into self.parameters. In case of the parameter is not a parameter.Parameter class
        it is checked if the parameter has parameters attribute as for Formula"""
        _parameters = []
        for parameter in parameters:
            if isinstance(parameter, Parameter):
                _parameters.append(parameter)

            else:
                if hasattr(parameter, 'parameters'):
                    _parameters.extend(parameter)

        return _parameters

    @property
    def value(self) -> Union[int, float, np.float64]:
        self.dictionary = {item.name: item.value for item in self.parameters}
        return ne.evaluate(self.formula, {}, self.dictionary)
    #"__builtins__": {"abs": np.abs, "nan": np.nan}
    @property
    def error(self) -> Union[int, float, np.float64, ufloat]:
        dictionary = {item.name: ufloat(item.value, item.error) for item in self.parameters}
        return eval(self.formula, dictionary).s

    def __call__(self, *args, **kwargs):
        return self.value()

    def __str__(self):
        return self.formula

    def __getitem__(self, item: str):
        return self.related[item]

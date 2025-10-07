import uuid
from typing import Tuple, Union
import numpy as np

from ..TCFit.error_messages import ErrorMessages
from ..TCFit.standards import Number, Tuples, String, Boolean

""" TODO:
  - Check if asym_errors matches error if both are set
"""


##########################################################################################################
#                                      Parameter Implementation                                          #
##########################################################################################################


class Parameter:
    """Defines a physical parameter and also its properties regarding the behavior in fitting
        
    Parameters
    ----------
        name  : str
            The name of the physical parameter. The `name` should have a minimal string length of `2`
            but not greater than 20, also the `name` will be converted to lowercase letters
        value : int, float, optional
            The value of the physical parameter
        limit : tuple, optional
            Limits the range of the values in `value`. The size of `limit` has to be `2`
        unit  : str, optional
            The unit of the physical parameter
        error : int, float, optional
            The uncertainty associated with the value
        asym_errors : tuple, optional
            1 sigma coverage of the uncertainty. The size of `asym_errors` has to be `2`
        fixed : bool, optional
            If the value should update in the fitting process or not

    Attributes
    ----------
        uuid : str
            The `uuid` uniquely defines the parameter, `name` of the parameter is used to compute
            the uuid, so parameters with the same names cannot exist
    """
#    name = String(min_size=1, max_size=20, lowercase=True)
#    value = Number()
#    limits = Tuples(size=2, only_numbers=True)
#    unit = String()
#    error = Number()
#    asym_errors = Tuples(size=2, only_numbers=True)
#    fixed = Boolean()

    def __init__(self, name: str, value: np.float64 = None,
                 limits: Tuple[np.float64, np.float64] = None,
                 unit: str = None, error: np.float64 = None,
                 asym_errors: Tuple[np.float64, np.float64] = None,
                 fixed: bool = False, latex:str = None):
        self.name = name
        self.value = np.float64(value)
        self.limits = (np.float64(limits[0]), np.float64(limits[1]))
        self.unit = unit
        if latex is not None:
            self.latex = latex
        else:
            if self.unit is not None:
                self.latex = f"{self.name} [{self.unit}]"
            else:
                self.latex = f"{self.name}"
                
        if error is None:
            self.error = np.float64(0.0)
        else:
            self.error = np.float64(error)
        if asym_errors is None and error is None:
            self.asym_errors = (np.float64(0.0), np.float64(0.0))
        elif asym_errors is None and error is not None:
            self.asym_errors = (np.float64(error), np.float64(error))
        elif asym_errors is not None:
            self.asym_errors = (np.float64(asym_errors[0]), np.float64(asym_errors[1]))
        self.fixed = fixed

        self.uuid = uuid.uuid5(uuid.NAMESPACE_DNS, name).hex
        ParameterDatabase.register_parameter(parameter=self)

        # set original values for reset
        self.original_value = self.value
        self.original_limits = self.limits
        self.original_unit = self.unit
        self.original_error = self.error
        self.original_asym_errors = self.asym_errors

    def __str__(self):
        return f'{self.name}: {self.value} \u00B1 {self.error} ({self.asym_errors[0]}, {self.asym_errors[1]})' \
               f' in the range of ({self.limits[0]}, {self.limits[1]})'

    def __call__(self):
        return self.value
    
    def __reset__(self):
        self.value = self.original_value
        self.limits = self.original_limits
        self.unit = self.original_unit
        self.error = self.original_error
        self.asym_errors = self.original_asym_errors
        
    def __vary__(self, scale):
        self.__reset__()
        var_lim = scale * min(self.original_value - self.limits[0], self.limits[1] - self.original_value)
        self.value += np.random.uniform(-var_lim, var_lim)
     
    def fix_value(self, fix_value):
        self.original_value = fix_value
        self.original_error = 0.0
        self.original_asym_errors = 0.0
        self.original_limits = (fix_value, fix_value)
        self.fixed = True
        self.__reset__()

##########################################################################################################
#                                   ParameterDatabase Implementation                                     #
##########################################################################################################
class ParameterDatabase:
    """Local database which manages the parameters, such that
    they are all unique at the same time

    Attributes
    ----------
        parameter_registry : dict
            Holds the parameter `uuid` as keys and the `name` of the parameter
            is stored as a value, this is the central place which defines which parameter
            is unique

    Methods
    -------
        register_parameter(parameter : Parameter)
            Registers a given uuid in the registry if it is unique and duplicates will lead
            to an error, because the name of the parameter is used to generate the uuid
    """
    parameter_registry = {}

    def __init__(self):
        """comment to be added"""
        pass

    @staticmethod
    def register_parameter(parameter: Parameter) -> None:
        """Registers the given `parameter` in the registry in order
        to bookmark unique defined parameters

        Parameters
        ----------
            parameter : Parameter
                The `parameter` contains a uuid and a name, the uuid is generated via the name
                and these information will be stored if they do not already exist

        Raises
        ------
        KeyError
            If the key of the parameter is already in the registry, it will lead to an error 
            since only unique parameters are allowed with unique names
        """
        if parameter.uuid in ParameterDatabase.parameter_registry:
            raise KeyError(
                ErrorMessages.get_error_message(8, f"{parameter.uuid}", "parameter_registry",
                                                               ", please do not duplicate the name of a parameter!"))
        ParameterDatabase.parameter_registry[parameter.uuid] = parameter.name

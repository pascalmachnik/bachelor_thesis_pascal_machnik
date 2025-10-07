import abc
from abc import abstractmethod
from typing import Tuple, List, Sequence, Union

import numpy as np

from ..TCFit.standards import String
from ..TCFit.error_messages import ErrorMessages


def check_data_compatibility(limits: Union[list, tuple, np.ndarray], data: Union[int, float, list, np.ndarray]):
    """Data has to be in the range of the given limit, all the data values
    which are not inside the limit boundaries, are getting removed

    Parameters
    ----------
        limits : array-like
            `limit` is an array-like object of size `2` which constrains the domain of acceptable values
        data  : numeric or array-like
            `data` represents the values which the pdf evaluates, data which lies outside of `limit` is excluded

    Raises
    ------
        ValueError
            This error is raised if the data is a numeric value and is not inside the `limit`
        TypeError
            This error is raised if the data is neither array-like or numeric
    """
    if isinstance(data, (int, float)):
        if limits[0] <= data <= limits[1]:
            return data
        raise ValueError(ErrorMessages.get_error_message(2))

    elif isinstance(data, (list, tuple, np.ndarray)):
        if type(data) == list:
            data = np.array(data)
        return data[(data >= limits[0]) & (data <= limits[1])]

    raise TypeError(
        ErrorMessages.get_error_message(3, "np.float64, float, int, np.ndarray, list"))


class IFDataset(metaclass=abc.ABCMeta):
    name = String(min_size=2, max_size=20, lowercase=True)

    def __init__(self, name: str, binned: bool = False):
        self.binned = binned
        self.name = name

    @abstractmethod
    def __add__(self, other):
        pass

    @abstractmethod
    def __set_data(self, *args, **kwargs):
        pass


class BinnedDataset(IFDataset):

    def __init__(self, name: str, data: Union[float, List[float], np.array, np.ndarray], errors: Union[float, List[float], np.array, np.ndarray]):
        super(BinnedDataset, self).__init__(name=name, binned=True)
        self.data == self.__set_data(data)

    def __set_data(self):
        pass

    def __add__(self, other):
        pass

    def __str__(self):
        pass


class UnbinnedDataset(IFDataset):

    def __init__(self, name: str,
                 data: Union[float, List[float], np.array],
                 weights: Union[float, List[float], np.array] = None,
                 column_names: Union[str,List[str]] = None,
                 category: str = None):
        super(UnbinnedDataset, self).__init__(name=name, binned=False)

        self.data, self.weights, self.column_names = self.__set_data(data, weights, column_names)


    def __set_data(self, data, weights, names):
        """Set the unbinned data structure. This reminds a ntuples with column names and optional weights columns

        Parameters
        ----------
        data    : float, List[float], np.array
        weights : weights associated to each entries (same lenght of ntuples)
        names   : column names
        """

        __data = np.array(data)
        """Check if dataset and weights have the same size"""
        if isinstance(names, str):
            names = [names]

        __names = dict(enumerate(names))
        __columns = {v: k for k, v in __names.items()}

        if data.ndim != len(weights):
            raise ValueError("Length of dataset not compatible with length of weights")

        __weights = np.array(weights)

        if names is not None:
            if data.ndim != len(__names):
                raise ValueError("Length of columns name ")


        return __data, __weights, __columns


    def __get_binned_set(self,**kwargs):
        pass

    def __getitem__(self,item):
        __items = []
        __columns = []
        if isinstance(item, str):
            __items = [item]
        elif isinstance(item, list):
            __items = item

        for __item in __items:
            if __item in self.column_names:
                __columns.append(self.column_names[__item])
            else:
                raise ValueError(f"{item} is not a column")


        return self.data[:,__columns]

    def __add__(self, other):
        pass

    def __str__(self):
        pass

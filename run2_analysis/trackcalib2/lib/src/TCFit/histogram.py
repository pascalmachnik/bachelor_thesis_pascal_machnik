import numpy as np
from typing import Union, Tuple, List, Sequence
from stats import eval_unc


def evaluate_uncertainties(values: list[Union[int,float]]) -> List[List[float]]:
    """Evaluate the uncertainties for a list of float"""

    uncertainties = []
    for element in values:
        if isinstance(element, int):
            uncertainties.append(list(eval_unc(element)))
        if isinstance(element, float):
            uncertainties.append([[element ** 0.5], [element ** 0.5]])
    return uncertainties


def check_shape(list1: Union[list, np.ndarray], list2: Union[list, np.ndarray]):
    """check if 2 list have the same shape"""

    if np.array(list1).shape == np.array(list2).shape:
        return True
    return False


class Histogram:

    """1D and 2D histogram class"""
    def __init__(self, name: str, **kwargs):

        self.name = name
        self.contents = None
        self.uncertainties = None
        self.n_bins = []
        self.ranges = None
        self.boundaries = None
        self.w2 = True

        if "w2" in kwargs:
            self.w2 = kwargs["w2"]

        if "n_bins" in kwargs and "ranges" in kwargs:
            __n_bins = kwargs["n_bins"]
            __ranges = kwargs["ranges"]
            self.n_bins, self.boundaries = self.__set_dimensions(n_bins=__n_bins, ranges=__ranges)
            self.contents = np.zeros(shape=self.n_bins)

        elif "contents" in kwargs and "boundaries" in kwargs:
            __uncertainties = None
            if "uncertainties" in kwargs:
                __uncertainties = kwargs['uncertainties']

            self.contents, self.uncertainties = self.__set_content(contents=kwargs['contents'],
                                                                   uncertainties=__uncertainties)

            self.n_bis, self.boundaries = self.__set_dimensions(contents=kwargs['contents'],
                                                                boundaries=kwargs['boundaries'])

    @staticmethod
    def __set_content(contents: Sequence[float], uncertainties: Sequence[float] = None) -> Tuple[Sequence[float],
                                                                                                 Sequence[float]]:
        __uncertainties = []
        if uncertainties == None:
            if isinstance(contents[0], float) or isinstance(contents[1], int):
                __uncertainties.extend(evaluate_uncertainties(contents))
            elif isinstance(contents[0], list):
                for c in contents:
                    __uncertainties.append(evaluate_uncertainties(c))

            return contents, __uncertainties

        else:
            if check_shape(contents, uncertainties):
                return contents, __uncertainties

            else:
                raise ValueError("Contents and Uncertainties lists don't have the same shape")

    @staticmethod
    def __set_dimensions(**kwargs):

        n_bins = []
        boundaries = []
        if "n_bins" in kwargs and "ranges" in kwargs:
            if len(kwargs['n_bins']) == len(kwargs['ranges']):
                for n, r in zip(kwargs['n_bins'], kwargs['ranges']):
                    boundaries.append(np.linspace(r[0], r[1], n + 1))
                    n_bins.append(n)
                return tuple(n_bins), boundaries
            else:
                raise TypeError("Number of ranges and bins should be the same")

        elif "contents" in kwargs and "boundaries" in kwargs:
            # list of list of contents and boundaries
            __contents = kwargs['contents']
            __boundaries = kwargs['boundaries']
            if len(__contents) == len(__boundaries):
                for c, b in zip(__contents, __boundaries):
                    if isinstance(c, list) and isinstance(b, list) and len(c) == len(b) - 1:
                        n_bins.append(len(c))
                        boundaries.append(b)
                    else:
                        raise TypeError("Size of contents not compatible with size of boundaries")

            elif len(__contents) == len(__boundaries) - 1:
                if isinstance(__contents, list) and isinstance(__boundaries, list):
                    n_bins.append(len(__contents))
                    boundaries.extend(__boundaries)
                else:
                    raise TypeError("Size of contents not compatible with size of boundaries")

            else:
                raise TypeError("Size of contents not compatible with size of boundaries")

        else:
            raise TypeError("Issue with initialisation of histogram")


    def scale(self):
        pass

    def mean(self):
        pass

    def median(self):
        pass

    def rms(self):
        pass
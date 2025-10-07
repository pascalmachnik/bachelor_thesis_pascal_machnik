import abc
from abc import abstractmethod

import numpy as np


class Validator(metaclass=abc.ABCMeta):
    def __set_name__(self, owner, name):
        self.name = name

    def __set__(self, instance, value):
        if value == None:
            return None
        self.validate(value)
        processed_value = self.process(value, self.options)
        if processed_value != None:
            value = processed_value
        instance.__dict__[self.name] = value

    def __get__(self, instance, owner):
        if self.name in instance.__dict__:
            return instance.__dict__[self.name]

    def __delete__(self, instance):
        del instance.__dict__[self.name]

    @abstractmethod
    def validate(self, value):
        pass

    @abstractmethod
    def process(self, value, **options):
        pass


class String(Validator):
    def __init__(self, min_size=None, max_size=None, **options):
        self.min_size = min_size
        self.max_size = max_size
        self.options = options

    def validate(self, value):
        if not isinstance(value, str):
            raise TypeError(f'Expected {value!r} to be an str')
        if self.min_size is not None and len(value) < self.min_size:
            raise ValueError(f'Expected {value!r} to be no smaller than {self.min_size!r}')
        if self.max_size is not None and len(value) > self.max_size:
            raise ValueError(f'Expected {value!r} to be no bigger than {self.max_size!r}')

    def process(self, value, options):
        for option_name, option_value in options.items():
            if option_name == "lowercase" and option_value is True:
                return value.lower()


class Number(Validator):
    def __init__(self, limit=None, **options):
        self.limit = limit
        self.options = options

    def validate(self, value):
        if type(value) not in [np.float64, float, int]:
            raise TypeError(f'Expected {value!r} to be a number of type: np.float64, float or int')
        if self.limit is not None and value not in range(self.limit[0], self.limit[1]):
            raise ValueError(f'Expected {self.name!r} to be in the range of {self.limit!r}')

    def process(self, value, options):
        """ This method has to be implemented since the interface requires it """
        pass

    def __get__(self, instance, owner):
        if self.name in instance.__dict__:
            return np.float64(instance.__dict__[self.name])
        return np.float64(0.0)

    
class Tuples(Validator):
    def __init__(self, size=None, only_numbers=None, **options):
        self.size = size
        self.only_numbers = only_numbers
        self.options = options

    def validate(self, value):
        if value is None:
            return
        if not isinstance(value, tuple):
            raise TypeError(f'Expected {value!r} to be a tuple')
        if self.size is not None and len(value) > self.size:
            raise ValueError(f'Expected {self.name!r} to be a tuple of size {self.size!r}')
        if self.only_numbers is not None:
            is_only_numbers = all(type(item) in [np.float64, float, int] for item in value)
            if is_only_numbers is not self.only_numbers:
                if self.only_numbers is True:
                    raise ValueError(f'Expected {self.name!r} to contain only numbers')


    def process(self, value, options):
        """ This method has to be implemented since the interface requires it """
        pass

    def __get__(self, instance, owner):
        if self.name in instance.__dict__:
            return instance.__dict__[self.name]
        return 0, 0

    
class Boolean(Validator):
    def __init__(self, **options):
        self.options = options

    def validate(self, value):
        if not isinstance(value, bool):
            raise TypeError(f'Expected {value!r} to be a boolean')

    def process(self, value, options):
        """ This method has to be implemented since the interface requires it """
        pass


class OrderedSet(Validator):
    def __init__(self, same_type=False, clones=False, **options):
        self.same_type = same_type
        self.clones = clones
        self.options = options

    def validate(self, value):
        if value == None:
            return
        if not isinstance(value, list) or not isinstance(value, tuple):
            raise TypeError(f"Expected {value!r} to be a list or tuple")
        if self.same_type == True:
            type_to_validate = type(value[0])
            is_same_type = np.all(type(item) == type_to_validate for item in value)
            if is_same_type == False:
                raise ValueError(f"Expected {self.name!r} contains clones")

    def process(self, value, options):
        """ This method has to be implemented since the interface requires it """
        pass

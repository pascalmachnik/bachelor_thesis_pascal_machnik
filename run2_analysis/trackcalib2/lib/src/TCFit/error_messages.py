base_messages = {
    1:  'The length of the corresponding data type: {0} should be not greater than {1}!',
    2:  'Data value has to be in the limit of the observable!',
    3:  'Data does not belong to the following types: {0}!',
    4:  'The optional argument limit needs to be a tuple!',
    5:  '{0} have to be the same!',
    6:  '{0} have to be of the same same type, namely {1}',
    7:  'The length of the array-like object {0} should be equal to the length of {1}!',
    8:  'The key {0} is already in the dictionary {1}{2}!',
    9:  'The array-like object {0} has to be of dimension {1}!',
    10: '{0} is not equal to {1} since {2}!',
    11: 'Instance {0} does not belong to class {1}!',
    12: 'Key {0} is missing in dictionary {1}!',
    13: 'Key {0} value has to be integer or an array of boundaries, None is passed!',
    14: 'Problem with argument {0}',
    15: 'The length of the array-like object {0} should be greater than {1}!',
    16: 'Bin ranges cannot be guessed from binned dataset',
}


class ErrorMessages:
    """ErrorMessages should unify the way error messages are thrown at the user, so 
    that these do not have to be redefined, although they remain the same. This class is 
    a static class only, thus one should not initialise this class but just use the 
    public methods.

    Methods
    -------
        @staticmethod
        get_error_message(identifier : int, *args)
            `identifier` uniquely determines an error message, some error messages can be customized by providing
            `*args`. The `*args` will be substituted into the error messages.
    """

    @staticmethod
    def get_error_message(identifier: int, *args):
        """Returns the requested error message

        Parameters
        ----------
            identifier : int
                `identifier` uniquely determines an error message
            *args : str
                `*args` should be strings which will then be substituted into the error message if the requested
                error message is customizable
        
        Returns
        -------
            str
                The requested error message is returned as a string
        """

        return base_messages[identifier].format(*args)
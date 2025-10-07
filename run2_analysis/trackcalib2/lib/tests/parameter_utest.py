import unittest

from src.TCFit.parameter import Parameter


class TestVariables(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        pass


    def test_parameter_s01(self):
        """SCENARIO
        Test the initialisation of the parameter (part 1)
        Full parameter list is given, check if it matches the string representation
        """
        
        """ PREPARATION """

        """ EXECUTION """
        parameter   = Parameter(name ="sigma", value = 0.25, limit = (0.1, 0.34), unit ="GeV", error = 0.05, asym_errors = (0.02, 0.07), fixed = True)
        test_output = str(parameter)

        """ VERIFICATION """
        actual_output = "sigma: 0.25 \u00B1 0.05 (0.02, 0.07) in the range of (0.1, 0.34)"

        self.assertEqual(first = actual_output, second = test_output)


    def test_parameter_s02(self):
        """SCENARIO
        Test the initialisation of the parameter (part 2)
        Part of the parameter list is given, check if it matches the string representation, here asym_error
        is left out, so asym_error should be (error, error)
        """

        """ PREPARATION """
        parameter = Parameter(name ="mu", value = 3200, limit = (2800, 3200), unit ="GeV", error = 200, fixed = False)

        """ EXECUTION """
        test_output = str(parameter)

        """ VERIFICATION """
        actual_output = "mu: 3200 \u00B1 200 (200, 200) in the range of (2800, 3200)"

        self.assertEqual(first = actual_output, second = test_output)


    def test_parameter_e01(self):
        """SCENARIO
        Test the faulty initialisation of the parameter when some input 
        is not according to the requirements, here the name is too long
        """

        """ PREPARATION """

        """ EXECUTION """

        """ VERIFICATION """
        with self.assertRaises(ValueError):
            parameter = Parameter(name ="This is a super good parameter name", value = 3200,
                                                      limit = (2800, 3200), unit = "GeV", error = 200, fixed = False)


    def test_parameter_e02(self):
        """SCENARIO
        Test the faulty initialisation of the parameter when some input 
        is not according to the requirements, here the tuple is not of range 2 but 1
        """

        """ PREPARATION """

        """ EXECUTION """

        """ VERIFICATION """
        with self.assertRaises(TypeError):
            parameter = Parameter(name ="mu", value = 3200, limit = (2800), unit ="GeV", error = 200, fixed = False)


    def test_parameter_e03(self):
        """SCENARIO
        Test the faulty initialisation of the parameter when some input
        has the wrong type, here the value has the wrong type
        """

        """ PREPARATION """

        """ EXECUTION """

        """ VERIFICATION """
        with self.assertRaises(TypeError):
            parameter = Parameter(name ="mu", value ="3200", limit = (2800, 3400), unit ="GeV", error = 200, fixed = False)


    def test_if_parameters_are_unique_e01(self):
        """SCENARIO
        Test whether parameters with two names raise an error because
        they are not unique
        """

        """ PREPARATION """

        """ EXECUTION """
        first_parameter  = Parameter(name ="gaussian mu")

        """ VERIFICATION """
        with self.assertRaises(KeyError):
            second_parameter = Parameter(name = "gaussian mu")


    def test_parameter_string_representation_s01(self):
        """SCENARIO
        Testing the string representation of parameters
        """

        """ PREPARATION """
        parameter = utils.parameter.Parameter(name = "str rep", value = 20)

        """ EXECUTION """
        print(parameter)

        """ VERIFICATION """
        self.assertTrue(True)

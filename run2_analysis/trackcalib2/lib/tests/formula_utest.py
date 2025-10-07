import unittest

import src.TCFit.formula
import src.TCFit.parameter


class TestFormula(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        pass

    
    def test_formula_initialisation_s01(self):
        """SCENARIO
        Checks whether the initialisation of formula works correctly
        """

        """ PREPARATION """
        parameter_one = src.TCFit.parameter.Parameter(name ="par sample", value = 2, limit = (0, 4), error = .5, fixed = False)
        parameter_two = src.TCFit.parameter.Parameter(name ="par sample 2", value = 4, limit = (0, 10), error = 1, fixed = False)

        """ EXECUTION """
        formula        = src.TCFit.formula.Formula(name ="form 1", formula ="@0 * @1", parameters = [parameter_one, parameter_two])
        act_name       = formula.name
        act_formula    = formula.formula
        act_parameters = formula.parameters
        act_value      = formula.value

        """ VERIFICATION """
        exp_name       = "form 1"
        exp_formula    = "@0 * @1"
        exp_parameters = [parameter_one, parameter_two]
        exp_value      = parameter_one.value * parameter_two.value

        self.assertEqual(first = act_name, second = exp_name)
        self.assertEqual(first = act_formula, second = exp_formula)
        self.assertEqual(first = act_parameters, second = exp_parameters)
        self.assertEqual(first = act_value, second = exp_value)


    def test_formula_evaluation_s01(self):
        """SCENARIO
        Check if the evaluation of the formula works correctly
        """

        """ PREPARATION """
        parameter_one = src.TCFit.parameter.Parameter(name ="p_one", value = 10, limit = (5, 15), unit ="m", error = 1, fixed = False)
        parameter_two = src.TCFit.parameter.Parameter(name ="p_two", value = 5, limit = (-5, 15), unit ="m", error = .5, fixed = True)
        formula       = src.TCFit.formula.Formula(name ="test", formula ="@0 + @1", parameters = [parameter_one, parameter_two])

        """ EXECUTION """
        act_result = formula.value

        """ VERIFICATION """
        exp_result = parameter_one.value + parameter_two.value

        self.assertEqual(first = act_result, second = exp_result)


    def test_formula_evaluation_s02(self):
        """SCENARIO
        Check the evaluation of formula for a realistic scenario, here we
        will use a typical example for defining a formula for fitting
        """

        """ PREPARATION """
        sample_size       = 1000
        signal_yield      = src.TCFit.parameter.Parameter("signal_yield", value =0.5 * sample_size, limit = (0, sample_size), error =0.05 * sample_size, fixed = False)
        efficiency_sig    = src.TCFit.parameter.Parameter("efficiency_sig", value = 0.5, limit = (0, 1), error = 0.05, fixed = False)
        signal_yield_fail = src.TCFit.formula.Formula("f_fail", "@0*(1.-@1)", [signal_yield, efficiency_sig])

        """ EXECUTION """
        act_result = signal_yield_fail.value

        """ VERIFICATION """
        exp_result = 0.5 * sample_size * (1 - 0.5)
        
        self.assertEqual(first = act_result, second = exp_result)


    def test_formula_updating_s01(self):
        """SCENARIO
        When a parameter value changes, the formula should also adjust accordingly
        when one tries to get the value of the formula
        """

        """ PREPARATION """
        parameter_one = src.TCFit.parameter.Parameter(name ="one", value = 10, limit = (5, 15), unit ="m", error = 1, fixed = False)
        parameter_two = src.TCFit.parameter.Parameter(name ="two", value = 5, limit = (-5, 15), unit ="m", error = .5, fixed = True)
        formula       = src.TCFit.formula.Formula(name ="test formula", formula ="@0 + @1", parameters = [parameter_one, parameter_two])

        """ EXECUTION """
        parameter_one.value = 13
        act_result          = formula.value

        """ VERIFICATION """
        exp_result = parameter_one.value + parameter_two.value

        self.assertEqual(first = act_result, second = exp_result)

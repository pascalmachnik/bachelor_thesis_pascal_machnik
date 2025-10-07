import unittest

import numpy as np

from src.TCFit.basic_pdfs import Exponential, Gaussian, CrystalBall, DoubleSidedCrystalBall
from src.TCFit.formula import Formula
from src.TCFit.parameter import Parameter


class TestStandardPdfs(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        """These parameters are needed for the subsequent tests, the pdfs need these parameters"""
        self.observable    = Parameter(name ="observable", limits = (2800, 3200))
        self.mu                 = Parameter(name ="mu", value = 3000)
        self.sigma              = Parameter(name ="sigma", value = 10)

        self.slope              = Parameter(name ="slope", value = 2)

        self.alpha              = Parameter(name ="alpha", value = 2)
        self.n                  = Parameter(name ="n_cryst", value = 1)

        self.alpha_l            = Parameter(name ="alpha_l", value = 2)
        self.alpha_r            = Parameter(name ="alpha_r", value = 1)
        self.n_l                = Parameter(name ="n_l", value = 1)
        self.n_r                = Parameter(name ="n_r", value = 1)

        self.data = np.array([0, 10, 20, 50, 100, 205.24, 525.2, 1000, 2000, 2800, 3400.52])

    
    def test_pdf_initialization_s01(self):
        """SCENARIO
        All the pdfs which are given as pdfs are initialised with correct parameters in order to check if
        some error occurs, if no errors are thrown, this test passes
        """

        """ PREPARATION """

        """ EXECUTION """
        gaus      = Gaussian(name ="gaus", observable = self.observable, mu = self.mu, sigma = self.sigma)
        exp       = Exponential(name ="exponential", observable = self.observable, slope = self.slope)
        crystal   = CrystalBall(name ="crystalball", observable = self.observable, mu = self.mu,
                                                     sigma = self.sigma, alpha = self.alpha, n = self.n)
        dscrystal = DoubleSidedCrystalBall(name ="dscrystalball", observable = self.observable, mu = self.mu,
                                                                sigma = self.sigma, alpha_l = self.alpha_l, n_l = self.n_l,
                                                                alpha_r = self.alpha_r, n_r = self.n_r)

        """ VERIFICATION """
        self.assertTrue(True)


    def test_is_valid_pdf_s01(self):
        """SCENARIO
        Check if the pdfs integrate to 1 when integrating over the whole range
        The following pdfs are checked:
            - Gaussian
            - Exponential
            - Crystalball
            - Double-Sided Crystalball (not yet correctly normalised)
        """

        """ PREPARATION """
        gaus      = Gaussian(name ="gaus", observable = self.observable, mu = self.mu, sigma = self.sigma)
        exp       = Exponential(name ="exponential", observable = self.observable, slope = self.slope)
        crystal   = CrystalBall(name ="crystalball", observable = self.observable, mu = self.mu,
                                                     sigma = self.sigma, alpha = self.alpha, n = self.n)
        """ dscrystal = pdfs.standard.DoubleSidedCrystalBall(name = "dscrystalball", observable = self.observable, mu = self.mu,
                                                         sigma = self.sigma, alpha_l = self.alpha_l, n_l = self.n_l,
                                                         alpha_r = self.alpha_r, n_r = self.n_r) """

        """ EXECUTION """
        act_gaus_area, _      = gaus.integrate(range = self.observable.limits)
        act_exp_area, _       = exp.integrate(range = self.observable.limits)
        act_crystal_area, _   = crystal.integrate(range = self.observable.limits)
        # act_dscrystal_area, _ = dscrystal.integrate(range = self.observable.limit)

        """ VERIFICATION """
        self.assertAlmostEqual(first = act_gaus_area, second = 1.00, delta = 0.01)
        self.assertAlmostEqual(first = act_exp_area, second = 1.00, delta = 0.01)
        self.assertAlmostEqual(first = act_crystal_area, second = 1.00, delta = 0.01)
        # self.assertAlmostEqual(first = act_dscrystal_area, second = 1.00, delta = 0.01)


    def test_parameter_list_s01(self):
        """SCENARIO
        The parameter list should contain the parameters of the pdfs, we want to check
        here if these are the same objects
        The pdfs on which we test are:
            - Gaussian
            - Crystalball
        """

        """ PREPARATION """
        gaus    = Gaussian(name ="gaus", observable = self.observable, mu = self.mu, sigma = self.sigma)
        crystal = CrystalBall(name ="crystalball", observable = self.observable, mu = self.mu,
                                                   sigma = self.sigma, alpha = self.alpha, n = self.n)

        """ EXECUTION """
        act_gaus_parameters = gaus.parameters
        act_crys_parameters = crystal.parameters

        """ VERIFICATION """
        self.assertEqual(first = act_gaus_parameters, second = [self.mu, self.sigma])
        self.assertEqual(first = act_crys_parameters, second = [self.mu, self.sigma, self.alpha, self.n])


    def test_evaluation_of_pdfs_s01(self):
        """SCENARIO
        All the pdfs are evaluated on a data array in order to check if these pdfs are correctly computed
        """

        """ PREPARATION """
        gaus      = Gaussian(name ="gaus", observable = self.observable, mu = self.mu, sigma = self.sigma)
        exp       = Exponential(name ="exponential", observable = self.observable, slope = self.slope)
        crys      = CrystalBall(name ="crystalball", observable = self.observable, mu = self.mu,
                                                     sigma = self.sigma, alpha = self.alpha, n = self.n)
        dscrystal = DoubleSidedCrystalBall(name ="dscrystalball", observable = self.observable, mu = self.mu,
                                                                sigma = self.sigma, alpha_l = self.alpha_l, n_l = self.n_l,
                                                                alpha_r = self.alpha_r, n_r = self.n_r)

        """ EXECUTION """
        """ actual_output_gaus  = gaus.evaluate(self.data)
        actual_output_exp   = exp.evaluate(self.data)
        actual_output_crys  = crys.evaluate(self.data)
        actual_output_dcrys = dscrystal.evaluate(self.data) """

        """ VERIFICATION """
        """ test_output_gaus   = scipy.stats.norm.pdf(2800, self.mu.value, self.sigma.value)
        test_output_exp    = scipy.stats.expon.pdf(self.data[:5], 0, 1/self.slope.value)
        
        z = (np.array(self.data[4:8]) - self.mu.value) / self.sigma.value
        test_output_crys  = scipy.stats.crystalball.pdf(z, self.n.value, self.alpha.value)

        self.assertAlmostEqual(first = actual_output_gaus[0], second = test_output_gaus, places = 7)
        for act_item, test_item in zip(actual_output_exp, test_output_exp):
            self.assertAlmostEqual(first = act_item, second = test_item, places = 7)
        for act_item, test_item in zip(actual_output_crys, test_output_crys):
            self.assertAlmostEqual(first = act_item, second = test_item, places = 7) """


    def test_pdf_evaluation_with_formula_s01(self):
        """SCENARIO
        Pdfs need also to work with formulas, in this test we use
        a simple Gaussian in order to test if it works and also if
        the parameters list is correctly set up
        """

        """ PREPARATION """
        factor_1 = Parameter(name ="mu factor 1", value = 30)
        factor_2 = Parameter(name ="mu factor 2", value = 100)
        mu       = Formula(name ="mu formula", formula ="@0 * @1", parameters = [factor_1, factor_2])

        gaussian = Gaussian(name ="gaus", observable = self.observable, mu = mu, sigma = self.sigma)

        """ EXECUTION """
        act_result     = gaussian.evaluate(self.data)
        act_parameters = gaussian.parameters

        """ VERIFICATION """
        exp_result     = Gaussian.func(self.data, mu = 3000, sigma = 10, limits = self.observable.limits)
        exp_parameters = [factor_1, factor_2, self.sigma]

        self.assertEqual(first = act_result, second = exp_result)
        self.assertEqual(first = act_parameters, second = exp_parameters)


    def test_pdf_evaluation_with_formula_s02(self):
        """SCENARIO
        We test the same scenario as in s01 but with an exponential pdf
        """

        """ PREPARATION """
        factor_1 = Parameter(name ="slope factor 1", value = 30)
        factor_2 = Parameter(name ="slope factor 2", value = 100)
        slope    = Formula(name ="slope formula", formula ="@0 / @1", parameters = [factor_1, factor_2])

        exponential = Exponential(name ="exp", observable = self.observable, slope = slope)

        """ EXECUTION """
        act_result     = exponential.evaluate(self.data)
        act_parameters = exponential.parameters

        """ VERIFICATION """
        exp_result     = Exponential.func(self.data, slope = 0.3, limits= self.observable.limits)
        exp_parameters = [factor_1, factor_2]

        for index in range(len(act_result)):
            self.assertEqual(first = act_result[index], second = exp_result[index])
        self.assertEqual(first = act_parameters, second = exp_parameters)
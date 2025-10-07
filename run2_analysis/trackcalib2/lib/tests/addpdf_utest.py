import unittest

import numpy as np

import src.TCFit.basic_pdfs
import src.TCFit.extensions
import src.TCFit.model
import src.TCFit.parameter


class TestAddPDF(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.gaus_observable    = src.TCFit.parameter.Parameter(name ="gaus_obs", limit = (2800, 3200))
        self.mu                 = src.TCFit.parameter.Parameter(name ="mu", value = 3000)
        self.sigma              = src.TCFit.parameter.Parameter(name ="sigma", value = 10)

        self.exp_observable     = src.TCFit.parameter.Parameter(name ="exp_obs", limit = (0, 100))
        self.slope              = src.TCFit.parameter.Parameter(name ="slope", value = 2)

        self.gaussian    = src.TCFit.basic_pdfs.Gaussian(name ="gaussian", observable = self.gaus_observable, mu = self.mu, sigma = self.sigma)
        self.exponential = src.TCFit.basic_pdfs.Exponential(name ="exponential", observable = self.gaus_observable, slope = self.slope)

        self.gaus_norm   = src.TCFit.parameter.Parameter(name ="Gaussian Normalisation", value = 100, limit = (80, 120), error = 5, fixed = False)
        self.exp_norm    = src.TCFit.parameter.Parameter(name ="Exp Normalisation", value = 500, limit = (200, 800), error = 20, fixed = False)

        self.ext_gaussian    = src.TCFit.extensions.ExtendedPDF(name ="extended gaussian", pdf = self.gaussian, normalisation = self.gaus_norm)
        self.ext_exponential = src.TCFit.extensions.ExtendedPDF(name ="extended exponential", pdf = self.exponential, normalisation = self.exp_norm)

        self.fraction_gaussian    = src.TCFit.parameter.Parameter(name ="gaussian frac", value = 0.4, limit = (0.3, 0.5), error = 0.02, fixed = True)

        self.data = np.array([0, 10, 20, 50, 100, 205.24, 525.2, 1000, 2000, 2800, 3400.52])


    def test_initialisation_s01(self):
        """SCENARIO
        We test the initialisation of AddPDF, the check includes
        whether the unextended pdfs are the same, the fractions and also the parameters
        """

        """ PREPARATION """

        """ EXECUTION """
        added_pdf        = src.TCFit.model.AddPDF(name ="gaus + exp", pdfs = [self.gaussian, self.exponential], fractions = [self.fraction_gaussian])
        act_pdfs         = added_pdf.pdfs
        act_parameters   = added_pdf.parameters
        act_fractions    = added_pdf.fractions

        """ VERIFICATION """
        # check whether pdfs are the same
        self.assertEqual(first = act_pdfs, second = [self.gaussian, self.exponential])
        # check whether parameters are the same
        self.assertEqual(first = act_parameters, second = [self.mu, self.sigma, self.slope])
        # check whether coefficients are as expected
        exp_normalisations = [self.fraction_gaussian]
        self.assertEqual(first = act_fractions, second = exp_normalisations)


    def test_initialisation_s02(self):
        """SCENARIO
        We test the initialisation of AddPDF regarding extended pdfs,
        no fractions are provided, but the comparison includes the pdfs
        and parameters and whether the fractions is None
        """

        """ PREPARATION """

        """ EXECUTION """
        added_pdf      = src.TCFit.model.AddPDF(name ="extended gaus + exp", pdfs = [self.ext_gaussian, self.ext_exponential], fractions = None)
        act_pdfs       = added_pdf.pdfs
        act_parameters = added_pdf.parameters
        act_fractions  = added_pdf.fractions

        """ VERIFICATION """
        self.assertEqual(first = act_pdfs, second = [self.ext_gaussian, self.ext_exponential])
        self.assertEqual(first = act_parameters, second = [self.mu, self.sigma, self.slope])
        self.assertEqual(first = act_fractions, second = None)


    def test_evaluation_s01(self):
        """SCENARIO
        We test if the AddPDF instance correctly evaluates the sum of pdfs for the unextended case
        """

        """ PREPARATION """
        added_pdf = src.TCFit.model.AddPDF(name ="gaus + exp", pdfs = [self.gaussian, self.exponential], fractions = [self.fraction_gaussian])

        """ EXECUTION """
        act_result = added_pdf.evaluate(self.data)

        """ VERIFICATION """
        exp_result = self.fraction_gaussian.value * self.gaussian.evaluate(self.data) + (1 - self.fraction_gaussian.value) * self.exponential.evaluate(self.data)

        self.assertEqual(first = act_result, second = exp_result)


    def test_evaluation_s02(self):
        """SCENARIO
        We test if the AddPDF instance correctly evaluates the sum of pdfs for the extended case
        """

        """ PREPARATION """
        added_pdf = src.TCFit.model.AddPDF(name ="extended gaus + exp", pdfs = [self.ext_gaussian, self.ext_exponential], fractions = None)

        """ EXECUTION """
        act_result = added_pdf.evaluate(self.data)

        """ VERIFICATION """
        exp_result = self.ext_gaussian.evaluate(self.data) + self.ext_exponential.evaluate(self.data)

        self.assertEqual(first = act_result, second = exp_result)
import unittest

import numpy as np

import src.TCFit.basic_pdfs
import src.TCFit.extensions
import src.TCFit.parameter


class TestExtendedPDF(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        """These parameters are needed for the subsequent tests, the pdfs need these parameters"""
        self.gaus_observable    = src.TCFit.parameter.Parameter(name ="gaus_obs", limit = (2800, 3200))
        self.mu                 = src.TCFit.parameter.Parameter(name ="mu", value = 3000)
        self.sigma              = src.TCFit.parameter.Parameter(name ="sigma", value = 10)

        self.exp_observable     = src.TCFit.parameter.Parameter(name ="exp_obs", limit = (0, 100))
        self.slope              = src.TCFit.parameter.Parameter(name ="slope", value = 2)

        self.gaussian           = src.TCFit.basic_pdfs.Gaussian(name ="gaussian", observable = self.gaus_observable, mu = self.mu, sigma = self.sigma)
        self.exponential        = src.TCFit.basic_pdfs.Exponential(name ="exponential", observable = self.exp_observable, slope = self.slope)

        self.gaus_norm          = src.TCFit.parameter.Parameter(name ="Gaussian Normalisation", value = 100, limit = (80, 120), error = 5, fixed = False)
        self.exp_norm           = src.TCFit.parameter.Parameter(name ="Exp Normalisation", value = 500, limit = (200, 800), error = 20, fixed = False)

        self.data = np.array([0, 10, 20, 50, 100, 205.24, 525.2, 1000, 2000, 2800, 3400.52])

    
    def test_extpdf_initialization_s01(self):
        """SCENARIO
        To test is the initialisation of the extended pdf if it works as expected. 
        For that a Gaussian pdf is used which should be extended, then compare the 
        pdf and normalisation that the extended pdf is holding with the input if these are the same
        """

        """ PREPARATION """

        """ EXECUTION """
        extended_gaussian = src.TCFit.extensions.ExtendedPDF(name ="extended gaussian", pdf = self.gaussian, normalisation = self.gaus_norm)
        act_pdf           = extended_gaussian.pdf
        act_norm          = extended_gaussian.normalisation

        """ VERIFICATION """
        self.assertEqual(first = act_pdf, second = self.gaussian)
        self.assertEqual(first = act_norm, second = self.gaus_norm)


    def test_extpdf_evaluation_s01(self):
        """SCENARIO 
        To test is the evaluation of the extended pdfs. It is evaluated on a simple data sample.
        """

        """ PREPARATION """
        extended_gaussian = src.TCFit.extensions.ExtendedPDF(name ="extended exponential", pdf = self.exponential, normalisation = self.exp_norm)

        """ EXECUTION """
        act_result = extended_gaussian.evaluate(self.data)

        """ VERIFICATION """
        exp_result = self.exp_norm.value * self.exponential.evaluate(self.data)

        self.assertTrue(expr = np.all([act_item == exp_item for act_item, exp_item in zip(act_result, exp_result)]))


    def test_rapid_extpdf_s01(self):
        """SCENARIO
        To test is the initialisation of multiple pdfs with the RapidExtendedPDFGenerator. As input is given
        two pdfs and two normalisations and one should correctly get two extended pdfs. The comparison is done
        for the name, pdf and normalisation.
        """

        """ PREPARATION """

        """ EXECUTION """
        extended_pdf_gen  = src.TCFit.extensions.RapidExtendedPDFGenerator(pdfs = [self.gaussian, self.exponential], normalisations = [self.gaus_norm, self.exp_norm])
        act_extended_pdfs = extended_pdf_gen.generate_extended_pdfs()

        """ VERIFICATION """
        exp_extended_gaussian    = src.TCFit.extensions.ExtendedPDF(name ="extended gaussian", pdf = self.gaussian, normalisation = self.gaus_norm)
        exp_extended_exponential = src.TCFit.extensions.ExtendedPDF(name ="extended exponential", pdf = self.exponential, normalisation = self.exp_norm)
        exp_extended_pdfs        = [exp_extended_gaussian, exp_extended_exponential]

        for act_extended_pdf, exp_extended_pdf in zip(act_extended_pdfs, exp_extended_pdfs):
            self.assertEqual(first = act_extended_pdf.name, second = exp_extended_pdf.name)
            self.assertEqual(first = act_extended_pdf.pdf, second = exp_extended_pdf.pdf)
            self.assertEqual(first = act_extended_pdf.normalisation, second = exp_extended_pdf.normalisation)
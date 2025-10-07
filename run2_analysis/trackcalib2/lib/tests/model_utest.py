import unittest

import numpy as np

import src.TCFit.basic_pdfs
import src.TCFit.extensions
import src.TCFit.model
import src.TCFit.parameter


class TestModel(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.observable   = src.TCFit.parameter.Parameter(name ="gaus_obs", limit = (2800, 3200))
        self.mu           = src.TCFit.parameter.Parameter(name ="mu", value = 3000)
        self.sigma_1      = src.TCFit.parameter.Parameter(name ="sigma", value = 10)
        self.sigma_2      = src.TCFit.parameter.Parameter(name ="sigma2", value = 5)
        self.slope        = src.TCFit.parameter.Parameter(name ="slope", value = 2)

        self.shared_norm  = src.TCFit.parameter.Parameter(name ="Gaussian Norm", value = 100, limit = (80, 120), error = 5, fixed = False)

        self.gaussian_1   = src.TCFit.basic_pdfs.Gaussian(name ="gaussian", observable = self.observable, mu = self.mu, sigma = self.sigma_1)
        self.gaussian_2   = src.TCFit.basic_pdfs.Gaussian(name ="gaussian 2", observable = self.observable, mu = self.mu, sigma = self.sigma_2)
        self.exponential  = src.TCFit.basic_pdfs.Exponential(name ="exponential", observable = self.observable, slope = self.slope)

        self.ext_gaussian    = src.TCFit.extensions.ExtendedPDF(name ="extended gaussian", pdf = self.gaussian_1, normalisation = self.shared_norm)
        self.ext_exponential = src.TCFit.extensions.ExtendedPDF(name ="extended exponential", pdf = self.exponential, normalisation = self.shared_norm)

        self.data = np.array([0, 10, 20, 50, 100, 205.24, 525.2, 1000, 2000, 2800, 3400.52])


    def test_model_initialisation_s01(self):
        """SCENARIO
        The model needs as input a name and pdfs,
        the model should make sure to extract correctly
        the pdf parameters and normalisations
        """

        """ PREPARATION """

        """ EXECUTION """
        model = src.TCFit.model.Model(name ="unique model", pdfs = [self.gaussian_1, self.ext_exponential])
        act_parameters     = model.parameters
        act_normalisations = model.normalisations

        """ VERIFICATION """
        exp_parameters     = [self.mu, self.sigma_1, self.slope]
        exp_normalisations = [self.shared_norm]
        self.assertEqual(first = act_parameters, second = exp_parameters)
        self.assertEqual(first = act_normalisations, second = exp_normalisations)

    
    def test_model_parameter_filtering_s01(self):
        """SCENARIO
        When the pdfs share parameters with each other, multiple
        parameters are equal to each other, thus they need to be 
        removed from the parameters list since they are not linearly
        independent
        Here we use two gaussians which share the parameter mu with each other
        """

        """ PREPARATION """
        model = src.TCFit.model.Model(name ="Another model", pdfs = [self.gaussian_1, self.gaussian_2])

        """ EXECUTION """
        act_parameters = model.parameters

        """ VERIFICATION """
        exp_parameters = [self.mu, self.sigma_1, self.sigma_2]
        self.assertEqual(first = act_parameters, second = exp_parameters)


    def test_model_normalisation_filtering_s01(self):
        """SCENARIO
        When the pdfs share the same normalisation with each other, multiple
        normalisations are equal, they need to be removed
        Here we test this on extended pdfs, one extended exponential pdf
        and one extended gaussian pdf
        """

        """ PREPARATION """
        model = src.TCFit.model.Model(name ="Extended Model", pdfs = [self.ext_exponential, self.ext_gaussian])

        """ EXECUTION """
        act_normalisations = model.normalisations

        """ VERIFICATION """
        exp_normalisations = [self.shared_norm]
        self.assertEqual(first = act_normalisations, second = exp_normalisations)


    def test_model_evaluation_s01(self):
        """SCENARIO
        The model implements the evaluate method, check if this method
        works properly. In this test, a gaussian and an exponential pdf
        are the reference pdfs for the model.
        """

        """ PREPARATION """
        model = src.TCFit.model.Model(name ="super model", pdfs = [self.gaussian_1, self.exponential])

        """ EXECUTION """
        act_result = model.evaluate(self.data)

        """ VERIFICATION """
        exp_result = self.gaussian_1.evaluate(self.data) + self.exponential.evaluate(self.data)
        self.assertEqual(first = act_result, second = exp_result)
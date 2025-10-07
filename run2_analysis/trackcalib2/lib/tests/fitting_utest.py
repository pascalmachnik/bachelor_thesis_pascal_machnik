import unittest

import iminuit
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats

import src.TCFit.basic_pdfs
import src.TCFit.extensions
import src.TCFit.formula
import src.TCFit.fit
import src.TCFit.model
import src.TCFit.parameter


class TestLikeLihoodFit(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        """ Define the parameters """
        self.observable    = src.TCFit.parameter.Parameter(name ="observable", limits = (100, 3000))
        self.observable_2  = src.TCFit.parameter.Parameter(name ="observable_2", limits = (0, 3000))

        # normalisations
        self.norm_exp      = src.TCFit.parameter.Parameter(name ="exponential norm", value = 2, limits = (0, 5), unit ="", error = 1, fixed = False)
        self.norm_cb       = src.TCFit.parameter.Parameter(name ="cb norm", value = 1800, limits = (1500, 2500), unit ="", error = 150, fixed = False)

        # for the gaussian
        self.mu            = src.TCFit.parameter.Parameter(name ="mu", value = 20, limits = (10, 50), error = 4, fixed = False)
        self.sigma         = src.TCFit.parameter.Parameter(name ="sigma", value = 5, limits = (0.2, 6.4), error = 0.1, fixed = False)

        # for the exponential
        self.slope         = src.TCFit.parameter.Parameter(name ="slope", value = 0.2, limits = (0.1, 0.5), error = 0.05, fixed = False)
        self.slope_2       = src.TCFit.parameter.Parameter(name ="slope 2", value = 10e-3, limits = (10e-5, 10e-2), unit ="", error = 10e-5, fixed = False)

        # for the crystalball
        self.mu_cb         = src.TCFit.parameter.Parameter(name ="cb location", value = 1000, limits = (900, 2000), unit ="", error = 5, fixed = False)
        self.sigma_cb      = src.TCFit.parameter.Parameter(name ="cb spread", value = 150, limits = (100, 300), unit ="", error = 5, fixed = False)
        self.alpha_cb      = src.TCFit.parameter.Parameter(name ="cb alpha", value = 1, limits = (0, 3), unit ="", error = 0, fixed = True)
        self.n_cb          = src.TCFit.parameter.Parameter(name ="cb n", value = 2, limits = (0, 4), unit ="", error = 0, fixed = True)

        """ Define the pdfs """
        self.gaussian      = src.TCFit.basic_pdfs.Gaussian(name ="gaussian", observable = self.observable, mu = self.mu, sigma = self.sigma)
        self.exponential   = src.TCFit.basic_pdfs.Exponential(name ="exponential", observable = self.observable, slope = self.slope)

        self.crystalball   = src.TCFit.basic_pdfs.CrystalBall(name ="crystalball", observable = self.observable_2, mu = self.mu, sigma = self.sigma,
                                                              alpha = self.alpha_cb, n = self.n_cb)
        self.exponential_2 = src.TCFit.basic_pdfs.Exponential(name ="exponential", observable = self.observable_2, slope = self.slope_2)

        """ Define the extended pdfs """
        self.ext_cb        = src.TCFit.extensions.ExtendedPDF(name ="extended cb", pdf = self.crystalball, normalisation = self.norm_cb)
        self.ext_exp       = src.TCFit.extensions.ExtendedPDF(name ="extended exp", pdf = self.exponential_2, normalisation = self.norm_exp)

        """ Define the addpdfs """
        self.sum_pdf       = src.TCFit.model.AddPDF(name ="ext cb + ext exp", pdfs = [self.ext_cb, self.ext_exp])

        """ Define the model """
        self.gaus_exp_model = src.TCFit.model.Model(name ="gaus + exp", pdfs = [self.gaussian, self.exponential])
        self.cb_exp_model   = src.TCFit.model.Model(name ="Model", pdfs = [self.sum_pdf])

    
    def test_fitting_s01(self):
        """SCENARIO
        The LikelihoodFit should be able to fit different scenarios
        of pdfs. The first scenario investigated is if it can fit
        a crystalball pdf + an exponential pdf
        TODO: wrong - to be fixed
        """

        """ PREPARATION """
        X        = np.arange(0, 3000, 1)
        mu_cb    = 1800
        sigma_cb = 200
        norm_cb  = 2000
        beta_cb  = 1
        m_cb     = 2
        Y_cb     = norm_cb * scipy.stats.crystalball.pdf(X, beta_cb, m_cb, mu_cb, sigma_cb)

        norm_exp = 3
        lamb     = 10e-4
        Y_exp    = norm_exp * np.exp(-lamb * X)

        data     = Y_cb + Y_exp

        """ EXECUTION """

        likelihood = src.TCFit.fit.LikelihoodFit(data = [data], model = self.cb_exp_model)
        minuit     = likelihood.minimise(hesse = False, minos = True)

        """ VERIFICATION """
        self.assertIsInstance(minuit, iminuit.Minuit)


    def test_fitting_s02(self):
        """SCENARIO
        A simple fit should be made with a single extended
        gaussian in this case, the result is compared with the made 
        plot
        """

        """ PREPARATION """
        observable = src.TCFit.parameter.Parameter(name ="another obs", limits = (0, 3000))
        mu         = src.TCFit.parameter.Parameter(name ="another mu", value = 1200, limits = (1000, 1800), error = 5, fixed = False)
        sigma      = src.TCFit.parameter.Parameter(name ="another sigma", value = 120, limits = (100, 400), error = 1, fixed = False)
        norm       = src.TCFit.parameter.Parameter(name ="another norm", value = 800, limits = (700, 1400), error = 20, fixed = False)

        gaus       = src.TCFit.basic_pdfs.Gaussian(name ="another gaus", observable = observable, mu = mu, sigma = sigma)
        ext_gaus   = src.TCFit.extensions.ExtendedPDF(name ="another ext gaus", pdf = gaus, normalisation = norm)

        model      = src.TCFit.model.Model(name ="another model", pdfs = [ext_gaus])

        x = np.linspace(0, 3000)
        Y = 1000 * scipy.stats.norm.pdf(x, loc = 1500, scale = 200)
        X = np.random.normal(loc = 1500, scale = 200, size = 1000)

        """ EXECUTION """
        fitting = src.TCFit.fit.LikelihoodFit(data = [X], model = model)
        minuit  = fitting.minimise(hesse = True, minos = True)

        """ VERIFICATION """
        _, bins, _ = plt.hist(X)
        bin_width  = np.diff(bins)[0]
        plt.scatter(x, Y)
        plt.plot(x, bin_width * src.TCFit.extensions.ExtendedPDF.func(pdf = gaus, normalisation = norm.value, x = x))
        plt.show()


    def test_fitting_s03(self):
        """SCENARIO
        This test checks whether the framework is able to fit a simultaneous fit in 
        a simple case. The case here is a double Gaussian which shares the parameter
        mu but have different normalisation and different sigma
        """

        """ PREPARATION """
        observable = src.TCFit.parameter.Parameter(name ="another obs 2", limits = (0, 3000))
        mu         = src.TCFit.parameter.Parameter(name ="another mu 2", value = 1200, limits = (1000, 1800), error = 5, fixed = False)
        sigma      = src.TCFit.parameter.Parameter(name ="another sigma 2", value = 120, limits = (100, 400), error = 1, fixed = False)
        sigma_2    = src.TCFit.parameter.Parameter(name ="another sigma 2.2", value = 350, limits = (200, 500), error = 2, fixed = False)
        norm       = src.TCFit.parameter.Parameter(name ="another norm 2", value = 800, limits = (700, 1400), error = 20, fixed = False)
        norm_2     = src.TCFit.parameter.Parameter(name ="another norm 2.2", value = 1100, limits = (1000, 1400), error = 30, fixed = False)

        gaussian   = src.TCFit.basic_pdfs.Gaussian(name ="gaus", observable = observable, mu = mu, sigma = sigma)
        gaussian_2 = src.TCFit.basic_pdfs.Gaussian(name ="gaus 2", observable = observable, mu = mu, sigma = sigma_2)

        ext_gaus   = src.TCFit.extensions.ExtendedPDF(name ="ext gaus", pdf = gaussian, normalisation = norm)
        ext_gaus_2 = src.TCFit.extensions.ExtendedPDF(name ="ext gaus 2", pdf = gaussian_2, normalisation = norm_2)

        model      = src.TCFit.model.Model(name ="first model", pdfs = [ext_gaus, ext_gaus_2])

        X_1        = np.random.normal(loc = 1350, scale = 250, size = 1000)
        X_2        = np.random.normal(loc = 1350, scale = 400, size = 1350)

        """ EXECUTION """
        fitting    = src.TCFit.fit.LikelihoodFit(data = [X_1, X_2], model = model, monitor_on = True)
        minuit     = fitting.minimise(hesse = True, minos = True)

        """ VERIFICATION """
        x = np.linspace(0, 3000)
        fig, ax       = plt.subplots(nrows = 1, ncols = 2, figsize = (12, 6))
        _, bins, _    = ax[0].hist(X_1, label = "Hist 1")
        _, bins_2, _  = ax[1].hist(X_2, label = "Hist 2") 
        bin_width     = np.diff(bins)[0]
        bin_width_2   = np.diff(bins_2)[0]
        ax[0].plot(x, bin_width * src.TCFit.extensions.ExtendedPDF.func(pdf = gaussian, normalisation = norm.value, x = x), label ="Fit 1")
        ax[1].plot(x, bin_width_2 * src.TCFit.extensions.ExtendedPDF.func(pdf = gaussian_2, normalisation = norm_2.value, x = x), label ="Fit 2")
        ax[0].legend()
        ax[1].legend()
        plt.show()


    def test_fitting_s04(self):
        """SCENARIO
        Until now, we tested only if the unbinned fit worked, now
        we want to make sure that also the binned fit is working properly.
        For that, we test again on two gaussians as in s03 but with binning on, but
        without providing any bins explicitly
        """

        """ PREPARATION """
        observable = src.TCFit.parameter.Parameter(name ="obs s04", limits = (0, 3000))
        mu         = src.TCFit.parameter.Parameter(name ="mu s04", value = 1200, limits = (1000, 1800), error = 5, fixed = False)
        sigma      = src.TCFit.parameter.Parameter(name ="sigma 2 s04", value = 120, limits = (100, 400), error = 1, fixed = False)
        sigma_2    = src.TCFit.parameter.Parameter(name ="sigma 2.2 s04", value = 350, limits = (200, 500), error = 2, fixed = False)
        norm       = src.TCFit.parameter.Parameter(name ="norm 2 s04", value = 800, limits = (700, 1400), error = 20, fixed = False)
        norm_2     = src.TCFit.parameter.Parameter(name ="norm 2.2 s04", value = 1100, limits = (1000, 1400), error = 30, fixed = False)

        gaussian   = src.TCFit.basic_pdfs.Gaussian(name ="gaus", observable = observable, mu = mu, sigma = sigma)
        gaussian_2 = src.TCFit.basic_pdfs.Gaussian(name ="gaus 2", observable = observable, mu = mu, sigma = sigma_2)

        ext_gaus   = src.TCFit.extensions.ExtendedPDF(name ="ext gaus", pdf = gaussian, normalisation = norm)
        ext_gaus_2 = src.TCFit.extensions.ExtendedPDF(name ="ext gaus 2", pdf = gaussian_2, normalisation = norm_2)

        model      = src.TCFit.model.Model(name ="first model", pdfs = [ext_gaus, ext_gaus_2])

        X_1        = np.random.normal(loc = 1350, scale = 250, size = 1000)
        X_2        = np.random.normal(loc = 1350, scale = 400, size = 1350)

        """ EXECUTION """
        fitting    = src.TCFit.fit.LikelihoodFit(data = [X_1, X_2], model = model, monitor_on = True, binning = True, bins = None)
        minuit     = fitting.minimise(hesse = True, minos = True)

        """ VERIFICATION """
        x = np.linspace(0, 3000)
        fig, ax       = plt.subplots(nrows = 1, ncols = 2, figsize = (12, 6))
        _, bins, _    = ax[0].hist(X_1, label = "Hist 1")
        _, bins_2, _  = ax[1].hist(X_2, label = "Hist 2") 
        bin_width     = np.diff(bins)[0]
        bin_width_2   = np.diff(bins_2)[0]
        ax[0].plot(x, bin_width * src.TCFit.extensions.ExtendedPDF.func(pdf = gaussian, normalisation = norm.value, x = x), label ="Fit 1")
        ax[1].plot(x, bin_width_2 * src.TCFit.extensions.ExtendedPDF.func(pdf = gaussian_2, normalisation = norm_2.value, x = x), label ="Fit 2")
        ax[0].legend()
        ax[1].legend()
        plt.show()


    def test_fitting_s05(self):
        """SCENARIO
        We reuse the scenario defined in s04 but with formula in the norm and sigma
        """

        """ PREPARATION """
        observable = src.TCFit.parameter.Parameter(name ="obs s05", limits = (0, 3000))
        mu         = src.TCFit.parameter.Parameter(name ="mu s05", value = 1200, limits = (1000, 1800), error = 5, fixed = False)

        factor_1   = src.TCFit.parameter.Parameter(name ="factor 1 sigma s05", value = 200, limits = (100, 400), error = 5, fixed = False)
        factor_2   = src.TCFit.parameter.Parameter(name ="factor 2 sigma s05", value = 0.4, limits = (0, 1), error = 0.01, fixed = False)
        sigma      = src.TCFit.formula.Formula(name ="sigma 2 s05", formula ="@0 * (1 - @1)", parameters = [factor_1, factor_2])
        sigma_2    = src.TCFit.parameter.Parameter(name ="sigma 2.2 s05", value = 350, limits = (200, 500), error = 2, fixed = False)

        signal_eff = src.TCFit.parameter.Parameter(name ="signal eff s05", value = 200, limits = (100, 350), error = 5, fixed = False)
        back_eff   = src.TCFit.parameter.Parameter(name ="back eff s05", value = 600, limits = (500, 750), error = 5, fixed = False)
        norm       = src.TCFit.formula.Formula(name ="norm 2 s05", formula ="@0 + @1", parameters = [signal_eff, back_eff])
        norm_2     = src.TCFit.parameter.Parameter(name ="norm 2.2 s05", value = 1100, limits = (1000, 1400), error = 30, fixed = False)

        gaussian   = src.TCFit.basic_pdfs.Gaussian(name ="gaus", observable = observable, mu = mu, sigma = sigma)
        gaussian_2 = src.TCFit.basic_pdfs.Gaussian(name ="gaus 2", observable = observable, mu = mu, sigma = sigma_2)

        ext_gaus   = src.TCFit.extensions.ExtendedPDF(name ="ext gaus", pdf = gaussian, normalisation = norm)
        ext_gaus_2 = src.TCFit.extensions.ExtendedPDF(name ="ext gaus 2", pdf = gaussian_2, normalisation = norm_2)

        model      = src.TCFit.model.Model(name ="first model", pdfs = [ext_gaus, ext_gaus_2])

        X_1        = np.random.normal(loc = 1350, scale = 250, size = 1000)
        X_2        = np.random.normal(loc = 1350, scale = 400, size = 1350)

        """ EXECUTION """
        fitting    = src.TCFit.fit.LikelihoodFit(data = [X_1, X_2], model = model, monitor_on = True, binning = True, bins = None)
        minuit     = fitting.minimise(hesse = True, minos = True)

        """ VERIFICATION """
        x = np.linspace(0, 3000)
        fig, ax       = plt.subplots(nrows = 1, ncols = 2, figsize = (12, 6))
        _, bins, _    = ax[0].hist(X_1, label = "Hist 1")
        _, bins_2, _  = ax[1].hist(X_2, label = "Hist 2") 
        bin_width     = np.diff(bins)[0]
        bin_width_2   = np.diff(bins_2)[0]
        ax[0].plot(x, bin_width * src.TCFit.extensions.ExtendedPDF.func(pdf = gaussian, normalisation = norm.value, x = x), label ="Fit 1")
        ax[1].plot(x, bin_width_2 * src.TCFit.extensions.ExtendedPDF.func(pdf = gaussian_2, normalisation = norm_2.value, x = x), label ="Fit 2")
        ax[0].legend()
        ax[1].legend()
        plt.show()
# -*- coding: utf-8 -*-
from typing import Optional, List, Tuple
import warnings
import numpy as np
from lmfit import Model, Parameters


class CurveFit:
    """
    Base class for fitting curves, subclasses have specific tools, plot and attributes
    """
    default_fit_kws = {
        'method': 'lstsq',
        'max_nfev': 5000
    }
    def_limit = (-np.inf, np.inf)

    # limits = defaultdict(default_factory=def_limit)
    def __init__(self, x, y, name=None, function=None, ic=None, limits=None, ):
        """
        :param x: 1D array
        :param y: 1D array
        :param name: str
        :param model: lmfit.models.Model (or None)
        """
        self.x = x
        self.y = y
        self.name = name
        self.model = Model(function)
        self.function = function
        self.params = Parameters()
        self.fitted_params = None

        if ic is None:
            self.ic = {}
        else:
            self.ic = ic

        if limits is None:
            self.limits = {}
        else:
            self.limits = limits

        self.result = None

    def guess_ic(self,guess=None):
        pass  # do nothing, can't guess. Subclass override

    def guess_limit(self):
        pass  # do nothing, can't guess. Subclass override

    def make_params(self):
        params = Parameters()
        for param, value in self.ic.items():
            limit = self.limits.get(param, self.def_limit)
            params.add(param, value=float(value), min=limit[0], max=limit[1],
                       vary=not isinstance(value, str)) # no vary if value is string
        self.params = params

    def fit_model(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            # Ignore UserWarning: Using UFloat objects with std_dev==0 may give unexpected results.
            result = self.model.fit(self.y, self.params, x=self.x, **self.default_fit_kws)
        self.result = result
        fitted_params = self.get_fit_params()
        return fitted_params, result

    def get_fit_params(self):
        fitted_params = {}
        for key, param in self.result.params.items():
            fitted_params[key] = param.value
        self.fitted_params = fitted_params
        return self.fitted_params

    def run_default(self):
        if not self.ic:
            self.guess_ic()
        if not self.ic:  # if above cannot make IC ...
            raise Exception('ic is required, no parameters to fit')
        if not self.limits:
            self.guess_limit()
        self.make_params()
        fitted_params, result = self.fit_model()
        return fitted_params, result

    def best_guess(self, guess_list=None, mode='r2'):
        if guess_list is None:
            return self.run_default()

        residuals = []
        fits = []
        results = []

        for guess in guess_list:
            self.guess_ic(guess)
            self.guess_limit()
            self.make_params()
            fitted_params, result = self.run_default()

            fits.append(fitted_params)
            results.append(result)
            if mode == 'r':  # closest area
                residuals.append(np.sum(result.residual))
            elif mode == 'r2':  # smallest error
                residuals.append(np.sum(result.residual * result.residual))
            ## something based on number of args?
        min_ind = np.argmin(residuals)
        self.fitted_params = fits[min_ind]
        self.result = results[min_ind]
        return self.fitted_params, self.result

    def __call__(self, x):
        if self.result is None:
            return None  # default param call?
        else:
            return self.result.eval(x=x)



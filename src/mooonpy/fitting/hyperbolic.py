# -*- coding: utf-8 -*-
from mooonpy.fitting.fitting import CurveFit
import numpy as np
from mooonpy.tools.math_utils import hyperbola, gaussian_turn, gaussian_quad_turn


class HyperbolaFit(CurveFit):
    def __init__(self, x, y, name=None, function=hyperbola, ic=None, limits=None):
        super().__init__(x, y, name, function, ic, limits)

    def guess_ic(self,guess=None):
        if self.function is hyperbola:
            slope_guess = (self.y[-1] - self.y[0]) / (self.x[-1] - self.x[0])
            if 'X0' not in self.ic:
                self.ic['X0'] = np.mean(self.x)
            if 'Y0' not in self.ic:
                self.ic['Y0'] = np.mean(self.y)
            if 'a' not in self.ic:
                self.ic['a'] = slope_guess
            if 'b' not in self.ic:
                self.ic['b'] = slope_guess
            if 'c' not in self.ic:
                self.ic['c'] = abs(self.x[-1] - self.x[0]) / 10
        else:
            raise Exception('Model not guessable')

    def guess_limit(self):
        slope_guess = (self.y[-1] - self.y[0]) / (self.x[-1] - self.x[0])
        if 'X0' not in self.limits:
            self.limits['X0'] = (min(self.x), max(self.x))
        if 'Y0' not in self.limits:
            self.limits['Y0'] = (min(self.y), max(self.y))
        if 'a' not in self.limits:
            self.limits['a'] = (-10 * abs(slope_guess), 10 * abs(slope_guess))
        if 'b' not in self.limits:
            self.limits['b'] = (-10 * abs(slope_guess), 10 * abs(slope_guess))
        if 'c' not in self.limits:
            if 'c' in self.ic:
                ic = float(self.ic['c'])
            else:
                ic = abs(self.x[-1] - self.x[0]) / 10
            self.limits['c'] = (0.1 * ic, 10 * ic)


class GaussianTurnFit(CurveFit):
    def __init__(self, x, y, name=None, function=gaussian_turn, ic=None, limits=None):
        super().__init__(x, y, name, function, ic, limits)

    def guess_ic(self,guess=None):
        slope_guess = (self.y[-1] - self.y[0]) / (self.x[-1] - self.x[0])
        if 'X0' not in self.ic:
            self.ic['X0'] = np.mean(self.x)
        if 'Y0' not in self.ic:
            self.ic['Y0'] = np.mean(self.y)
        if 'a' not in self.ic:
            self.ic['a'] = slope_guess
        if 'b' not in self.ic:
            self.ic['b'] = slope_guess
        if 's' not in self.ic:
            self.ic['s'] = abs(self.x[-1] - self.x[0]) / 10
        if self.function is gaussian_quad_turn:
            if 'c' not in self.ic:
                self.ic['c'] = 0.0
            if 'd' not in self.ic:
                self.ic['d'] = 0.0

    def guess_limit(self):
        slope_guess = (self.y[-1] - self.y[0]) / (self.x[-1] - self.x[0])
        if 'X0' not in self.limits:
            self.limits['X0'] = (min(self.x), max(self.x))
        if 'Y0' not in self.limits:
            self.limits['Y0'] = (min(self.y), max(self.y))
        if 'a' not in self.limits:
            self.limits['a'] = (-10 * abs(slope_guess), 10 * abs(slope_guess))
        if 'b' not in self.limits:
            self.limits['b'] = (-10 * abs(slope_guess), 10 * abs(slope_guess))
        if 's' not in self.limits:
            if 's' in self.ic:
                ic = float(self.ic['s'])
            else:
                ic = abs(self.x[-1] - self.x[0]) / 10
            self.limits['s'] = (0.1 * ic, 100 * ic)

        if self.function is gaussian_quad_turn:
            quad_limit = slope_guess / np.mean(self.limits['X0'])
            if 'c' not in self.limits:
                self.limits['c'] = (-quad_limit,quad_limit)
            if 'd' not in self.limits:
                self.limits['d'] = (-quad_limit,quad_limit)

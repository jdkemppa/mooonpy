# -*- coding: utf-8 -*-
from mooonpy.fitting.fitting import CurveFit
from scipy.stats import norm
import numpy as np
from matplotlib import pyplot as plt
from lmfit import Parameters

def symmetric_von_mises(x, center, sigma, amplitude=1, rad=False):
    if not rad:
        sigma = np.deg2rad(sigma)
    k = 1 / (sigma * sigma)

    if rad is True:
        B = np.i0(k) * 4 * np.pi
        y = (np.exp(k * np.cos(x - center))+np.exp(k * np.cos(x + center))) / B
    else:
        B = np.i0(k) * 4 * 180
        y = (np.exp(k * np.cos(np.deg2rad(x - center))) +np.exp(k * np.cos(np.deg2rad(x + center))))/ B

    # y += y[::-1]

    return y * amplitude


def symmetric_gaussian_periodic(x, center, sigma, amplitude, rad=False):
    """
    Symmetric Gaussian with proper periodicity handling.
    Creates Gaussians at ±center and their periodic images for bleed-through.
    amplitude = ?
    https://en.wikipedia.org/wiki/Wrapped_normal_distribution
    """
    y = np.zeros_like(x)

    y += amplitude / 2 * norm.pdf(x, center, sigma)
    y += amplitude / 2 * norm.pdf(x, -center, sigma)

    # Add periodic images for bleed-through
    if rad:
        y += amplitude / 2 * norm.pdf(x, center - np.pi * 2, sigma)
        y += amplitude / 2 * norm.pdf(x, -center + np.pi * 2, sigma)
    else:
        y += amplitude / 2 * norm.pdf(x, center - 360, sigma)
        y += amplitude / 2 * norm.pdf(x, -center + 360, sigma)

    return y


def multi_symmetric_gaussian(x, **params):
    """
    Sum of multiple symmetric Gaussians with periodicity
    """
    y = np.zeros_like(x)

    i = 1
    while f'center_{i}' in params:  # change this to not assume 1?
        mu = params[f'center_{i}']
        sigma = params[f'sigma_{i}']
        amplitude = params[f'amp_{i}']

        y += symmetric_gaussian_periodic(x, mu, sigma, amplitude)
        i += 1

    return y

def multi_symmetric_von_mises(x, **params):
    """
    Sum of multiple symmetric Von Mises
    """
    y = np.zeros_like(x)

    i = 1
    while f'center_{i}' in params:  # change this to not assume 1?
        mu = params[f'center_{i}']
        sigma = params[f'sigma_{i}']
        amplitude = params[f'amp_{i}']

        y += symmetric_von_mises(x, mu, sigma, amplitude)
        i += 1

    return y

class DihedralDist(CurveFit):
    """
    Operators on dihedral distributions for modeling energy minima
    """
    default_fit_kws = {
        'method': 'powell',
        'max_nfev': 5000
    }

    def __init__(self, phi_angles, bin_scale=10, name=None, function=multi_symmetric_gaussian, ic=None, limits=None):
        self.phi_angles = phi_angles
        self.bin_scale = bin_scale

        self.bins = round(len(self.phi_angles) / self.bin_scale)  # count
        self.hist, self.bin_edges = np.histogram(self.phi_angles, bins=self.bins, density=False)
        self.bin_width = self.bin_edges[1] - self.bin_edges[0]

        self.bin_centers = (self.bin_edges[:-1] + self.bin_edges[1:]) / 2
        self.bin_edges = np.repeat(self.bin_edges, 2)  # overwrite with both edges
        self.hist_edges = np.concatenate([[0], np.repeat(self.hist, 2), [0]])  # copy to make step histogram
        self.amp2count = self.bin_width * len(self.phi_angles)
        super().__init__(self.bin_centers, self.hist, name, function, ic, limits)  # makes x and y attributes

    def guess_ic(self, guess=None):
        """
        Basic guessing function.
        could make smart with peak detection but this is good enough for now
        """
        if self.function is multi_symmetric_gaussian or self.function is multi_symmetric_von_mises:
            if guess == 1:  # could do is int and some spacing
                self.ic['center_1'] = 90
                self.ic['sigma_1'] = 30
                self.ic['amp_1'] = self.amp2count
            elif guess == 2:
                self.ic['center_1'] = 0
                self.ic['sigma_1'] = 30
                self.ic['amp_1'] = self.amp2count / 2
                self.ic['center_2'] = 180
                self.ic['sigma_2'] = 30
                self.ic['amp_2'] = self.amp2count / 2
            elif guess == 3 or guess is None:
                self.ic['center_1'] = 0
                self.ic['sigma_1'] = 30
                self.ic['amp_1'] = self.amp2count / 3
                self.ic['center_2'] = 180
                self.ic['sigma_2'] = 30
                self.ic['amp_2'] = self.amp2count / 3
                self.ic['center_3'] = 90
                self.ic['sigma_3'] = 30
                self.ic['amp_3'] = self.amp2count / 3
            elif isinstance(guess, dict):
                for key, value in guess.items():
                    self.ic[key] = value
                    if key.startswith('center'):
                        N = int(key[-1])  # number
                        if f'amp_{N}' not in guess:
                            self.ic[f'amp_{N}'] = self.amp2count / len(guess)
                        if f'sigma_{N}' not in guess:
                            self.ic[f'sigma_{N}'] = 30

            else:
                pass

    def guess_limit(self):
        for key in self.ic.keys():
            if key.startswith('center'):
                N = int(key[-1])  # number
                if key not in self.limits:
                    self.limits[key] = (0, 180)
                if f'amp_{N}' not in self.limits:
                    self.limits[f'amp_{N}'] = (0, self.amp2count * 1.1)
                if f'sigma_{N}' not in self.limits:
                    self.limits[f'sigma_{N}'] = (5, 180)

    def make_params(self):
        """
        Override base make_params function with extra amplitude constraint
        """
        params = Parameters()
        amps = [param for param in self.ic.keys() if param.startswith('amp_')]
        for param, value in self.ic.items():
            limit = self.limits.get(param, self.def_limit)
            if len(amps) > 1 and param == amps[-1]:  # constraint mode for total area on last param
                expr = f'{self.amp2count}'
                for amp in amps[:-1]:
                    expr += f' - {amp}'
                # print(expr) # amp_3 = 10000 - amp_1 - amp_2
                params.add(param, min=limit[0], max=limit[1], expr=expr)
                ## last amp is dependent such that area is preserved from ic
            else:
                params.add(param, value=float(value), min=limit[0], max=limit[1],
                           vary=not isinstance(value, str))  # no vary if value is string
        self.params = params

    def plot_mode_decomposition(self, ax=None, figsize=(10, 6)):
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)
            ax.set(xlabel='Phi Angle [Degrees]', ylabel='Count', xlim=(-180, 180), ylim=(0, max(self.hist) * 1.01),
                   title=self.name)

        ax.plot(self.bin_edges, self.hist_edges, color='gray',
                label=f'Data: {len(self.phi_angles)} Dihedrals, {self.bins} Bins')

        bin_centers = self.bin_centers

        y_fit = self.result.eval(x=bin_centers)
        N = len(self.fitted_params.keys()) // 3  # 3 params per mode
        if N > 1:
            ax.plot(bin_centers, y_fit, 'r-', linewidth=2, label=f'Total Fit: Area={round(np.sum(y_fit))}')

        colors = ['blue', 'green', 'orange', 'purple', 'brown']
        for ii in range(1, N + 1):
            center = self.fitted_params[f'center_{ii}']
            sigma = self.fitted_params[f'sigma_{ii}']
            amplitude = self.fitted_params[f'amp_{ii}']
            if self.function is multi_symmetric_gaussian:
                y_component = symmetric_gaussian_periodic(bin_centers, center, sigma, amplitude)
            elif self.function is multi_symmetric_von_mises:
                y_component = symmetric_von_mises(bin_centers, center, sigma, amplitude)
            else:
                y_component = self.function(bin_centers, center, sigma, amplitude)

            label = f"Mode {ii}: ±{center:.1f}° (σ={sigma:.1f}° Area={round(np.sum(y_component))}) )"

            ax.plot(bin_centers, y_component, '--', color=colors[ii - 1 % len(colors)],
                    linewidth=1.5, label=label)

        ax.legend()

        return ax


if __name__ == '__main__':
    import mooonpy as mp
    mol = mp.Molspace('..\\..\\..\\examples\\EPON_862\\lmp_dump\\small_epon_densify.data') # included in git
    bond_hist, angle_hist, dihedral_hist, improper_hist = mol.compute_BADI_by_type()

    # test_dist = dihedral_hist['cp-cp-cp-hc']
    test_dist = []
    for list_ in dihedral_hist.values():
        test_dist += list_
    dihed_dist = DihedralDist(test_dist)
    dihed_dist.run_default()
    dihed_dist.plot_mode_decomposition()

    vm_dist = DihedralDist(test_dist,function=multi_symmetric_von_mises)
    # vm_dist.amp2count = vm_dist.amp2count/np.pi # close to the right value, 3 and pi are close but not exact?
    vm_dist.guess_ic()
    # vm_dist.ic.update({'center_1':'0','center_2':'180'}) # hardcode 2, last is free
    vm_dist.run_default()
    vm_dist.plot_mode_decomposition()
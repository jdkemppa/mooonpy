# -*- coding: utf-8 -*-
"""
Dihedral angle distribution analysis and fitting tools
Cite PBZ paper by Tristan Muzzy
"""
from mooonpy.fitting.fitting import CurveFit
import numpy as np
from mooonpy.tools.misc_utils import calculate_axis_ticks
import warnings

def symmetric_von_mises(x, center, sigma, amplitude=1.0, rad=False):
    """
    Periodic von Mises fit, from -pi to pi or -180 to 180
    Symmetry of dihedral angle PDF is enforced

    Note under ~2.3 deg sigma uses gaussain to prevent float overflow
    """
    if not rad:
        sigma = np.deg2rad(sigma)
    k = 1 / (sigma * sigma)
    if k > 625: # float overflow limit of bessel function
        warnings.warn('symmetric von Mises fit cannot be computed with high kappa, using gaussian eval')
        if not rad:
            sigma = np.rad2deg(sigma)
        return symmetric_gaussian_periodic(x, center, sigma, amplitude, rad)

    try:
        if rad is True:
            B = np.i0(k) * 4 * np.pi
            y = (np.exp(k * np.cos(x - center)) + np.exp(k * np.cos(x + center))) / B
        else:
            B = np.i0(k) * 4 * 180
            y = (np.exp(k * np.cos(np.deg2rad(x - center))) + np.exp(k * np.cos(np.deg2rad(x + center)))) / B
    except:
        print(k)
        raise Exception('foo')
    return y * amplitude


def symmetric_gaussian_periodic(x, center, sigma, amplitude, rad=False):
    """
    Symmetric Gaussian with proper periodicity handling.
    Creates Gaussian's at ±center and their first periodic images for bleed-through.
    https://en.wikipedia.org/wiki/Wrapped_normal_distribution
    """
    from scipy.stats import norm
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
    Sum of multiple symmetric Gaussian's with periodicity
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
    Sum of multiple symmetric von Mises
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
    Fit and analyze the dihedral angle distributions with histograms and probability density functions.
    using fit_model curve fits histogram bins
    using fit_mle uses MLE (recommended/default)
    """
    default_fit_kws = {
        'method': 'powell',
        'max_nfev': 5000,

    }
    default_limits = {
        'sigma': (5, 90),
        'mu': (0, 180)
    }

    def __init__(self, phi_angles, bin_scale=10, name=None, function=multi_symmetric_gaussian, ic=None, limits=None):
        self.phi_angles = np.array(phi_angles)
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
        could make smart with peak detection or BIC but this is good enough for now
        Guess may be int 1-3 for number of peaks to use a default, or a dict of custom ICs
        """
        if self.function is multi_symmetric_gaussian or self.function is multi_symmetric_von_mises:
            if guess == 1:  # could do isint(guess) and some spacing thing to compact
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
            elif isinstance(guess, dict):  # unpack custom guesses
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
        """
        Add parameter limits. May be overridden outside
        """
        for key in self.ic.keys():
            if key.startswith('center'):
                N = int(key[-1])  # number
                if key not in self.limits:
                    self.limits[key] = self.default_limits['mu']
                if f'amp_{N}' not in self.limits:
                    self.limits[f'amp_{N}'] = (0, self.amp2count * 1.1)
                if f'sigma_{N}' not in self.limits:
                    self.limits[f'sigma_{N}'] = self.default_limits['sigma']

    def make_params(self):
        """
        Override base class make_params function with extra amplitude constraint
        Initial area under the curve and sum of amplitudes=1 is preserved
        Note: this is only relevant for fitting histograms, MLE constrains using scaling.
        """
        from lmfit import Parameters
        params = Parameters()
        amps = [param for param in self.ic.keys() if param.startswith('amp_')]
        for param, value in self.ic.items():
            limit = self.limits.get(param, self.def_limit)
            if param == amps[-1]:  # constraint mode for total area on last param
                expr = f'{self.amp2count}'
                for amp in amps[:-1]:
                    expr += f' - {amp}'
                # print(expr) # amp_3 = 10000 - amp_1 - amp_2
                params.add(param, min=limit[0], max=limit[1], expr=expr)
                ## ^ last amp is dependent such that area is preserved from ic
            else:
                params.add(param, value=float(value), min=limit[0], max=limit[1],
                           vary=not isinstance(value, str))  # no vary if value is string
        self.params = params

    def plot_mode_decomposition(self, ax=None, figsize=(10, 6)):
        """
        Make pretty plots
        """
        import matplotlib.pyplot as plt
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)

        else:
            fig = None

        ax.plot(self.bin_edges, self.hist_edges, color='gray',
                label=f'Data: {len(self.phi_angles)} Dihedrals, {self.bins} Bins')

        bin_centers = self.bin_centers

        y_fit = self.function(bin_centers, **self.fitted_params)

        N = len(self.fitted_params.keys()) // 3  # 3 params per mode
        if N > 1:
            ax.plot(bin_centers, y_fit, 'k-', linewidth=2, label=f'Total Fit: Area={round(np.sum(y_fit))}')

        colors = ['#FFB300', 'brown', 'blue', 'green', 'orange', 'purple']
        color_cyle = 0
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

            if N > 1:
                label = f"Mode {ii}: ±{center:.1f}°; σ={sigma:.1f}°; Area={round(np.sum(y_component))}→{round(np.sum(y_component) / np.sum(y_fit) * 100)}%"
            else:
                label = f"Mode {ii}: ±{center:.1f}°; σ={sigma:.1f}°; Area={round(np.sum(y_component))}"
            if center == 0:
                color = '#00CED1'  # turquoise
            elif center == 180:
                color = '#D946A6'  # magenta
            else:
                color = colors[color_cyle]  # amber
                color_cyle += 1

            # ax.plot(bin_centers, y_component, '--', color=colors[ii - 1 % len(colors)],
            #         linewidth=1.5, label=label)
            ax.plot(bin_centers, y_component, '--', color=color,
                    linewidth=1.5, label=label)

            ax.set(xlabel='Phi Angle [Degrees]', ylabel='Count', xlim=(-180, 180), xticks=np.arange(-180, 181, 60),
                   ylim=(0, max(self.hist) * 1.01),
                   title=self.name)
            try:
                yticks, ylim, spacing = calculate_axis_ticks(0, max(self.hist))
                ax.set(ylim=[0, ylim], yticks=yticks)
            except:
                pass

        ax.legend()

        return fig, ax

    def fit_mle(self):
        """
        Fit using Maximum Likelihood Estimation on raw angle data
        Returns fitted parameters in the same format as fit_model()
        Note: lmfit.minimize doesn't support expr constraints, so we normalize in post-processing
        """
        from lmfit import minimize
        # Use lmfit.minimize to minimize negative log-likelihood
        result = minimize(self._negative_log_likelihood, self.params, **self.default_fit_kws)

        # Store result (lmfit.MinimizerResult instead of ModelResult)
        self.result = result

        # Extract fitted parameters and normalize amplitudes to sum to amp2count
        fitted_params_raw = {name: result.params[name].value for name in result.params}

        # Find all amplitude parameters and normalize them
        amp_keys = sorted([k for k in fitted_params_raw.keys() if k.startswith('amp_')])
        raw_amp_sum = sum(fitted_params_raw[k] for k in amp_keys)

        # Rescale amplitudes to sum to amp2count
        fitted_params = {}
        for name, value in fitted_params_raw.items():
            if name.startswith('amp_'):
                fitted_params[name] = value * (self.amp2count / raw_amp_sum) if raw_amp_sum > 0 else value
            else:
                fitted_params[name] = value

        self.fitted_params = fitted_params

        return self.fitted_params, result

    def run_default(self):
        """
        Override to use MLE by default for DihedralDist
        """
        if not self.ic:
            self.guess_ic()
        if not self.ic:  # if above cannot make IC ...
            raise Exception('ic is required, no parameters to fit')
        if not self.limits:
            self.guess_limit()
        self.make_params()
        fitted_params, result = self.fit_mle()
        return fitted_params, result

    def _convert_amps_to_weights(self, params_dict):
        """
        Convert amplitude parameters to weights (sum to 1) for PDF calculation
        Enforces constraint by normalizing all amplitudes by their sum
        """
        # Find all amp parameters
        amp_keys = sorted([k for k in params_dict.keys() if k.startswith('amp_')])
        amp_sum = sum(params_dict[k] for k in amp_keys)

        # Normalize to get weights that sum to 1
        weight_params = {}
        for key, value in params_dict.items():
            if key.startswith('amp_'):
                weight_params[key] = value / amp_sum
            else:
                weight_params[key] = value

        return weight_params

    def _negative_log_likelihood(self, params):
        """
        Negative log-likelihood for MLE optimization
        """
        # Convert Parameters object to dict and then to weights
        params_dict = {name: params[name].value for name in params}
        weight_params = self._convert_amps_to_weights(params_dict)

        # Calculate PDF for all angles using existing functions with weights as amplitudes
        pdf_values = self.function(self.phi_angles, **weight_params)

        # Avoid log(0) by adding small epsilon
        pdf_values = np.maximum(pdf_values, 1e-10)

        # Return negative log-likelihood
        nll = -np.sum(np.log(pdf_values))

        return nll


if __name__ == '__main__':
    import mooonpy as mp

    mol = mp.Molspace('..\\..\\..\\examples\\EPON_862\\lmp_dump\\small_epon_densify.data')  # included in git
    bond_hist, angle_hist, dihedral_hist, improper_hist = mol.compute_BADI_by_type()

    # test_dist = dihedral_hist['cp-cp-cp-hc']
    test_dist = []
    for list_ in dihedral_hist.values():
        test_dist += list_
    dihed_dist = DihedralDist(test_dist)
    dihed_dist.run_default()
    fig1, ax1 = dihed_dist.plot_mode_decomposition()
    ax1.grid()

    vm_dist = DihedralDist(test_dist, function=multi_symmetric_von_mises)
    # vm_dist.amp2count = vm_dist.amp2count/np.pi # close to the right value, 3 and pi are close but not exact?
    vm_dist.guess_ic()
    # vm_dist.ic.update({'center_1':'0','center_2':'180'}) # hardcode 2, last is free
    vm_dist.run_default()
    fig2, ax2 = vm_dist.plot_mode_decomposition()
    ax2.grid()

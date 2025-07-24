# -*- coding: utf-8 -*-


import numpy as np
from scipy.signal import find_peaks

from numbers import Number
from typing import Union, Tuple, Optional

Array1D = Union[np.ndarray, list]


def aggregate_fun(fun_name: str, vector: Array1D) -> Number:
    """
    Apply generic aggregate functions to a vector, similar to the bottom row of a spreadsheet.
    Modes control values to be ignored in operation, by slicing the vector before the operation.

    Note that operating on an empty sclice will return a np.nan.

    **Currently supported Operators: fun_name must start with one of these substrings**
        - sum: Summation
        - avg: Arithmetic mean
        - std: Standard Deviation

    **Currently supported Modes: if fun_name ends with one of these substrings, the rule is applied**
        -0: Ignores 0's
        -None: Ignores None or nan

    :param fun_name: Selection of Function and Mode.
    :type fun_name: str
    :param vector: row or column vector to operate on.
    :type vector: list or numpy array
    :return: Result of operation.
    :rtype: Number

    :Example:
        >>> from mooonpy.tools import aggregate_fun
        >>> MyVect = [0,1,2]
        >>> avg_no0 = aggregate_fun('avg0', MyVect)
        1.5
        >>> avg_all = aggregate_fun('avg', MyVect)
        1.0
    """
    # Do sclicing for Mode rules
    if fun_name.endswith('0'):
        slice_ = np.nonzero(vector)
        if len(slice_) == 0:
            return np.nan
    elif fun_name.endswith('None'):
        slice_ = np.logical_not(np.isnan(vector))
        if len(slice_) == 0:
            return np.nan
    else:
        slice_ = True

    # Compute aggregate
    if fun_name.startswith('sum'):
        return np.sum(vector[slice_])
    elif fun_name.startswith('avg'):
        return np.mean(vector[slice_])
    elif fun_name.startswith('std'):
        return np.std(vector[slice_])
    else:
        raise Exception(f'ERROR: fun_name {fun_name} not recognized')


def find_peaks_and_valleys(xdata: Array1D, ydata: Array1D, prominence: Optional[Number] = None) -> Tuple[
    np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute the peaks and valleys of a curve using a specified prominence for cutoff.

    :param xdata: Array of x values. Must be increasing monotonically.
    :type xdata: Array1D
    :param ydata: Array of y values.
    :type ydata: Array1D
    :param prominence: Prominence for cutoff
    :type prominence: Number or None

    :return: X peaks, Y peaks, X valleys, Y valleys
    :rtype: Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]

    .. TODO::
        make an example image
    """
    # Find peaks
    peaks, properties = find_peaks(ydata, prominence=prominence)
    xpeaks = xdata[peaks]
    ypeaks = ydata[peaks]

    # Find valleys
    xvalleys, yvalleys = [], []
    if len(xpeaks) >= 2 and len(ypeaks) >= 2:
        for i in range(len(peaks) - 1):
            lo = peaks[i]
            hi = peaks[i + 1]
            between_peaksx = xdata[lo:hi]
            between_peaksy = ydata[lo:hi]
            minimum_index = np.min(np.where(between_peaksy == between_peaksy.min())[0])
            xvalleys.append(between_peaksx[minimum_index])
            yvalleys.append(between_peaksy[minimum_index])
    xvalleys = np.array(xvalleys)
    yvalleys = np.array(yvalleys)
    return xpeaks, ypeaks, xvalleys, yvalleys


def compute_derivative(xdata: Array1D, ydata: Array1D) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Function to compute the 1st and 2nd order central derivatives.
    Edges are not considered, so trimmed x array is returned.

    :param xdata: Array of x values.
    :type xdata: Array1D
    :param ydata: Array of y values.
    :type ydata: Array1D
    :return: Trimmed x, 1st derivative, 2nd derivative
    :rtype: Tuple[np.ndarray, np.ndarray, np.ndarray

    .. TODO::
        make an example image
    """
    dxn = xdata[1:-1]  # ignore first and last point
    dxn = dxn.astype(float)
    dy1 = np.zeros_like(dxn)
    dy2 = np.zeros_like(dxn)
    if len(xdata) == len(ydata):
        for i in range(1, len(xdata) - 1):
            dx = (xdata[i + 1] - xdata[i - 1]) / 2
            if dx == 0:
                print(
                    'WARNING finite difference dx was zero at x={}. Derivative was set to zero to avoid infinite derivative.'.format(
                        xdata[i]))
            else:
                dy1[i - 1] = (ydata[i + 1] - ydata[i - 1]) / (2 * dx)
                dy2[i - 1] = (ydata[i + 1] - 2 * ydata[i] + ydata[i - 1]) / (dx * dx)
    else:
        print('ERROR (compute_derivative) inconsistent number of data points between X and Y arrays')
    return dxn, dy1, dy2


def compute_fringe_slope(strain: np.ndarray, stress: np.ndarray, min_strain: Optional[Number] = None,
                         max_strain: Optional[Number] = None, direction: str = 'forward') -> Tuple[
    np.ndarray, np.ndarray]:
    """
    Compute fringe slope

    .. TODO::
        Most concise writeup and a link to the paper

    :param strain: Strain
    :type strain: Array1D
    :param stress: Stress
    :type stress: Array1D
    :param min_strain: Minimum strain, defaults to None to use minimum index
    :type min_strain: Optional(Number)
    :param max_strain: Maximum strain, defaults to None to use maximum index
    :type max_strain: Optional(Number)
    :param direction: Direction, Must be 'forward' or 'backward', flips both vectors, defaults to 'forward'.
    :type direction: str

    :return: Fringe slope X and Y
    :rtype: Tuple[np.ndarray, np.ndarray]
    """
    # Set direction
    if direction == 'forward':
        strain = strain.copy()
        stress = stress.copy()
    elif direction == 'reverse':
        strain = np.flip(strain.copy())
        stress = np.flip(stress.copy())
    else:
        raise Exception(
            f'ERROR direction={direction} is not supported. Supported directions are "forward" or "reverse"')

    # Set defaults if min_strain or max_strain are None
    if min_strain is None: min_strain = min(strain)
    if max_strain is None: max_strain = max(strain)

    ## Vectorized is much quicker but there's an off by 1 error
    #
    # n = np.arange(0, len(strain))
    # slice_ = (min_strain <= strain) & (strain <= max_strain) & (n >= 4)
    #
    # sum_xi = np.cumsum(strain)
    # sum_yi = np.cumsum(stress)
    # sum_xi_2 = np.cumsum(strain*strain)
    # sum_yi_2 = np.cumsum(stress*stress)
    # sum_xi_yi = np.cumsum(strain*stress)
    #
    # sum_xi = sum_xi[slice_]
    # sum_yi = sum_yi[slice_]
    # sum_xi_2 = sum_xi_2[slice_]
    # sum_yi_2 = sum_yi_2[slice_]
    # sum_xi_yi = sum_xi_yi[slice_]
    # n = n[slice_]
    # fringe = strain[slice_]
    #
    # SSxy = sum_xi_yi - (sum_xi * sum_yi / n)
    # SSxx = sum_xi_2 - (sum_xi * sum_xi / n)
    # slopes = SSxy / SSxx
    #
    # return fringe, slopes

    ## Original non-vectorized
    # Start the walked linear regression method
    slopes, fringe = [], []
    sum_xi, sum_yi, sum_xi_2, sum_yi_2, sum_xi_yi, n = 0, 0, 0, 0, 0, 0
    for x, y in zip(strain, stress):
        # Compute cumulative linear regression parameters
        sum_xi += x
        sum_yi += y
        n += 1
        sum_xi_2 += x * x
        sum_yi_2 += y * y
        sum_xi_yi += x * y

        # Need at least 2 points to perform linear regression
        if n <= 3: continue

        # Only compute outputs if x is in the desired range
        if min_strain <= x <= max_strain:
            SSxy = sum_xi_yi - (sum_xi * sum_yi / n)
            SSxx = sum_xi_2 - (sum_xi * sum_xi / n)
            b1 = SSxy / SSxx

            slopes.append(b1)
            fringe.append(x)

    return np.array(fringe), np.array(slopes)


def first_value_cross(xdata: Array1D, ydata: Array1D, cross: Optional[Number] = None):
    """
    Find x location of the first time the y data crosses a value of the y data.

    :param xdata: X data
    :type xdata: Array1D
    :param ydata: Y data
    :type ydata: Array1D
    :param cross: Value to compare against, defaults to None, which uses the mean value of ydata
    :type cross: Number
    :return: First x location
    :rtype: Number
    """
    if cross is None:
        cross = np.mean(ydata)
    x_cross = xdata[np.min(np.where(ydata < cross)[0])]
    return x_cross

def hyperbola(x, X0, Y0, a, b, c):
    """
    Patrone Hyperbola fit method

    :param x: X data
    :type x: Array1D
    :param X0: center x coordinate
    :type X0: Number
    :param Y0: center y coordinate
    :type Y0: Number
    :param a: slope left of turn
    :type a: Number
    :param b: slope increase after turn
    :type b: Number
    :param c: slope curvature parameter
    :type c: Number
    """
    dx = x - X0
    h0 = 0.5 * dx + np.sqrt((dx * dx) / 4 + np.exp(c))
    p = Y0 + a * dx + b * h0
    return p

class MixingRule():
    def __init__(self, style='add'):
        """
        Initialize the mixing rule with a specific method name.

        :param style: The name of the addition method to use
        """
        self.style = style.lower()

        self._method_map = {
            'add': self._arithmetic,
            'diameter': self._diameter,
            'arithmetic': self._arithmetic,
            'geometric': self._geometric,
            'sixth_epsilon': self._sixth_epsilon,
            'sixth_sigma': self._sixth_sigma,
        }
    def __call__(self,*args, **kwargs):
        return self._method_map[self.style](*args, **kwargs)

    def _add(self,*args):
        return sum(args)

    def _diameter(self,*args):
        return sum(args)/2

    def _arithmetic(self,*args):
        return sum(args)/len(args)

    def _geometric(self,*args):
        return np.power(np.prod(args),1/len(args))

    def _sixth_epsilon(self,*args):
        """
        LAMMPS mix sixthpower function for epsilon.
        Assumes ep1, sig1 ep2, sig2 ordering
        """
        return 2*np.sqrt(args[0],args[2])*args[1]**3*args[3]**3/(args[1]**6+args[3]**6)
    def _sixth_sigma(self,*args):
        if len(args) == 2:
            sig1, sig2 = args
        elif len(args) == 4: # order from above function
            sig1 = args[1]
            sig2 = args[3]
        else:
            raise Exception()
        return np.power(0.5*(sig1**6+sig2**6),1/6)
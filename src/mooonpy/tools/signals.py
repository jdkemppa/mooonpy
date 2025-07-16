from numbers import Number
from typing import Tuple, Union, Optional

import numpy as np
import scipy as sp

from mooonpy.tools.math_utils import first_value_cross

Array1D = Union[np.ndarray, list]


def compute_PSD(xdata: Array1D, ydata: Array1D) -> Tuple[Array1D, Array1D]:
    """
    Compute the Power Spectral Density (PSD) with the X-values being the normalized frequencies.

    Returns the normalized frequencies (Omega N / Wn) and power spectral density (PSD) as a tuple of arrays.

    :param xdata: Array of x values.
    :type xdata: Array1D
    :param ydata: Array of y values.
    :type ydata: Array1D
    :return: wn's, PDS's
    :rtype: Tuple(Array1D, Array1D)
    """
    # Define sampling rate and number of data points
    dx = np.mean(np.abs(np.diff(xdata)))
    if dx != 0:
        fs = 1 / dx  # sampling rate
    else:
        fs = xdata.shape[0] / (np.max(xdata) - np.min(xdata))
    N = xdata.shape[0]  # number of data points
    d = 1 / fs  # sampling space

    # Perform one sided FFT
    fft_response = np.fft.rfft(ydata, axis=0, norm='backward')
    x_fft = np.fft.rfftfreq(N, d=d)
    y_fft = fft_response

    # Compute the final PSD and normalized cutoff frequencies
    psd = np.real((y_fft * np.conjugate(y_fft)) / N)
    wns = x_fft / (0.5 * fs)
    return wns, psd


def power_to_db(power: Array1D, ref_power: Number = 1) -> np.ndarray:
    """
    Convert power to decibels (dB).

    :param power: Power to convert.
    :type power: Array1D
    :param ref_power: Reference power.
    :type ref_power: Number
    :return: Converted power.
    :rtype: Array1D
    """
    power_in_dB = 10 * np.log10(power / ref_power)
    return power_in_dB


def _butter(ydata: Array1D, wn: float, order: Number) -> np.ndarray:
    """
    Alias for scipy's butterworth low-pass filter.

    :param ydata: Y data
    :type ydata: Array1D
    :param wn: Normalized Cutoff frequency, defaults to 'PSD' optimization
    :type wn: Number
    :param order: Order of filter, defaults to 2
    :type order: Number
    """
    sos = sp.signal.butter(order, wn, btype='low', analog=False, output='sos', fs=None)
    y_filt = sp.signal.sosfiltfilt(sos, ydata, padtype=None)
    return y_filt


def butter_lowpass(xdata: Optional[Array1D], ydata: Array1D, wn: Union[Number, str] = 'PSD', quadrant: str = 'msr',
                   order: Number = 2) -> Tuple[np.ndarray,Number,str]:
    """
    Call butterworth low-pass filter with Quadrant Mirroring padding.
    The original image region is returned fully filtered.

    Cutoff frequency may be input as a number, or can be chosen using an optimizer style

    If the quadrant is 'msr' the lowest residual low and high data_extension is used,
    or the quadrants are specified with a comma separated string of '2,3' for left and right extension.

    :param xdata: Array of x values. Only used for PSD optimization
    :type xdata: Array1D
    :param ydata: Y data to filter.
    :type ydata: Array1D
    :param wn: Normalized Cutoff frequency, defaults to 0.1.
    :type wn: Number
    :param quadrant: Quadrants, defaults to 'msr' optimization.
    :type quadrant: str
    :param order: Order of filter, defaults to 2
    :type order: Number

    :return: Filtered data, wn used, qm string used
    """
    # Scrub inputs
    # ------------------------------------
    if xdata is None:
        xdata = np.arange(len(ydata))
    else:
        xdata = xdata.copy()
    ydata = ydata.copy()

    if wn == 'PSD':
        wns, psd = compute_PSD(xdata, ydata)
        wn = first_value_cross(wns, psd)
        if wn == 0: wn=wns[1] # use the lowest allowed value

    if quadrant == 'msr':
        quad_lo, quad_hi = determine_mirroring_locations(xdata, ydata, wn, order)
    elif len(quadrant) == 3 and ',' in quadrant:
        splits = quadrant.split(',')
        quad_lo = int(splits[0])
        quad_hi = int(splits[1])
    elif len(quadrant) == 2:
        quad_lo, quad_hi = int(quadrant[0]), int(quadrant[1])  # works for tuple or string
    else:
        raise Exception(f'quadrant {quadrant} must be either "msr" or comma-separated quadrants i.e. "1,3"')

    # do filtering
    # ------------------------------------
    xdata, ydata, lo_trim, hi_trim = data_extension(xdata, ydata, quad_lo, quad_hi)

    y_butter = _butter(ydata, wn, order)
    y_filt = y_butter[lo_trim:hi_trim]
    return y_filt, wn, f'{quad_lo},{quad_hi}'


def determine_mirroring_locations(xdata: Optional[Array1D], ydata: Array1D, wn: Union[Number, str] = 'PSD',
                                  order: Number = 2) -> Tuple[int, int]:
    """
    Compare 4 mirror images at both ends of dataset.
    The combination with the lowest residual^2 after filtering is returned.

    :param xdata: X data
    :type xdata: Array1D
    :param ydata: Y data
    :type ydata: Array1D
    :param wn: Normalized Cutoff frequency, defaults to 'PSD' optimization.
    :type wn: Number or str
    :param order: Order of filter, defaults to 2
    :type order: Number

    :return: Lowest residual and Highest residual.
    :rtype: Tuple(int, int)
    """
    if xdata is None:
        xdata = np.arange(len(ydata))
    # Determine half_data to only check for residuals
    # either from lo-half_data or half_data-hi
    half_data = int(xdata.shape[0] / 2)

    # Determine optimal wn with PSD method. Used as constant between all mirrors
    if wn == 'PSD':
        wns, psd = compute_PSD(xdata, ydata)
        wn = first_value_cross(wns, psd, cross=None)  # use mean with default

    # ------------------------------------
    # First: Optimize the "lo" end
    lo_quads2test = [1, 2, 3, 4];
    lo_summed_residuals2 = {}  # {quadrant_mirror:sum-of-residuals-squared}
    for quad in lo_quads2test:
        x_quad, y_quad, lo_ind, hi_ind = data_extension(xdata, ydata, quad, 1)
        y_butter = _butter(y_quad, wn, order)
        y_filt = y_butter[lo_ind:hi_ind]  # get original slice
        residuals = ydata - y_filt
        residuals = residuals[:half_data]  # we only care about the first half fit
        lo_summed_residuals2[quad] = np.sum(residuals ** 2)

    # Find minimized sum of residuals squared
    lo = min(lo_summed_residuals2, key=lo_summed_residuals2.get)

    # ------------------------------------
    # Second: Optimize the "hi" end
    hi_quads2test = [1, 2, 3, 4];
    hi_summed_residuals2 = {}  # {quadrant_mirror:sum-of-residuals-squared}
    for quad in hi_quads2test:
        x_quad, y_quad, lo_ind, hi_ind = data_extension(xdata, ydata, 1, quad)
        y_butter = _butter(y_quad, wn, order)
        y_filt = y_butter[lo_ind:hi_ind]  # get original slice
        residuals = ydata - y_filt
        residuals = residuals[half_data:]  # we only care about the last half fit
        hi_summed_residuals2[quad] = np.mean(residuals ** 2)  # account for 1 being shorter

    # Find minimized sum of residuals squared
    hi = min(hi_summed_residuals2, key=hi_summed_residuals2.get)

    # ------------------------------------#
    # return optimal quadrant_mirror
    return lo, hi


def data_extension(xdata: Array1D, ydata: Array1D, lo: Number = 1, hi: Number = 1) -> Tuple[
    np.ndarray, np.ndarray, int, int]:
    """
    Function to extend data at the "lo" and "hi" end. This function assumes xdata is increasing and
    ydata is the dependent variable. This function replaces the usage of scipy's "signal.filtfilt"
    or "signal.sosfiltfilt" padtype=<'odd' or 'even' or 'const' or None>, with greater functionality.

    For example the padtype=<'odd' or 'even' or 'const' or None> implementation assumes that the
    padding on the "lo X-end" will be the same as the padding that should be used on the "hi X-end",
    however for stress-strain data, this is usually not the case.

    lo/hi integer meanings can be thought of as quadrant numbers for "lo" and have the same meaning
    for "hi", but at the tail end of the data. For example:

        - 1 = means shut off data mirroring
        - 2 = means mirror data (lo=mirrored into quadrant 2)
        - 3 = means apply 180-degree rotation (lo=mirrored into quadrant 3)
        - 4 = means flip the Y-data and append at the end of the Y-data

    :param xdata: X data
    :type xdata: Array1D
    :param ydata: Y data
    :type ydata: Array1D
    :param lo: Left side quadrant setting, defaults to 1.
    :type lo: Number
    :param hi: Right side quadrant setting, defaults to 1.
    :type hi: Number
    :return: Extended data as X / Y, and lo_trim / hi_trim slicing index to access original data within the extension.
    :rtype: Tuple[np.ndarray, np.ndarray, int, int]

    .. TODO:: Make example images
    """
    # Setup number of data points based on index mirroring
    index = 0  # mirroring index (0=means mirror at first data point at low and end last data point on hi end)
    ndata = xdata.shape[0] - (index + 1)

    # Perform lo padding operations
    if lo == 2:
        lo_xdata = min(xdata) + xdata[index] - xdata[::-1][index + 1:]
        lo_ydata = ydata[::-1][index + 1:]
    elif lo == 3:
        lo_xdata = min(xdata) + xdata[index] - xdata[::-1][index + 1:]
        lo_ydata = ydata[index] - ydata[::-1][index + 1:]
    elif lo == 4:
        lo_xdata = min(xdata) + xdata[index] - xdata[::-1][index + 1:]
        lo_ydata = ydata[index + 1:] - (ydata[-(index + 1)] - ydata[index])
    else:
        lo_xdata = np.array([])
        lo_ydata = np.array([])

    # Perform hi padding operations
    if hi == 2:
        hi_xdata = -min(xdata) + max(xdata) + xdata[index + 1:]
        hi_ydata = ydata[::-1][index + 1:]
    elif hi == 3:
        hi_xdata = -min(xdata) + max(xdata) + xdata[index + 1:]
        hi_ydata = ydata[index] - ydata[::-1][index + 1:] + 2 * ydata[-(index + 1)] + ydata[index]
    elif hi == 4:
        hi_xdata = -min(xdata) + max(xdata) + xdata[index + 1:]
        hi_ydata = ydata[-(index + 1)] + ydata[index + 1:]
    else:
        hi_xdata = np.array([])
        hi_ydata = np.array([])

    # Assemble data
    if lo in [2, 3, 4] and hi in [2, 3, 4]:
        xdata = np.concatenate((lo_xdata, xdata, hi_xdata), axis=0)
        ydata = np.concatenate((lo_ydata, ydata, hi_ydata), axis=0)
        lo_trim = ndata
        hi_trim = -ndata
    elif lo in [2, 3, 4] and hi == 1:
        xdata = np.concatenate((lo_xdata, xdata), axis=0)
        ydata = np.concatenate((lo_ydata, ydata), axis=0)
        lo_trim = ndata
        hi_trim = xdata.shape[0]
    elif lo == 1 and hi in [2, 3, 4]:
        xdata = np.concatenate((xdata, hi_xdata), axis=0)
        ydata = np.concatenate((ydata, hi_ydata), axis=0)
        lo_trim = 0
        hi_trim = -ndata
    else:
        lo_trim = 0
        hi_trim = xdata.shape[0]
    return xdata, ydata, lo_trim, hi_trim

def sparse_eye_diff(m, d, format='csc'):
    diagonals = np.zeros(2*d + 1)
    diagonals[d] = 1
    for i in range(d):
        diagonals = diagonals[:-1] - diagonals[1:]
    offsets = np.arange(d+1)
    return sp.sparse.diags(diagonals, offsets, shape=(m-d, m), format=format)

def whitt_smooth(ydata, order, lambda_):
    # Base Whittaker Eilers smoother
    m = len(ydata)
    E = sp.sparse.eye(m, dtype=int, format='csc')
    D = sparse_eye_diff(m, order, format='csc')
    C = E + lambda_ * D.conj().T.dot(D)
    z = sp.sparse.linalg.spsolve(C, ydata)
    return z

"""
Functions for file resampling.

"""
from typing import Union

import numpy as np


def log_resample(
    wavelength, sampling: Union[int, float], resolution: Union[int, float]
) -> np.ndarray:
    """Resample spectrum with a given sampling per resolution element.

    Parameters
    ----------
    wavelength: np.ndarray
        Wavelength array.
    sampling: int, float
        Points to sample per resolution element
    resolution: int, float
        Instrumental resolution

    Uses faster method using log and powers of a base.
    The base is (1.0 + 1.0/(sampling*resolution).

    Almost equivalent to using
    np.logspace(np.log(wavelength)/np.log(base), np.log(wavelength)/np.log(base),
        np.log(wavelength_end / wavelength_start) / np.log(base), base).
    """
    wavelength_start = np.nanmin(wavelength)
    wavelength_end = np.nanmax(wavelength)

    # Create grid using logarithms with base of (1.0 + 1.0/(sampling*resolution))
    base = 1.0 + 1.0 / (float(sampling) * float(resolution))  # type : float
    return my_logspace(wavelength_start, wavelength_end, base, end_point=True)


def my_logspace(start, stop, base, end_point: bool = False):
    """Like np.logspace but start and stop in wavelength units.

    Parameters
    ----------
    start: float
        Starting point (in real units)
    stop: float
        End point (in real units)
    base: float
        Logarithmic base to jump between points.
    end_point: bool
        Make sure to include/go past the end point

    Returns
    -------
    logspace: ndarray
        Array of points with a spacing such that x[ii+1] = x[ii] * base
        between start and stop (or stop*base if end_point = True).

    """
    n = np.log(stop / start) / np.log(base)
    # to use the end point
    if end_point:
        n = n + 1
    powers = np.arange(np.ceil(n))
    return start * base ** powers


def log_chunks(wavelength, percent):
    """Define the bounds at which $(\Delta \lambda)/\lambda = X\%$.

    Allows spectrum to be split into chunks in which the size is X% of the given wavelength.
    This takes logarithmic steps with a base of (1+X/100).
    """
    base = 1 + percent / 100.
    wl_min = np.nanmin(wavelength)
    wl_max = np.nanmax(wavelength)

    return my_logspace(wl_min, wl_max, base, end_point=True)

"""
Auxiliary functions for eniric
"""
import collections
import errno
import os
from typing import Any, List, Optional, Sequence, Tuple, Union

import astropy.constants as const
import numpy as np
from numpy import ndarray

import eniric

c = const.c
# Band limits.
bands_ = {
    "VIS": (0.38, 0.78),
    "GAP": (0.78, 0.83),
    "Z": (0.83, 0.93),
    "Y": (1.0, 1.1),
    "J": (1.17, 1.33),
    "H": (1.5, 1.75),
    "K": (2.07, 2.35),
    "CONT": (0.45, 1.05),
    "NIR": (0.83, 2.35),
}

bands_.update(eniric.custom_bands)


def band_selector(wav: ndarray, flux: ndarray, band: str) -> Tuple[ndarray, ndarray]:
    """Select a specific wavelength band.

    Parameters
    ----------
    wav: array-like
        Wavelength values.
    flux: array-like
        Flux values.
    band: str
        Band letter to select, upper or lower case is supported. Options
        are ("ALL" or ""), "VIS", "GAP", "Z", "Y", "J", "H", "K".
    """
    band = band.upper()

    if band in ["ALL", ""]:
        return wav, flux
    else:
        band_min, band_max = band_limits(band)
        return wav_selector(wav, flux, band_min, band_max)


def band_limits(band: str) -> Tuple[float, float]:
    """ Get wavelength limits of band in microns.

    Parameters
    ----------
    band: str
        Band letter to get wavelength range for.

    Returns
    -------
    wav_min: float
        Lower wavelength bound of band in microns
    wav_max: float
        Upper wavelength bound of band in microns
    """
    if not isinstance(band, str):
        raise AttributeError("Band name must be a string")
    else:
        band = band.upper()

    if band in bands_:
        return bands_[band]
    else:
        raise ValueError("The band {0} requested is not a valid option".format(band))


def band_middle(band: str) -> float:
    """Calculate band mid-point.

    Input
    -----
    band: str
        Band label

    Return
    -------
    middle: float
        Wavelength at middle band.
    """
    band_min, band_max = band_limits(band)

    return (band_min + band_max) / 2


def wav_selector(
    wav: Union[ndarray, List[float]],
    flux: Union[ndarray, List[float]],
    wav_min: float,
    wav_max: float,
) -> Tuple[ndarray, ndarray]:
    """
    function that returns wavelength and flux within a giving range

    Parameters
    ----------
    wav: array-like
        Wavelength array.
    flux: array-like
        Flux array.
    wav_min: float
        Lower bound wavelength value.
    wav_max: float
        Upper bound wavelength value.

    Returns
    -------
    wav_sel: array
        New wavelength array within bounds wav_min, wav_max
    flux_sel: array
        New wavelength array within bounds wav_min, wav_max
        """
    wav = np.asarray(wav, dtype=float)
    flux = np.asarray(flux, dtype=float)

    mask = mask_between(wav, wav_min, wav_max)
    flux_sel = flux[mask]
    wav_sel = wav[mask]

    return wav_sel, flux_sel


def mask_between(x: ndarray, xmin: float, xmax: float) -> ndarray:
    """Create boolean mask of x between xmin and xmax."""
    return (x >= xmin) & (x < xmax)


def silent_remove(filename: str) -> None:
    """Remove file without failing when it doesn't exist."""
    try:
        os.remove(filename)
    except OSError as e:
        if e.errno != errno.ENOENT:  # errno.ENOENT = no such file or directory
            raise  # re-raise exception if a different error occurred


####################################################
def issequenceforme(obj):
    if isinstance(obj, str):
        return False
    return isinstance(obj, collections.Sequence)


def resolutions2ints(resolution: Sequence[Any]) -> List[int]:
    """List of resolutions to list of integer resolutions.

    Convert from ["100k", "30000"] to [100000, 30000].
    """
    if not issequenceforme(resolution):
        raise TypeError(
            "resolution was a {} but needs to be a Sequence".format(type(resolution))
        )

    res_list = []
    for res in resolution:
        res_list.append(res2int(res))
    return res_list


def res2int(res: Any) -> int:
    """Convert from "100k" or "100000" to 100000."""
    if issequenceforme(res):
        raise TypeError("res was a {} but needs to be a non-Sequence".format(type(res)))

    if isinstance(res, (np.int, np.float, int, float)):
        value = res
    elif isinstance(res, str):
        if res.lower().endswith("k"):
            value = float(res[:-1]) * 1000
        else:
            value = float(res)
    else:
        raise TypeError("Resolution name Type error of type {}".format(type(res)))

    return int(value)


def resolutions2strs(resolution: Sequence[Any]) -> List[str]:
    """List of resolutions to list of string resolutions.

    Convert from ["100000", 10000] to ["100k", "10k"].
    """
    if not issequenceforme(resolution):
        raise TypeError(
            "resolution was a {} but needs to be a Sequence".format(type(resolution))
        )

    res_list = []
    for res in resolution:
        res_list.append(res2str(res))
    return res_list


def res2str(res: Any) -> str:
    """Convert from "100000" or 100000 to "100k"."""
    if issequenceforme(res):
        raise TypeError(
            "resolution was a {} but needs to be a non-Sequence".format(type(res))
        )

    if isinstance(res, (np.int, np.float)):
        value = res / 1000
    elif isinstance(res, str):
        if res.lower().endswith("k"):
            value = res[:-1]
        else:
            value = float(res) / 1000
    else:
        raise TypeError("Resolution name TypeError of type {}".format(type(res)))

    return "{0}k".format(int(value))


#################################
def rv_cumulative(rv_vector: Union[List, ndarray], single: bool = False) -> List[float]:
    """Function that calculates the cumulative RV vector weighted_error."""
    if single:
        # Include 1st value for reference
        return [
            weighted_error(rv_vector[0]),
            weighted_error(rv_vector[:2]),
            weighted_error(rv_vector[:3]),
            weighted_error(rv_vector[:4]),
            weighted_error(rv_vector),
        ]

    else:
        return [
            weighted_error(rv_vector[:2]),
            weighted_error(rv_vector[:3]),
            weighted_error(rv_vector[:4]),
            weighted_error(rv_vector),
        ]


def rv_cumulative_full(rv_vector: Union[List, ndarray]) -> ndarray:
    """Function that calculates the cumulative RV vector weighted_error. In both directions."""
    assert len(rv_vector) == 5, "This hardcoded solution only meant for 5 bands."

    cumulation = np.asarray(
        [
            weighted_error(rv_vector[0]),  # First
            weighted_error(rv_vector[:2]),
            weighted_error(rv_vector[:3]),
            weighted_error(rv_vector[:4]),
            weighted_error(rv_vector[:]),  # All
            weighted_error(rv_vector[1:]),
            weighted_error(rv_vector[2:]),
            weighted_error(rv_vector[3:]),
            weighted_error(rv_vector[4]),  # Last
        ],
        dtype=float,
    )
    return cumulation


def weighted_error(rv_vector: Union[List[float], ndarray]) -> float:
    """Function that calculates the average weighted error from a vector of errors."""
    rv_vector = np.asarray(rv_vector)
    rv_value = 1.0 / (np.sqrt(np.sum((1.0 / rv_vector) ** 2.0)))

    return rv_value


def moving_average(x: ndarray, window_size: Union[int, float]) -> ndarray:
    """Moving average."""
    window = np.ones(int(window_size)) / float(window_size)
    return np.convolve(x, window, "same")


#################################
def load_aces_spectrum(
    params: Union[ndarray, List[float]],
    photons: bool = True,
    air: bool = False,
    wl_range: List[float] = [3000, 54000],
):
    """Load a Phoenix spectrum from the phoenix library using STARFISH.

    Parameters
    ----------
    params: ndarray
         [temp, logg, metallicity(, alpha)]
    photons: bool
        Necessary conversions into photons for precisions.
    air: bool
       Convert to air wavelengths (default = False).
    wl_range: [float, float]
        Min/Max wavelength range to load. Default = [3000, 54000] Angstrom.

    Returns
    -------
    wav_micron: ndarray
        Wavelength in microns
    flux_micron: ndarray
        Photon counts per (cm**2 s) or SED/micron (within a multiplicative constant 1/(h*c)).

    Spectra available from http://phoenix.astro.physik.uni-goettingen.de
    """
    base = eniric.paths["phoenix_raw"] + os.sep

    if len(params) == 3:  # Only 3 parameters given
        params = [params[0], params[1], params[2], 0]  # Set alpha=0

    if len(params) == 4:
        from Starfish.grid_tools import PHOENIXGridInterface as PHOENIX

        phoenix_grid = PHOENIX(base=base, air=air, norm=False, wl_range=wl_range)

    else:
        raise ValueError("Number of parameters is incorrect")

    wav = phoenix_grid.wl
    flux, hdr = phoenix_grid.load_flux(params)

    # Convert wavelength Angstrom to micron
    wav_micron = wav * 10 ** -4
    # Convert SED from /cm  to /micron
    flux_micron = flux * 10 ** -4

    if photons:
        # Convert to photons
        # The energy units of Phoenix fits files is erg/s/cm**2/cm
        # PHOENIX ACES gives the Spectral Energy Density (SED)
        # We transform the SED into photons by
        # multiplying the flux by the wavelength (lambda)
        #
        #     Flux_photon = Flux_energy/Energy_photon
        # with
        #     Energy_photon = h*c/lambda
        # Flux_photon = Flux_energy * lambda / (h * c)

        flux_micron = flux_micron * wav_micron
    return wav_micron, flux_micron


def load_btsettl_spectrum(
    params: Union[ndarray, List[float]],
    photons: bool = True,
    air: bool = False,
    wl_range: List[float] = [3000, 30000],
):
    """Load a BT-Settl spectrum from the CIFIST2011 library using STARFISH.

    Parameters
    ----------
    params: ndarray
         [temp, logg]. Metallicity = 0, alpha = 0
    photons: bool
        Necessary conversions into photons for precisions.
    air: bool
       Convert to air wavelengths (default = False).
    wl_range: [float, float]
        Min/Max wavelength range to load. Default = [3000, 30000] Angstrom.

    Returns
    -------
    wav_micron: ndarray
        Wavelength in microns
    flux_micron: ndarray
        Photon counts per (cm**2 s) or SED/micron. (within a multiplicative constant 1/(h*c))

    Notes
    -----
    From BT-SETTL readme:
        CIFIST2011_2015: published version of the BT-Settl grid (Baraffe et al. 2015,
        Allard et al. 2015. This grid will be the most complete
        of the CIFIST2011 grids above, but currently: Teff = 1200 - 7000K, logg = 2.5 - 5.5,
        [M/H] = 0.0.

        Available from https://phoenix.ens-lyon.fr/Grids/BT-Settl/CIFIST2011_2015/FITS/
    """
    from Starfish.grid_tools import CIFISTGridInterface as BTSETTL

    if (2 < len(params)) and (len(params) <= 4):
        assert params[2] == 0
        assert params[-1] == 0  # Checks index 3 when present.
        params = params[0:2]  # Only allow 2 params

    base = eniric.paths["btsettl_raw"] + os.sep

    btsettl_grid = BTSETTL(base=base, air=air, norm=False, wl_range=wl_range)

    wav = btsettl_grid.wl
    # Convert wavelength from Angstrom to micron
    wav_micron = wav * 10 ** -4

    # CIFIST flux is  W/m**2/um
    flux, hdr = btsettl_grid.load_flux(params)

    flux_micron = flux * 10 ** -7  # Convert W/m**2/um to ergs/s/m**2/um)
    flux_micron *= 10 ** -4  # Convert 1/m**2 to 1/cm**2

    if photons:
        # Convert to photon counts:
        # The energy units of CIFIST fits files is W/m**2/um
        # We have converted it to ergs/(s cm**2 um)
        # BT-Settl gives the Spectral Energy Density (SED)
        # We transform the SED into photons by
        # multiplying the flux by the wavelength (lambda)
        #
        #     Flux_photon = Flux_energy/Energy_photon
        # with
        #     Energy_photon = h*c/lambda
        # Flux_photon = Flux_energy * lambda / (h * c)
        flux_micron = flux_micron * wav_micron
    return wav_micron, flux_micron


#####################################################
def doppler_shift_wav(wavelength: ndarray, vel: float):
    r"""Doppler shift wavelength by a given velocity (non-relativistic).

    Apply Doppler shift to the wavelength values of the spectrum
    using the velocity value provided and the relation
    \(\Delta\lambda / \lambda = v / c\)

    Parameters
    ----------
    wavelength: ndarray
        Wavelength vector
    vel : float
        Velocity to Doppler shift by in km/s.

    Notes
    -----
    The Doppler shift is calculated using the relation
       \Delta\lambda / \lambda\] = \[v / c
    Where RV is the radial velocity (in km/s), \(\lambda_0\)`
    is the rest wavelength and \(\Delta\lambda\) is the wavelength
    shift, \(\lambda_{shift} - \lambda_0\)
    """

    if not np.isfinite(vel):
        ValueError("The velocity is not finite.")

    shifted_wavelength = wavelength * (1 + (vel * 1000 / c.value))
    return shifted_wavelength


def doppler_shift_flux(
    wavelength: ndarray, flux: ndarray, vel: float, new_wav: Optional[ndarray] = None
):
    r"""Doppler shift flux by a given velocity, return it at the original wavelengths (non-relativistic).

    Apply Doppler shift to the wavelength values of the spectrum
    using the velocity value provided and the relation
    \(\Delta\lambda / \lambda = v / c\)

    Then linearly interpolate the flux with the new wavelength to the old wavelengths.

    Parameters
    ----------
    wavelength: ndarray
        Wavelength vector
    flux: ndarray
        Flux vector
    vel : float
        Velocity to Doppler shift by in km/s.
    new_wav: Optional[ndarray]
        New wavelength array to evaluate the doppler shifted flux at.
        If None then defaults to new_wav=wavelength.
    Returns
    -------
    new_flux: ndarray
        Doppler-shifted flux evaluated at new_wav.
    """
    shifted_wavelength = doppler_shift_wav(wavelength, vel)

    if new_wav is None:
        new_wav = wavelength
    new_flux = np.interp(new_wav, shifted_wavelength, flux)
    return new_flux


def doppler_limits(rvmax, wmin, wmax):
    """Calculate wavelength limits to apply if preforming doppler shift.

    To avoid any edge effects within wmin and wmax after doppler shift.

    Parameters
    ----------
    rvmax: float
        Maximium absolute RV offset in km/s. Uses np.abs() to constrain as absolute.
    wmin: float
        Lower wavelength limit.
    wmax: float
        Upper wavelength limit.
    Returns
    -------
    new_wmin: float
      Lower wavelength bound shifted by -rvmax
    new_wmax: float
       Lower wavelength bound shifted by +rvmax
    """
    c_km = const.c.value / 1000  # c in km/s
    doppler_minus, doppler_plus = (1 - np.abs(rvmax) / c_km), (1 + np.abs(rvmax) / c_km)

    return wmin * doppler_minus, wmax * doppler_plus

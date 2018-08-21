"""
Auxiliary functions for eniric
"""
import collections
import errno
import os
from typing import Any, List, Sequence, Tuple, Union

import numpy as np
from numpy import ndarray
from Starfish.grid_tools import (
    CIFISTGridInterface as BTSETTL,
    PHOENIXGridInterface as PHOENIX,
    PHOENIXGridInterfaceNoAlpha as PHOENIXNoAlpha,
)

import eniric

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


def band_middle(band):
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
    wav = np.asarray(wav, dtype="float64")
    flux = np.asarray(flux, dtype="float64")

    mask = mask_between(wav, wav_min, wav_max)
    flux_sel = flux[mask]
    wav_sel = wav[mask]

    return wav_sel, flux_sel


def mask_between(x, xmin, xmax):
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
def load_aces_spectrum(params, photons=True, air=False):
    """Load a Phoenix spectrum from the phoenix library using STARFISH.

    Parameters
    ----------
    params: ndarray
         [temp, logg, metallicity(, alpha)]
    photons: bool
        Necessary conversions into photons for precisions.

    Returns
    -------
    wav_micron: ndarray
        Wavelength in microns
    flux_micron: ndarray
        Photon counts or SED/micron
    """
    base = eniric.paths["phoenix_raw"] + os.sep

    if params[3] == 0:  # Alpha value
        params = params[:-1]
        assert len(params) == 3
        phoenix_grid = PHOENIXNoAlpha(base=base, air=air, norm=False)
    elif len(params) == 4:
        print("USING ALPHA in PHOENIX LOADING")
        phoenix_grid = PHOENIX(
            base=base, air=air, norm=False
        )  # , param_names = ["temp", "logg", "Z", "alpha"])
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
        """The energy units of Phoenix fits files is erg/s/cm**2/cm
        PHOENIX ACES gives the Spectral Energy Density (SED)
        We transform the SED into photons by
        multiplying the flux by the wavelength (lambda)

            Flux_photon = Flux_energy/Energy_photon
        with
            Energy_photon = h*c/lambda
        Flux_photon = Flux_energy * lambda / (h * c)

        Here we convert the flux into erg/s/cm**2/\mum by multiplying by 10**-4 cm/micron
        Flux_e(erg/s/cm**2/\mum)  = Flux_e(erg/s/cm**2/cm) * (1 cm) / (10000 \mum)
        """
        # Turn into photon counts
        flux_micron = flux_micron * wav_micron
    return wav_micron, flux_micron


def load_btsettl_spectrum(params, photons=True, air=False):
    """Load a BT-Settl spectrum from the CIFIST2011 library using STARFISH.

    Parameters
    ----------
    params: ndarray
         [temp, logg]. Metallicity = 0, alpha = 0
    photons: bool
        Necessary conversions into photons for precisions.

    Returns
    -------
    wav_micron: ndarray
        Wavelength in microns
    flux_micron: ndarray
        Photon counts or SED/micron

    Notes:
    From BT-SETTL readme -
        CIFIST2011_2015: published version of the BT-Settl grid (Baraffe et al. 2015,
        Allard et al. 2015. This grid will be the most complete
        of the CIFIST2011 grids above, but currently: Teff = 1200 - 7000K, logg = 2.5 - 5.5,
        [M/H] = 0.0.
    """
    if (2 < len(params)) and (len(params) <= 4):
        assert params[2] == 0
        assert params[-1] == 0  # Checks index 3 when present.
        params = params[0:2]  # Only allow 2 params

    base = eniric.paths["btsettl"] + os.sep

    btsettl_grid = BTSETTL(base=base, air=air, norm=False, wl_range=[3000, 24000])

    wav = btsettl_grid.wl
    flux, hdr = btsettl_grid.load_flux(params)

    # Convert wavelength Angstrom to micron
    wav_micron = wav * 10 ** -4
    # Convert SED from /cm  to /micron
    # TODO check the conversion fraction.
    flux_micron = flux * 10 ** -4

    if photons:
        # Convert to photons
        """The energy units of CIFIST fits files is erg/s/cm**2/cm
        BT-Settl gives the Spectral Energy Density (SED)
        We transform the SED into photons by
        multiplying the flux by the wavelength (lambda)

            Flux_photon = Flux_energy/Energy_photon
        with
            Energy_photon = h*c/lambda
        Flux_photon = Flux_energy * lambda / (h * c)

        Here we convert the flux into erg/s/cm**2/\mum by multiplying by 10**-4 cm/micron
        Flux_e(erg/s/cm**2/\mum)  = Flux_e(erg/s/cm**2/cm) * (1 cm) / (10000 \mum)
        """
        # Turn into photon counts
        flux_micron = flux_micron * wav_micron
    return wav_micron, flux_micron

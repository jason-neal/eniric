"""
Auxiliary functions for nIRanalysis

"""
import collections
import errno
import os
import re
from typing import Any, List, Optional, Tuple, Union, Sequence

import numpy as np
from Starfish.grid_tools import PHOENIXGridInterface as PHOENIX
from Starfish.grid_tools import PHOENIXGridInterfaceNoAlpha as PHOENIXNoAlpha
from numpy import ndarray

import eniric
import eniric.IOmodule as io


def read_spectrum(spec_name: str) -> Tuple[ndarray, ndarray]:
    """Function that reads a flux spectra from the database!.

    If a energy flux spectra is read then it converts it to photons.

    Parameters
    ----------
    spec_name: str
        Location and name of model spectrum file.

    Returns
    -------
    wav: array-like, float64
        Wavelength in microns.
    flux_photons: array-like, float64
        Photon flux.

    """
    if "_res" in spec_name or "_vsini" in spec_name:
        raise ValueError("Using wrong function to load resampled spectrum.")

    if "photon" in spec_name:
        wav_micron, flux_photons = io.pdread_2col(spec_name)
    else:
        wav, flux = io.pdread_2col(spec_name)

        wav_micron = wav * 1.0e-4  # Conversion to microns
        flux_photons = flux * wav_micron  # Convert to photons

    return wav_micron, flux_photons


def get_spectrum_name(
        startype: str,
        logg: Union[float, int] = 4.50,
        feh: Union[float, int] = 0.0,
        alpha: Optional[Union[int, float]] = None,
        org: bool = False,
        flux_type: str = "photon",
) -> str:
    """Return correct phoenix spectrum filename for a given spectral type.

    Based off phoenix_utils module.

    Parameters
    ----------
    flux_type: str
        Indicate which file type to try find. e.g. "photon", "wave", ("fits" Not Implemented yet)

    Returns
    -------
    spectrum_name: str
        The name of spectrum with chosen Parameters


    Ability to select logg and metallicity (feh) later on.
    org = original locations (without Z folder option)
    """
    if feh == 0:
        feh = -0.0  # make zero negative to signed integer.

    temps = {"M0": 3900, "M3": 3500, "M6": 2800, "M9": 2600}
    if (flux_type == "photon") and (not org):
        base = "PHOENIX-ACES-AGSS-COND-2011-HiRes_wave_photon.dat"
    else:
        base = "PHOENIX-ACES-AGSS-COND-2011-HiRes_wave.dat"

    # noinspection SpellCheckingInspection
    if startype in temps.keys():
        if org:
            phoenix_name = "lte{0:05d}-{1}-{2}.{3}".format(
                temps[startype], "4.50", "0.0", base
            )
        elif (alpha is not None) and (alpha != 0.0):
            if abs(alpha) > 0.2:
                raise ValueError(
                    "Warning! Alpha is outside acceptable range -0.2->0.2. (for current science case)"
                )

            phoenix_name = os.path.join(
                "Z{0:+4.1f}.Alpha={1:+5.2f}".format(feh, alpha),
                "lte{0:05d}-{1:4.2f}{2:+4.1f}.Alpha={3:+5.2f}.{4:s}".format(
                    temps[startype], logg, feh, alpha, base
                ),
            )
        else:
            phoenix_name = os.path.join(
                "Z{0:+4.1f}".format(feh),
                "lte{0:05d}-{1:4.2f}{2:+4.1f}.{3:s}".format(
                    temps[startype], logg, feh, base
                ),
            )

        spectrum_name = phoenix_name
    elif re.match(r"^[OBAFGKML][0-9]$", startype):  # Valid spectral types
        raise NotImplementedError(
            "The spectral type '{0:s}' is not implemented yet.".format(startype)
        )
    else:
        raise ValueError("'{0:s}' is not a valid spectral type.".format(startype))

    return spectrum_name


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

    # bands = {"VIS": (0.38, 0.78), "GAP": (0.78, 0.83), "Z": (0.83, 0.93),
    #         "Y": (1.0, 1.1), "J": (1.17, 1.33), "H": (1.5, 1.75),
    #         "K": (2.07, 2.35), "CONT": (0.45, 1.05), "NIR": (0.83, 2.35)}
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

    bands = {
        "VIS": (0.38, 0.78),
        "GAP": (0.78, 0.83),
        "Z": (0.83, 0.93),
        "Y": (1.0, 1.1),
        "J": (1.17, 1.33),
        "H": (1.5, 1.75),
        "K": (2.07, 2.35),
        "CONT": (0.45, 1.05),
        "NIR": (0.83, 2.35),
        "CARMENES_NIR": (0.96, 1.71),
        "CARMENES_VIS": (0.52, 0.96),
    }

    if band in bands:
        return bands[band]
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


def unitary_gaussian(
        x: Union[range, int, ndarray],
        center: Union[float, int, str],
        fwhm: Union[float, int, str],
) -> ndarray:
    """Gaussian function of area = 1.

    Parameters
    ----------
    x: array-like
        Position array
    center: float
        Central position of Gaussian
    fwhm: float
        Full Width at Half Maximum

    Returns
    -------
    result: array-like
        Result of gaussian function sampled at x values.
    """
    if not isinstance(fwhm, (np.float, np.int)):
        raise TypeError("The fwhm value is not a number, {0}".format(type(fwhm)))
    if not isinstance(center, (np.float, np.int)):
        raise TypeError("The center value is not a number, {0}".format(type(center)))
    if not isinstance(x, np.ndarray):
        raise TypeError

    sigma = np.abs(fwhm) / (2 * np.sqrt(2 * np.log(2)))
    amp = 1.0 / (sigma * np.sqrt(2 * np.pi))
    tau = -((x - center) ** 2) / (2 * (sigma ** 2))
    result = amp * np.exp(tau)

    return result


def rotation_kernel(
        delta_lambdas: ndarray, delta_lambda_l: float, vsini: float, epsilon: float
) -> ndarray:
    """Calculate the rotation kernel for a given wavelength

    Parameters
    ----------
    delta_lambdas: array
        Wavelength values selected within delta_lambda_l around central value. (check)
    delta_lambda_l: float
        FWHM of rotational broadening. (check)
    vsini: float
        Projected rotational velocity [km/s]
    epsilon: float
        Linear limb-darkening coefficient (0-1).

    Returns
    -------
        Rotational kernel

    Notes:
    Equations * from .... book.

    """
    denominator = np.pi * vsini * (1.0 - epsilon / 3.0)
    lambda_ratio_sqr = (delta_lambdas / delta_lambda_l) ** 2.0

    c1 = 2.0 * (1.0 - epsilon) / denominator
    c2 = 0.5 * np.pi * epsilon / denominator

    return c1 * np.sqrt(1.0 - lambda_ratio_sqr) + c2 * (1.0 - lambda_ratio_sqr)


def oned_circle_kernel(x, center, fwhm):
    """Calculate the convolution kernel for a circular fiber.

    Parameters
    ----------
    x: array
        Value to evaluate kernel at.
    center: float
        Center of kernel.
    fwhm: float
        FWHM of desired kernel.

    Returns
    -------
        Collapsed circle kernel

    Notes:
    Tries to represent the broadening by the fiber of a fiber feed spectrograph.

    Artigau 2018 - stated mathematically equivalent to a cosine between -pi/2 and pi/2. This is what has tried to be created.
    """
    fwhm_scale = 2.0943951  # Numerically derived

    A = 1  # Amplitude
    B = fwhm_scale / fwhm  # Scale to give specific fwhm

    result = A * np.cos(B * (x - center))

    # Limit to main cos lobe only
    upper_xi = center + np.pi / 2 / B
    lower_xi = center - np.pi / 2 / B
    mask = mask_between(x, lower_xi, upper_xi)
    result[~mask] = 0

    return result


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
        raise TypeError("resolution was a {} but needs to be a Sequence".format(type(resolution)))

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
        raise TypeError("resolution was a {} but needs to be a Sequence".format(type(resolution)))

    res_list = []
    for res in resolution:
        res_list.append(res2str(res))
    return res_list


def res2str(res: Any) -> str:
    """Convert from "100000" or 100000 to "100k"."""
    if issequenceforme(res):
        raise TypeError("resolution was a {} but needs to be a non-Sequence".format(type(res)))

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

def load_aces_spectrum(params, photons=True):
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
        phoenix_grid = PHOENIXNoAlpha(base=base)
    elif len(params) == 4:
        print("USING ALPHA in PHOENIX LOADING")
        phoenix_grid = PHOENIX(
            base=base
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


# TODO: Use BT-Settl also
def load_btsettl_spectrum(params, photons=True):
    raise NotImplementedError("Need to include BT-SETTL")

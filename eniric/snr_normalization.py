"""Automated snr normalization.

The photon flux can be scaled by multiplicative constant which
affects the SNR of the spectra and the radial velocity precision.

For consistency and valid comparision we normalize each spectra
to achieve a consistent SNR at a specific location.
"""

from typing import Optional, Union

import numpy as np
from numpy import ndarray

import eniric.utilities as utils


def snr_constant_band(
    wav: ndarray,
    flux: ndarray,
    snr: Union[int, float] = 100,
    band: str = "J",
    sampling: Union[int, float] = 3.0,
) -> float:
    """Determine the normalization constant to achieve a SNR in the middle of a given band.

    SNR estimated by the square root of the number of photons in a resolution element.

    Parameters
    ----------
    wav: ndarray
        Wavelength array (microns)
    flux: ndarray
        Photon flux array (photons/s/cm**2)
    snr: int,  default = 100
        SNR to normalize to.
    band: str, default = "J"
        Band to use for normalization.
    sampling: int or float
       Number of pixels per resolution element.

    Returns
    -------
    normalization_value: float
        Normalization value to divide spectrum by to achieve a
        signal-to-noise level of snr within an resolution element
        in the middle of the band.

    Note
    ----
    If sampling is a float it will rounded to the nearest integer for indexing.
    """
    band_middle = utils.band_middle(band)

    if not (wav[0] < band_middle < wav[-1]):
        raise ValueError("Band center not in wavelength range.")

    norm_constant = snr_constant_wav(
        wav, flux, wav_ref=band_middle, snr=snr, sampling=sampling
    )

    return norm_constant


def snr_constant_wav(
    wav: ndarray,
    flux: ndarray,
    wav_ref: float,
    snr: Union[int, float] = 100,
    sampling: Union[int, float] = 3.0,
) -> float:
    """Determine the normalization constant to achieve a SNR at given wavelength.

    SNR estimated by the square root of the number of photons in a resolution element.

    Parameters
    ----------
    wav: ndarray
        Wavelength in micron
    flux: ndarray
        Photon flux array (photons/s/cm**2)
    wav_ref: float
        Wavelength to set the SNR per resolution element.
    snr: int,  default = 100
        SNR to set.
    sampling: int or float
       Number of pixels per resolution element.

    Returns
    -------
    norm_value: float
        Normalization value to divide the flux by to achieve the desired SNR "SNR"
        in resolution element (defined by "sampling") around the wavelength "wav_ref".

    Notes
    -----
    We want to be consistent for each spectra. If we specify the middle of J band
    as the reference it will need to be used for all bands of that spectra.

    If sampling is a float it will rounded to the nearest integer for indexing.
    """
    if wav_ref < wav[0] or wav_ref > wav[-1]:
        raise ValueError("Reference wavelength is outside of the wavelength bounds")
    index_ref = np.searchsorted(wav, [wav_ref])[0]  # Searching for the closest index

    indexes = sampling_index(
        index_ref, sampling=np.round(sampling), array_length=len(wav)
    )

    snr_estimate = np.sqrt(np.sum(flux[indexes]))

    print(
        "\tSanity Check: The reference S/N at {1:3.02f} was of {0:4.2f}.".format(
            snr_estimate, wav_ref
        )
    )
    norm_value = (snr_estimate / snr) ** 2
    return norm_value


def sampling_index(
    index: int, sampling: int = 3, array_length: Optional[int] = None
) -> ndarray:
    """Get a small number of index values around the given index value.

    Parameters
    ----------
    index: int
        The index value which to select values around.
    sampling: int
        Number of index values to return (sampling per resolution element). Must be an integer.
    array_length: int or None, default = None
        Length of array the indexes will be used in. To not exceed array length.

    Returns
    -------
    indexes: ndarray of int64
        The index values.

    """
    import math

    half_sampling = math.floor(sampling / 2)
    if sampling % 2 == 0:  # even sampling
        # index values must be integer
        indexes = np.arange(index - half_sampling, index + half_sampling, dtype=int)
        assert len(indexes) % 2 == 0  # confirm even
        assert len(indexes) == sampling
    else:
        indexes = index + np.arange(-half_sampling, sampling - half_sampling, dtype=int)
        assert len(indexes) % 2 != 0  # confirm odd
        assert len(indexes) == sampling

    if array_length is not None:
        if np.any(indexes >= array_length):
            # This may need fixed up in the future.
            raise ValueError("Indexes has values greater than the length of array.")

    if np.any(indexes < 0):
        raise ValueError("Indexes has values less than 0.")

    return indexes

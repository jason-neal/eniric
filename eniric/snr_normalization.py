"""Automated snr normalization.

The photon flux can be scaled by multiplicative constant which
affects the SNR of the spectra and the radial velocity precision.

For consistency and valid comparision we normalize each spectra
to achieve a consistent SNR at a specific location.
"""

# Normalize to SNR 100 in middle of J band 1.25 micron!
import re
import os
import numpy as np

import eniric
import eniric.IOmodule as Io
import eniric.utilities as utils

file_error_to_catch = getattr(__builtins__, 'FileNotFoundError', IOError)

resampled_dir = eniric.paths["resampled"]


def normalize_spectrum(*args, **kwargs):
    raise NotImplementedError("Use normalize_flux")

def normalize_flux2(*args, **kwargs):
    raise NotImplementedError("Use normalize_flux")


def normalize_flux(flux, id_string, new=True, snr=100, ref_band="J"):
    """Normalize flux to have SNR of 100 in middle of J band.

    Parameters
    ----------
    flux_stellar: ndarray
        Photon flux.
    id_string: str
        Idenitifing sting for spectra.
    new: bool default=True
        Choose between new and old constant for testing.
    snr: int default=100
        SNR to normalize to, .
    ref_band: str default="J"
        References band to normalize to.

    Returns
    -------
    normalized_flux: ndarray
        Flux normalized to a S/N of SNR in the middle of the ref_band.

    """
    # print("Starting norm of {}".format(id_string))
    if new:
        wav_ref, flux_ref = get_reference_spectrum(id_string, ref_band)
        norm_const = snr_constant_band(wav_ref, flux_ref, snr=snr, band=ref_band)
    else:
        if ref_band != "J" or snr != 100:
            raise ValueError("The old way does not work with these reference values.")
        norm_const = old_norm_constant(id_string) * 1e4  # Input flux offset

    print("{0:s} normalization constant = {1:f}".format(id_string, norm_const))

    return flux / norm_const


def old_norm_constant(id_string):
    """Normalization constants for Figueira et al 2016.

    These are the manual values to achieve a SNR of 100 at 1.25 micro
    for the set of parameter combinations presented.
    """
    if "M0" in id_string and ("1.0"in id_string or "5.0"in id_string or "10.0" in id_string):
        norm_constant = 1607
    elif "M3" in id_string and ("1.0"in id_string or "5.0"in id_string or "10.0" in id_string):
        norm_constant = 1373
    elif "M6" in id_string and "1.0" in id_string:
        norm_constant = 933
    elif "M6" in id_string and "5.0" in id_string:
        norm_constant = 967
    elif "M6" in id_string and "10.0" in id_string:
        norm_constant = 989
    elif "M9" in id_string and "1.0" in id_string:
         norm_constant = 810
    elif "M9" in id_string and "5.0" in id_string:
        norm_constant = 853
    elif "M9" in id_string and "10.0" in id_string:
        norm_constant = 879
    else:
        print("Constant not defined. Aborting...")
        raise ValueError("Bad ID string")

    return (norm_constant / 100.0) ** 2.0


def get_reference_spectrum(id_string, ref_band="J"):
    """From the id_string find the correct Spectrum to calculate norm_constant from."""
    # TODO: add option for Alpha into ID-String
    # TODO: Add metalicity and logg into id string
    # TODO: Add metalicity folder
    ref_band = ref_band.upper()
    # Determine the correct reference file to use.
    if ("Alpha=" in id_string) or ("smpl" in id_string):
        raise NotImplementedError
    else:
        try:
            # Use band for when ref_band is "self"
            star, band, vel, res = re.search(r"(M\d)-(\w)-(\d{1,2}\.0)-(\d{2,3}k)", id_string).groups()

            if band not in eniric.bands["all"]:
                raise ValueError
        except:
            raise ValueError("Id-string {} is not valid for normalization.".format(id_string))

        smpl = 3.0   # Fixed value atm
    if ref_band == "SELF":
        ref_band = band

    file_to_read = ("Spectrum_{0}-PHOENIX-ACES_{1}band_vsini{2}_R{3}"
                    "_res{4:1.0f}.txt").format(star, ref_band, vel, res, smpl)

    try:
        wav_ref, flux_ref = Io.pdread_2col(os.path.join(eniric.paths["resampled"], file_to_read))
    except file_error_to_catch:
        print("The reference spectra in {0:s} band was not found for id {1:s}".format(ref_band, id_string))
        raise file_error_to_catch
    return wav_ref, flux_ref


def snr_constant_band(wav, flux, snr=100, band="J"):
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

    Returns
    -------
    normalization_value: float
        Normalization value to divide spectrum by to achieve a
        signal-to-noise level of snr within an resolution element
        in the middle of the band.

    """
    band_min, band_max = utils.band_limits(band)

    band_middle = (band_min + band_max) / 2

    if not (wav[0] < band_middle < wav[-1]):
        raise ValueError("Band center not in wavelength range.")

    norm_constant = snr_constant_wav(wav, flux, wav_ref=band_middle, snr=snr)

    return norm_constant


def snr_constant_wav(wav, flux, wav_ref, snr=100, sampling=3):
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
    sampling: int
       Number of pixels per resolution element.

    Returns
    -------
    norm_value: float
        Normalization value to divide the flux by to acheive the desired SNR "SNR"
        in resolution element (defined by "sampling") around the wavelength "wav_ref".

    Notes
    -----
    We want to be consistent for each spectra. If we specify the middle of J band
    as the reference it will need to be used for all bands of that spectra.

    """
    index_ref = np.searchsorted(wav, [wav_ref])[0]  # Searching for the closest index

    indexes = sampling_index(index_ref, sampling=sampling, array_length=len(wav))

    snr_estimate = np.sqrt(np.sum(flux[indexes]))

    print("\tSanity Check: The S/N for the reference model was of {:4.2f}.".format(snr_estimate))
    norm_value = (snr_estimate / snr)**2
    return norm_value


def sampling_index(index, sampling=3, array_length=None):
    """Get a small number of index values around the given index value.

    Parameters
    ----------
    index: int
        The index value which to select values around.
    sampling: int
        Number of index values to return (sampling per resolution element)
    array_length: int or None, default = None
        Length of array the indexes will be used in. To not exceed array length.

    Returns
    -------
    indexes: ndarray of int64
        The index values.

    """
    import math
    half_sampling = math.floor(sampling / 2)
    if sampling % 2 == 0:    # even sampling
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

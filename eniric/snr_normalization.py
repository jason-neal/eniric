
""" The photon flux can be scaled by multiplicative constant which affects the
SNR of the spectra and the radial velocity precision.

For consistancy and valid comparision we normalize each spectra to acheive a
consistant SNR at a specific location.
"""

# Normaize to SNR 100 in middle of J band 1.25 micron!
import re
import numpy as np
import eniric.utilities as utils
import eniric.IOmodule as IO

file_error_to_catch = getattr(__builtins__, 'FileNotFoundError', IOError)


def normalize_flux(flux_stellar, id_string):
    """Normalize flux to have SNR of 100 in middle of J band.

    This is the original values to acheive a SNR of 100 in the middle of the J band at 1.25 micron.
    """

    if("M0" in id_string):
        norm_constant = 1607

    elif("M3" in id_string):
        norm_constant = 1373

    elif("M6" in id_string):
        if("1.0" in id_string):
            norm_constant = 933
        elif("5.0" in id_string):
            norm_constant = 967
        else:
            norm_constant = 989

    elif("M9" in id_string):
        if("1.0" in id_string):
            norm_constant = 810
        elif("5.0" in id_string):
            norm_constant = 853
        else:
            norm_constant = 879
    else:
        print("Constant not defined. Aborting...")
        exit()

    return flux_stellar / ((norm_constant / 100.0)**2.0)


def get_reference_spectrum(id_string, ref_band="J", resampled_dir="data/resampled/"):
    """ From the id_string find the correct Spectrum to
    calculate norm_constant from"""
    # TODO: add option for Alpha into ID-String
    # TODO: Add metalicity and logg into id string
    # TODO: Add metalicity folder

    # Determine the corrrect reference file to use.
    if ("Alpha=" in id_string) or ("smpl" in id_string):
        raise NotImplementedError
    else:
        try:
            # don't need the band value as will use ref_band
            star, vel, res = re.search(r"(M\d)-[ZYHJK]-(\d{1,2}\.0)-(\d{2,3}k)", id_string).groups()
        except:
            raise ValueError("Id-string {} is not valid for normalization.".format(id_string))

        smpl = 3.0   # Fixed value atm
    file_to_read = ("Spectrum_{0}-PHOENIX-ACES_{1}band_vsini{2}_R{3}"
                    "_res{4:1.0f}.txt").format(star, ref_band, vel, res, smpl)

    try:
        wav_ref, flux_ref = IO.pdread_2col(resampled_dir + file_to_read)
    except file_error_to_catch:
        print("The reference spectra in {0:s} band was not found for id {1:s}".format(ref_band, id_string))
        raise
    return wav_ref, flux_ref


def normalize_spectrum(id_string, wav, flux, snr=100, ref_band="J", resampled_dir="data/resampled/"):
    """Normalize spectrum flux with to have a SNR of snr in middle of ref_band band."""
    if ref_band in id_string:
        # Dont need to load different spectrum as already have the reference spectrum
        wav_ref, flux_ref = wav, flux
    else:
        wav_ref, flux_ref = get_reference_spectrum(id_string, ref_band, resampled_dir)

    norm_constant = snr_constant_band(wav_ref, flux_ref, snr=100, band=ref_band)

    normalized_flux = flux / norm_constant
    print("{0:s} norm_constant = {1:f}".format(id_string, norm_constant))
    return normalized_flux


def snr_constant_band(wav, flux, snr=100, band="J"):
    """ Determine the normalization constant to acheive a SNR in the middle of a given band.

    SNR estimated by the square root of the number of photons in a resolution element.

    Parameters
    ----------
    id_string: str
        Identifying string for spectrum we want to normalize.
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
        Normalization value to divide spectrum by to achive a signal-to-noise level of snr within an resolution element in the middle of the band.

    """

    band_min, band_max = utils.band_limits(band)

    band_middle = (band_min + band_max) / 2

    if band_middle < wav[0] or band_middle > wav[-1]:
        # not in range
        pass
    # Option to specify own wavelength?
    norm_constant = snr_constant_wav(wav, flux, band_middle, snr=snr)

    # Test it
    # flux2 = flux / norm_constant

    return norm_constant


def snr_constant_wav(wav, flux, wav_ref, snr=100, sampling=3):
    """ Determine the normalization constant to acheive a SNR at given wavelength.

    SNR estimated by the square root of the number of photons in a resolution element.

    Parameters
    ----------
    wav: ndarray
        Wavlength in micron
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

    index_ref = np.searchsorted(wav, wav_ref)  # Searching for the closest index

    indexes = sampling_index(index_ref, sampling=sampling, array_length=len(wav))

    snr_estimate = np.sqrt(np.sum(flux[indexes]))

    # print("\tSanity Check: The S/N for the reference model was of {:4.2f}.".format(snr_estimate))
    norm_value = (snr_estimate / snr)**2
    return norm_value


def sampling_index(index, sampling=3, array_length=None):
    """ Get a small number of index values around the given index value.

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
    if sampling % 2 == 0:    # even sampling
        # index values must be integer
        indexes = np.arange(index - sampling / 2, index + sampling / 2, dtype=int)
        assert len(indexes) % 2 == 0  # confirm even
        assert len(indexes) == sampling
    else:
        indexes = index + np.arange(-np.floor(sampling / 2), sampling - np.floor(sampling / 2), dtype=int)
        assert len(indexes) % 2 != 0  # confirm odd
        assert len(indexes) == sampling

    if array_length is not None:
        if np.any(indexes >= array_length):
            # This may need fixed up in the future.
            raise ValueError("Indexes has values greater than the length of array.")

    if np.any(indexes < 0):
        raise ValueError("Indexes has values less than 0.")

    return indexes

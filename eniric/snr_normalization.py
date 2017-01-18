
""" The photon flux can be scaled by multiplicative constant which affects the SNR of the spectra and the radial velocity precision.

For consistancy and valid comparision we normalize each spectra to acheive a consistant SNR at a specific location. """

# Normaize to SNR 100 in middle of J band 1.25 micron!
import numpy as np
from eniric.utilities import band_limits


def normalize_flux(flux_stellar, id_string):
    """Normalize flux to have SNR of 100 in middle of J band."""

    """ This is the original values to acheive a SNR of 100 in the middle of the J band at 1.25 micron """


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




def snr_constant_band(wav, flux, SNR=100, band="J"):
    """ Determine the constant to acheive a SNR in the middle of a given band."

    SNR estimated by the square root of the number of photons in a resolution element.

    Parameters
    ----------
    wav: ndarray (micron)
    flux: ndarray (photons/s/cm**2)
    SNR: int,  default = 100
    band: str, default = "J"

    Returns
    -------
    normalization_value: float

    """

    band_min, band_max = band_limits(band)

    band_middle = (band_min + band_max) / 2

    if band_middle < wav[0] or band_middle > wav[-1]:
        # not in range
        pass
    # Option to specify own wavelength?
    norm_constant = snr_constant_wav(wav, flux, band_middle, SNR=SNR)

    # Test it
    # flux2 = flux / norm_constant

    return norm_constant


def snr_constant_wav(wav, flux, wav_ref, SNR=100, sampling=3):
    """ Determine the constant to acheive a SNR at given wavelength."

    SNR estimated by the square root of the number of photons in a resolution element.

    Parameters
    ----------
    wav: ndarray
        Wavlength in micron
    flux: ndarray
        Photon flux array (photons/s/cm**2)
    wav_ref: float
        Wavelength to set the SNR per resolution element.
    SNR: int,  default = 100
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

    SN_estimate = np.sqrt(np.sum(flux[indexes]))

    # print("\tSanity Check: The S/N for the reference model was of {:4.2f}.".format(SN_estimate))
    norm_value = (SN_estimate / SNR)**2
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
        indexes = np.arange(index - sampling/2, index + sampling/2, dtype=int)
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

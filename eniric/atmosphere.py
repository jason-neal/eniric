"""Functions to deal with the atmosphere models.

Mainly the barycentric shifting of the absorption mask.
"""

from typing import Tuple

import matplotlib.pyplot as plt
import numpy as np
from astropy.constants import c
from numpy import ndarray

import eniric.IOmodule as io


def prepare_atmosphere(atmmodel: str) -> Tuple[ndarray, ndarray, ndarray, ndarray]:
    """Read in atmospheric model and prepare."""
    wav_atm, flux_atm, std_flux_atm, mask_atm = io.pdread_4col(atmmodel)
    # pandas already returns numpy arrays
    wav_atm = wav_atm / 1000.0  # conversion from nanometers to micrometers
    mask_atm = np.array(mask_atm, dtype=bool)
    return wav_atm, flux_atm, std_flux_atm, mask_atm


def barycenter_shift(wav_atm: ndarray, mask_atm: ndarray, rv_offset: float = 0.0,
                     consecutive_test: bool = True) -> ndarray:
    """Calculating impact of Barycentric movement on mask...

    Extends the masked region to +-30 km/s due to the barycentric motion of the earth.
    """
    # Mask values to the left and right side of mask_atm.
    # To avoid indexing errors have padded with first and last values.
    mask_iminus1 = np.concatenate(([mask_atm[0], mask_atm[0]], mask_atm[:-2]))  # padding with first value
    mask_iplus1 = np.concatenate((mask_atm[2:], [mask_atm[-1], mask_atm[-1]]))  # padding with last value

    pixels_total = len(mask_atm)
    masked_start = pixels_total - np.sum(mask_atm)

    barycenter_rv = 3e4  # 30 km/s in m/s
    offset_rv = rv_offset * 1.0e3  # Convert to m/s

    # Doppler shift  applied to the vectors
    delta_lambdas = wav_atm * barycenter_rv / c.value
    offset_lambdas = wav_atm * offset_rv / c.value  # offset lambda

    # Doppler shift limits of each pixel
    wav_lower_barys = wav_atm + offset_lambdas - delta_lambdas
    wav_upper_barys = wav_atm + offset_lambdas + delta_lambdas

    mask_atm_30kms = np.empty_like(mask_atm, dtype=bool)

    for i, (wav_value, mask_val) in enumerate(zip(wav_atm, mask_atm)):
        """If there are 3 consecutive zeros within +/-30km/s then make the value 0."""

        # Offset_RV is the offset applied for the star RV.
        if (((mask_val is False) and (mask_iminus1[i] is False) and
             (mask_iplus1[i] is False) and (rv_offset == 0))):  # If the mask is false and the offset is equal to zero
            """If the value and its 2 neighbours are already zero don't do the barycenter shifts"""
            mask_atm_30kms[i] = False
        else:
            # np.searchsorted is faster then the boolean masking wavelength range
            # It returns index locations to place the min/max doppler-shifted values
            slice_limits = np.searchsorted(wav_atm, [wav_lower_barys[i], wav_upper_barys[i]])
            slice_limits = [index if (index < len(wav_atm)) else len(wav_atm) - 1
                            for index in slice_limits]  # Fix searchsorted index

            mask_atm_slice = mask_atm[slice_limits[0]:slice_limits[1]]

            mask_atm_slice = np.asarray(mask_atm_slice, dtype=bool)  # Insure boolean dtype

            if consecutive_test:
                # Make mask value False if there are 3 or more consecutive zeros in slice.
                len_consec_zeros = consecutive_truths(~mask_atm_slice)
                if np.all(~mask_atm_slice):  # All pixels of slice is zeros (shouldn't get here)
                    mask_atm_30kms[i] = False
                elif np.max(len_consec_zeros) >= 3:
                    mask_atm_30kms[i] = False
                else:
                    mask_atm_30kms[i] = True
                    if np.sum(~mask_atm_slice) > 3:
                        print("There were {0}/{1} zeros in this barycentric shift but None were 3 consecutive!".format(
                            np.sum(~mask_atm_slice), len(mask_atm_slice)))
            else:
                mask_atm_30kms[i] = np.product(mask_atm_slice)
    masked_end = pixels_total - np.sum(mask_atm_30kms)
    print(("New Barycentric impact affects the number of masked pixels by {0:04.1%} due to the atmospheric"
           " spectrum").format((masked_end - masked_start) / pixels_total))
    print(("Masked Pixels start = {1}, masked_pixel_end = {0}, Total = {2}").format(masked_end, masked_start,
                                                                                    pixels_total))
    return mask_atm_30kms


def consecutive_truths(condition: ndarray) -> ndarray:
    """Length of consecutive true values in an bool array.

    Parameters
    ----------
    condition: ndarray of bool
        True False array of a condition.

    Returns
    -------
    len_consecutive: ndarray of ints
        Array of lengths of consecutive true values of the condition.

    Notes
    -----
    Solution found at {http://stackoverflow.com/questions/24342047/
                       count-consecutive-occurences-of-values-varying-in-length-in-a-numpy-array}
    """
    if not np.any(condition):  # No match to condition
        return np.array([0])
    else:
        unequal_consec = np.concatenate(([condition[0]], condition[:-1] != condition[1:], [True]))
        where_changes = np.where(unequal_consec)[0]  # indices where condition changes
        len_consecutive = np.diff(where_changes)[::2]  # step through every second to get the "True" lenghts.
    return len_consecutive


def plot_atm_masks(wav, flux, mask, old_mask=None, new_mask=None, block=True):
    """Plot spectrum with absorption masks."""
    plt.figure()
    plt.plot(wav, flux, label="Spectrum")
    plt.plot(wav, mask, label="Absorption Mask")

    if old_mask is not None:
        plt.plot(wav, old_mask + 0.01, ".-", label="Old code Mask")

    if new_mask is not None:
        plt.plot(wav, new_mask + 0.02, "+-", label="New code Mask")

    plt.legend(loc=0)
    plt.show(block=block)

    return 0


def bugged_old_barycenter_shift(wav_atm: ndarray, mask_atm: ndarray, rv_offset: float = 0.0):
    """Buggy Old version Calculating impact of Barycentric movement on mask...

    Extends the masked region to +-30 km/s due to the barycentric motion of the earth.
    """
    pixels_total = len(mask_atm)
    masked_start = pixels_total - np.sum(mask_atm)

    mask_atm_30kms = []
    for value in zip(wav_atm, mask_atm):
        if (value[1] is False) and (rv_offset == 666.0):  # if the mask is false and the offset is equal to zero
            mask_atm_30kms.append(value[1])

        else:
            delta_lambda = value[0] * 3.0e4 / c.value
            starting_lambda = value[0] * rv_offset * 1.0e3 / c.value
            indexes_30kmslice = np.searchsorted(wav_atm, [starting_lambda + value[0] - delta_lambda,
                                                          starting_lambda + value[0] + delta_lambda])
            indexes_30kmslice = [index if (index < len(wav_atm)) else len(wav_atm) - 1 for index in indexes_30kmslice]

            mask_atm_30kmslice = np.array(mask_atm[indexes_30kmslice[0]:indexes_30kmslice[1]],
                                          dtype=bool)  # selecting only the slice in question

            mask_atm_30kmslice_reversed = [not i for i in mask_atm_30kmslice]

            # Found suspected culprit code
            clump = np.array_split(mask_atm_30kmslice, np.where(np.diff(mask_atm_30kmslice_reversed))[0] + 1)[
                    ::2]  # This code is the bug!

            tester = True
            for block in clump:
                if len(clump) >= 3:  # This is the bug! should read len(block)
                    tester = False
                    break

            mask_atm_30kms.append(tester)

    mask_atm = np.array(mask_atm_30kms, dtype=bool)
    masked_end = pixels_total - np.sum(mask_atm)
    print(("Old Barycentric impact affects number of masked pixels by {0:04.1%} due to the atmospheric"
           " spectrum").format((masked_end - masked_start) / pixels_total))
    print(("Pedros Pixels start = {1}, Pixel_end = {0}, Total = {2}").format(masked_end, masked_start, pixels_total))
    return mask_atm


def old_barycenter_shift(wav_atm: ndarray, mask_atm: ndarray, rv_offset: float = 0.0) -> ndarray:
    """Bug Fixed old version Calculating impact of Barycentric movement on mask...

    Extends the masked region to +-30 km/s due to the barycentric motion of the earth.

    This fixed version is used to compare with the new version.
    """
    pixels_total = len(mask_atm)
    masked_start = pixels_total - np.sum(mask_atm)

    mask_atm_30kms = []
    for value in zip(wav_atm, mask_atm):
        if (value[1] is False) and (rv_offset == 666.0):  # if the mask is false and the offset is equal to zero
            mask_atm_30kms.append(value[1])

        else:
            delta_lambda = value[0] * 3.0e4 / c.value
            starting_lambda = value[0] * rv_offset * 1.0e3 / c.value
            indexes_30kmslice = np.searchsorted(wav_atm, [starting_lambda + value[0] - delta_lambda,
                                                          starting_lambda + value[0] + delta_lambda])
            indexes_30kmslice = [index if (index < len(wav_atm)) else len(wav_atm) - 1 for index in indexes_30kmslice]

            mask_atm_30kmslice = np.array(mask_atm[indexes_30kmslice[0]:indexes_30kmslice[1]],
                                          dtype=bool)  # selecting only the slice in question

            # mask_atm_30kmslice_reversed = [not i for i in mask_atm_30kmslice] # Unneeded

            # Found suspected culprit code
            # clump = np.array_split(mask_atm_30kmslice,
            #                        np.where(np.diff(mask_atm_30kmslice_reversed))[0] + 1)[::2]  # <<<<<< Bug!
            if mask_atm_30kmslice[0]:  # A true first value
                clump = np.array_split(mask_atm_30kmslice, np.where(np.diff(mask_atm_30kmslice))[0] + 1)[1::2]
            else:
                clump = np.array_split(mask_atm_30kmslice, np.where(np.diff(mask_atm_30kmslice))[0] + 1)[::2]

            tester = True
            for block in clump:
                if True in block:
                    raise ValueError("There is a true value in False blocks!")

                if len(block) >= 3:  # This is the bug! should read len(block)
                    tester = False
                    break

            mask_atm_30kms.append(tester)

    mask_atm = np.array(mask_atm_30kms, dtype=bool)
    masked_end = pixels_total - np.sum(mask_atm)
    print(("Old Barycentric impact affects number of masked pixels by {0:04.1%} due to the atmospheric"
           " spectrum").format((masked_end - masked_start) / pixels_total))
    print(("Pedros Pixels start = {1}, Pixel_end = {0}, Total = {2}").format(masked_end, masked_start, pixels_total))
    return mask_atm


def dopplershift(wav, flux):
    """Doppler shift the flux of spectrum."""
    pass
    # return newflux

def atm_mask(flux, cutoff=0.98):
    """Mask flux below the cutoff value."""
    if not isinstance(flux, np.ndarray):
        flux = np.asarary(flux)
    return flux < cutoff

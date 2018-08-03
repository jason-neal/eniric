"""Functions to deal with the atmosphere models.

Mainly the barycentric shifting of the absorption mask.
"""

from typing import List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
from astropy.constants import c
from numpy import ndarray

import eniric.IOmodule as io


class Atmosphere(object):
    """Atmospheric transmission object.

    Stores wavelength and atmospheric transmission arrays.
    """

    def __init__(self, wavelength, transmission, std=None, mask=None):
        assert len(wavelength) == len(
            transmission
        ), "Wavelength and transmission do not match length."
        self.wl = np.asarray(wavelength)
        self.transmission = np.asarray(transmission)
        if std is None:
            self.std = np.zeros_like(wavelength)
        else:
            self.std = np.asarray(std)
        if mask is None:
            self.mask = np.ones_like(wavelength, dtype=bool)
        else:
            self.mask = np.asarray(mask, dtype=bool)
        self.shifted = False

    @classmethod
    def from_file(cls, atmmodel: str):
        """Read in atmospheric model and prepare.

        Alternate constructor for Atmosphere.

        Parameters
        ----------
        atmmodel: str
            Name of atmosphere file.
        """
        wav_atm, flux_atm, std_flux_atm, mask_atm = io.pdread_4col(atmmodel)
        wav_atm = wav_atm / 1e3  # conversion from nanometers to micrometers
        mask_atm = np.array(mask_atm, dtype=bool)
        # We do not use the std from the year atmosphere.
        return cls(
            wavelength=wav_atm, transmission=flux_atm, std=std_flux_atm, mask=mask_atm
        )

    def to_file(
        self, new_atmmodel: str, header: Optional[List[str]] = None, fmt: str = "%11.8f"
    ):
        """Save the atmospheric model to new_atmmodel file.

        Converts micron back into nanometers to be consistent with from_file().

        Parameters
        ----------
        new_atmmodel: str
            Name of atmosphere file to save to.
        header:
            Header lines to add.
        fmt: str
             String formatting
        """
        if header is None:
            header = ["# atm_wav(nm)", "atm_flux", "atm_std_flux", "atm_mask"]
        io.pdwrite_cols(
            new_atmmodel,
            self.wl * 1000,
            self.transmission,
            self.std,
            self.mask,
            header=header,
            float_format=fmt,
        )

    def mask_transmission(self, depth: float) -> None:
        """Mask the transmission below given depth. e.g. 3%

        Parameters
        ----------
        depth : float
            Telluric line depth percentage to mask out.
            E.g. depth=2 will mask transmission deeper than 2%.

        Returns
        -------

        """
        cutoff = 1 - depth / 100.0
        self.mask = self.transmission < cutoff

    def bary_shift_mask(self, rv: float = 30.0, consecutive_test: bool = False):
        """RV shift mask symmetrically.

        Parameters
        ----------
        rv: float (default=30 km/s)
            Barycentric RV to extend masks in km/s. (Default=30 km/s)
        consecutive_test: bool (default False)
            Checks for 3 consecutive zeros to mask out transmission.

        """
        rv_mps = rv * 1e3  # Convert from km/s into m/s

        shift_amplitudes = self.wl * rv_mps / c.value
        # Operate element wise
        blue_shifts = self.wl - shift_amplitudes
        red_shifts = self.wl + shift_amplitudes

        bary_mask = []
        for (blue_wl, wl, red_wl, mask) in zip(
            blue_shifts, self.wl, red_shifts, self.mask
        ):
            if mask == 0:
                this_mask_value = False
            else:
                # np.searchsorted is faster then the boolean masking wavelength range
                # It returns index locations to place the min/max doppler-shifted values
                slice_limits = np.searchsorted(self.wl, [blue_wl, red_wl])
                slice_limits = [
                    index if (index < len(self.wl)) else len(self.wl) - 1
                    for index in slice_limits
                ]  # Fix searchsorted end index

                mask_slice = self.mask[slice_limits[0] : slice_limits[1]]

                if consecutive_test:
                    # Make mask value False if there are 3 or more consecutive zeros in slice.
                    len_consec_zeros = consecutive_truths(~mask_slice)
                    if np.all(
                        ~mask_slice
                    ):  # All pixels of slice is zeros (shouldn't get here)
                        this_mask_value = False
                    elif np.max(len_consec_zeros) >= 3:
                        this_mask_value = False
                    else:
                        this_mask_value = True
                        if np.sum(~mask_slice) > 3:
                            print(
                                "There were {0}/{1} zeros in this barycentric shift but None were 3 consecutive!".format(
                                    np.sum(~mask_slice), len(mask_slice)
                                )
                            )
                else:
                    this_mask_value = np.bool(
                        np.product(mask_slice)
                    )  # Any 0s will make it 0

                # Checks

                if not this_mask_value:
                    assert np.any(mask_slice == False)
                else:
                    if not consecutive_test:
                        assert np.all(mask_slice)
            bary_mask.append(this_mask_value)
        # print("bary mask", bary_mask)
        self.mask = np.asarray(bary_mask, dtype=np.bool)

    def broaden(self, resolution: float):
        """Broaden atmospheric transmission profile.

        This does not change any created masks.

        Parameters
        ----------
        resolution: float
            Instrumental resolution/resolving power
        """
        from eniric.broaden import resolution_convolution

        self.transmission = resolution_convolution(
            self.wl, self.transmission, R=resolution
        )


def prepare_atmosphere(atmmodel: str) -> Tuple[ndarray, ndarray, ndarray, ndarray]:
    """Read in atmospheric model and prepare."""
    wav_atm, flux_atm, std_flux_atm, mask_atm = io.pdread_4col(atmmodel)
    # pandas already returns numpy arrays
    wav_atm = wav_atm / 1000.0  # conversion from nanometers to micrometers
    mask_atm = np.array(mask_atm, dtype=bool)
    return wav_atm, flux_atm, std_flux_atm, mask_atm


def barycenter_shift(
    wav_atm: ndarray,
    mask_atm: ndarray,
    rv_offset: float = 0.0,
    consecutive_test: bool = True,
) -> ndarray:
    """Calculating impact of Barycentric movement on mask...

    Extends the masked region to +-30 km/s due to the barycentric motion of the earth.
    """
    # Mask values to the left and right side of mask_atm.
    # To avoid indexing errors have padded with first and last values.
    mask_iminus1 = np.concatenate(
        ([mask_atm[0], mask_atm[0]], mask_atm[:-2])
    )  # padding with first value
    mask_iplus1 = np.concatenate(
        (mask_atm[2:], [mask_atm[-1], mask_atm[-1]])
    )  # padding with last value

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
        if (
            (mask_val is False)
            and (mask_iminus1[i] is False)
            and (mask_iplus1[i] is False)
            and (rv_offset == 0)
        ):  # If the mask is false and the offset is equal to zero
            """If the value and its 2 neighbours are already zero don't do the barycenter shifts"""
            mask_atm_30kms[i] = False
        else:
            # np.searchsorted is faster then the boolean masking wavelength range
            # It returns index locations to place the min/max doppler-shifted values
            slice_limits = np.searchsorted(
                wav_atm, [wav_lower_barys[i], wav_upper_barys[i]]
            )
            slice_limits = [
                index if (index < len(wav_atm)) else len(wav_atm) - 1
                for index in slice_limits
            ]  # Fix searchsorted index

            mask_atm_slice = mask_atm[slice_limits[0] : slice_limits[1]]

            mask_atm_slice = np.asarray(
                mask_atm_slice, dtype=bool
            )  # Insure boolean dtype

            if consecutive_test:
                # Make mask value False if there are 3 or more consecutive zeros in slice.
                len_consec_zeros = consecutive_truths(~mask_atm_slice)
                if np.all(
                    ~mask_atm_slice
                ):  # All pixels of slice is zeros (shouldn't get here)
                    mask_atm_30kms[i] = False
                elif np.max(len_consec_zeros) >= 3:
                    mask_atm_30kms[i] = False
                else:
                    mask_atm_30kms[i] = True
                    if np.sum(~mask_atm_slice) > 3:
                        print(
                            "There were {0}/{1} zeros in this barycentric shift but None were 3 consecutive!".format(
                                np.sum(~mask_atm_slice), len(mask_atm_slice)
                            )
                        )
            else:
                mask_atm_30kms[i] = np.product(mask_atm_slice)
    masked_end = pixels_total - np.sum(mask_atm_30kms)
    print(
        (
            "New Barycentric impact affects the number of masked pixels by {0:04.1%} due to the atmospheric"
            " spectrum"
        ).format((masked_end - masked_start) / pixels_total)
    )
    print(
        ("Masked Pixels start = {1}, masked_pixel_end = {0}, Total = {2}").format(
            masked_end, masked_start, pixels_total
        )
    )
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
        unequal_consec = np.concatenate(
            ([condition[0]], condition[:-1] != condition[1:], [True])
        )
        where_changes = np.where(unequal_consec)[0]  # indices where condition changes
        len_consecutive = np.diff(where_changes)[
            ::2
        ]  # step through every second to get the "True" lenghts.
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

import os
import re
from genericpath import isfile
from os.path import join
from typing import Optional, Tuple, Union

import numpy as np
from astropy import constants as const
from numpy.core.multiarray import ndarray

import eniric
import eniric.obsolete.IOmodule
from eniric import io_module as io
from eniric.atmosphere import consecutive_truths
from eniric.resample import log_resample


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
    delta_lambdas = wav_atm * barycenter_rv / const.c.value
    offset_lambdas = wav_atm * offset_rv / const.c.value  # offset lambda

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


def read_spectrum(spec_name: str) -> Tuple[ndarray, ndarray]:
    """Function that reads a flux spectra from the database!.

    If a energy flux spectra is read then it converts it to photons.

    Parameters
    ----------
    spec_name: str
        Location and name of model spectrum file.

    Returns
    -------
    wav: array-like, float
        Wavelength in microns.
    flux_photons: array-like, float
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


results_dir = eniric.paths["results"]
resampled_dir = eniric.paths["resampled"]


def resample_allfiles(
    results_dir: Optional[str] = None, resampled_dir: Optional[str] = None
) -> int:
    """Resample all files inside results_dir folder.

    Parameters
    ----------
    results_dir: str
        Directory containing results to resample.
    resampled_dir: str
        Directory to save resampled results.
    """
    if results_dir is None:
        results_dir = eniric.paths["results"]  # type: ignore
    if resampled_dir is None:
        resampled_dir = eniric.paths["resampled"]  # type: ignore
    # Getting a list of all the files
    onlyfiles = [f for f in os.listdir(results_dir) if isfile(join(results_dir, f))]

    for spectrum_file in onlyfiles:
        if spectrum_file.endswith(".txt"):
            resampler(
                spectrum_file, results_dir=results_dir, resampled_dir=resampled_dir
            )

    return 0


def resampler(
    spectrum_name: str = "Spectrum_M0-PHOENIX-ACES_Yband_vsini1.0_R60k.txt",
    results_dir: str = results_dir,
    resampled_dir: str = resampled_dir,
    sampling: Union[int, float] = 3.0,
) -> int:
    """Resample a spectrum file by interpolation onto a grid with a
    sampling of 3 pixels per resolution element.

    Parameters
    ----------
    spectrum_name: str
        Name of spectrum.
    results_dir: str
        Directory to find the spectrum to load.
    resample_dir: str
        Directory to save the results.
    sampling: float (default=3.0)
        Sampling per pixel.
    """
    os.makedirs(resampled_dir, exist_ok=True)
    read_name = os.path.join(results_dir, spectrum_name)

    match = re.search("_R(\d{2,3})k", spectrum_name)
    if match:
        resolution = int(match.group(1)) * 1000  # type: int
    else:
        raise Exception("Did not match Resolution")

    wavelength, __, spectrum = io.pdread_3col(read_name, noheader=True)
    wav_grid = log_resample(wavelength, sampling, resolution)

    interpolated_flux = np.interp(wav_grid, wavelength, spectrum)

    output_path = [
        resampled_dir,
        "{0}_res{1:3.01f}.txt".format(spectrum_name[:-4], float(sampling)),
    ]
    filetowrite = os.path.join(*output_path)
    eniric.obsolete.IOmodule.write_e_2col(
        filetowrite, wav_grid[1:-2], interpolated_flux[1:-2]
    )  # [1:-2] for border effects

    return 0

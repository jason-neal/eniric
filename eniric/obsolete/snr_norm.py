import os
import re
from typing import Tuple, Union

from numpy.core.multiarray import ndarray

import eniric
from eniric.io_module import pdread_2col
from eniric.snr_normalization import snr_constant_band

resampled_dir = eniric.paths["resampled"]


def normalize_flux(
    flux: ndarray,
    id_string: str,
    new: bool = True,
    snr: Union[int, float] = 100.0,
    ref_band: str = "J",
    sampling: Union[int, float] = 3.0,
) -> ndarray:
    """Normalize flux to have SNR of 100 in middle of reference band.

    Parameters
    ----------
    flux: ndarray
        Photon flux.
    id_string: str
        Identifying string for spectra.
    new: bool default=True
        Choose between new and old constant for testing.
    snr: int, float default=100
        SNR to normalize to, .
    ref_band: str default="J"
        References band to normalize to.
    sampling: int or float
       Number of pixels per resolution element.

    Returns
    -------
    normalized_flux: ndarray
        Flux normalized to a S/N of SNR in the middle of the ref_band.

    """
    # print("Starting norm of {}".format(id_string))
    if new:
        if ref_band.upper() == "SELF":
            __, ref_band, __, __ = decompose_id_string(id_string)

        wav_ref, flux_ref = get_reference_spectrum(id_string, ref_band)
        norm_const = snr_constant_band(
            wav_ref, flux_ref, snr=snr, band=ref_band, sampling=sampling
        )
    else:
        if (ref_band.upper() != "J") or (snr != 100):
            raise ValueError("The old way does not work with these reference values.")
        norm_const = old_norm_constant(id_string) * 1e4  # Input flux offset

    print("{0:s} normalization constant = {1:f}".format(id_string, norm_const))

    return flux / norm_const


def old_norm_constant(id_string: str) -> float:
    """Normalization constants for Figueira et al 2016.

    These are the manual values to achieve a SNR of 100 at 1.25 micro
    for the set of parameter combinations presented.
    """
    if "M0" in id_string and (
        "1.0" in id_string or "5.0" in id_string or "10.0" in id_string
    ):
        norm_constant = 1607
    elif "M3" in id_string and (
        "1.0" in id_string or "5.0" in id_string or "10.0" in id_string
    ):
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


def decompose_id_string(id_string: str) -> Tuple[str, str, str, str]:
    """Get the values back out of the id-string."""
    match = re.search(r"(M\d)-(\w{1,4})-(\d{1,2}\.0)-(\d{2,3}k)", id_string)
    if match:
        star, band, vel, res = match.groups()
    else:
        raise ValueError(
            "Id-string {0} is not valid for normalization.".format(id_string)
        )

    if band not in eniric.bands["all"]:
        raise ValueError(
            "The band '{}' found in id-string is not valid.".format(id_string)
        )

    return star, band, vel, res


def get_reference_spectrum(
    id_string: str, ref_band: str = "J"
) -> Tuple[ndarray, ndarray]:
    """From the id_string find the correct Spectrum to calculate norm_constant from."""
    # TODO: add option for Alpha into ID-String
    # TODO: Add metallicity and logg into id string
    # TODO: Add metallicity folder

    # Determine the correct reference file to use.
    if ("Alpha=" in id_string) or ("smpl" in id_string):
        raise NotImplementedError
    else:
        star, band, vel, res = decompose_id_string(id_string)

    smpl = 3.0  # Fixed value atm

    ref_band = ref_band.upper()
    if ref_band == "SELF":
        ref_band = band

    file_to_read = (
        "Spectrum_{0}-PHOENIX-ACES_{1}band_vsini{2}_R{3}" "_res{4:3.01f}.txt"
    ).format(star, ref_band, vel, res, float(smpl))

    try:
        wav_ref, flux_ref = pdread_2col(
            os.path.join(eniric.paths["resampled"], file_to_read)
        )
    except FileNotFoundError as e:
        print(
            "The reference spectra in {0:s} band was not found for id {1:s}".format(
                ref_band, id_string
            )
        )
        raise e
    return wav_ref, flux_ref

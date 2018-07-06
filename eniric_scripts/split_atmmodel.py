#!/usr/bin/env python
"""Split the large atmospheric model transmission spectra into the separate bands.
To be able to include the separate files and to speed up performances for
calculations on individual bands only.
"""

import argparse
import os
import sys
from typing import List, Optional

import numpy as np
from astropy.constants import c

import eniric
import eniric.IOmodule as io
import eniric.utilities as utils


def _parser():
    """Take care of all the argparse stuff.

    :returns: the args
    """

    parser = argparse.ArgumentParser(description="Band separate out atmospheric model.")

    parser.add_argument(
        "-m", "--model", help="Model name", type=str, default="Average_TAPAS_2014.txt"
    )
    parser.add_argument(
        "-b",
        "--bands",
        type=str,
        default=None,
        nargs="+",
        choices=["ALL", "VIS", "GAP", "Z", "Y", "J", "H", "K"],
        help="Wavelength band to select, Default='All'",
    )
    parser.add_argument(
        "-d",
        "--data_dir",
        help="Telluric model data directory",
        type=str,
        default=eniric.paths["atmmodel"],
    )
    parser.add_argument(
        "--new_name",
        default=None,
        type=str,
        help="Base name for new files. Default is the original model name.",
    )
    parser.add_argument(
        "--rv_extend",
        default=100,
        type=check_positive,
        help="Doppler RV (km/s) to extend the wavelength limits of the band. Default=100 km/s",
    )

    return parser.parse_args()


def check_positive(value: str) -> float:
    """Function to check if input is positive.

    http://stackoverflow.com/questions/14117415/in-python-using-argparse-allow-only-positive-integers.

    Parameters
    ----------
    value: "str"
        A input string from argparse to check if it is a positive number.

    Returns
    -------
    float_value: float
        The value if it is positive as a float.

    Raises
    ------
    ArgumentTypeError:
        If value is not a positive number.
    """
    if not isinstance(value, str):
        raise ValueError("Input value is not a string.")

    float_value = float(value)
    if float_value <= 0:
        raise argparse.ArgumentTypeError(
            "{0:s} is an invalid positive value".format(value)
        )
    return float_value


def main(
    model: str = "Average_TAPAS_2014.txt",
    bands: Optional[List[str]] = None,
    new_name=None,
    data_dir=None,
    rv_extend: float = 100,
):
    """Split the large atmospheric model transmission spectra into the separate bands.

    Keeps wavelength of atmosphere model as nanometers.
    """
    if bands is None:
        bands = ["All"]
    if new_name is None:
        new_name = model.split(".")[0]
    if data_dir is None:
        data_dir = eniric.paths["atmmodel"]

    model_name = os.path.join(data_dir, model)

    atm_wav, atm_flux, atm_std_flux, atm_mask = io.pdread_4col(model_name)

    # Return value from saving each band
    write_status = np.empty_like(bands, dtype=int)

    for i, band in enumerate(bands):
        band_name = "{0}_{1}.txt".format(new_name, band)
        band_min, band_max = utils.band_limits(band)

        # Doppler shift values to extend saved wavelengths
        band_min = band_min * (1 - rv_extend / c.value)
        band_max = band_max * (1 + rv_extend / c.value)

        # Convert band limits (micron) into nanometers (Keeps datafiles cleaner)
        band_min, band_max = band_min * 1e3, band_max * 1e3

        band_wav, band_flux = utils.wav_selector(atm_wav, atm_flux, band_min, band_max)
        __, band_std_flux = utils.wav_selector(
            atm_wav, atm_std_flux, band_min, band_max
        )
        __, band_mask = utils.wav_selector(atm_wav, atm_mask, band_min, band_max)
        assert (
            (len(band_wav) == len(band_flux))
            & (len(band_std_flux) == len(band_mask))
            & (len(band_flux) == len(band_mask))
        )  # Check lengths are the same

        band_mask = np.asarrya(band_mask, dtype=bool)

        # Save the result to file
        filename = os.path.join(data_dir, band_name)
        header = ["# atm_wav(nm)", "atm_flux", "atm_std_flux", "atm_mask"]

        write_status[i] = io.pdwrite_cols(
            filename,
            band_wav,
            band_flux,
            band_std_flux,
            band_mask,
            sep="\t",
            header=header,
            float_format="%10.8f",
        )

        return np.sum(write_status)  # If any extracts fail they will turn up here.


if __name__ == "__main__":
    args = vars(_parser())
    opts = {k: args[k] for k in args}
    sys.exit(main(**opts))

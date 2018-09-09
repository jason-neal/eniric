#!/usr/bin/env python
"""Script to split the large atmospheric model transmission spectra into the separate bands.
This create smaller files to load for each band for individual bands only.
"""
import argparse
import os
import sys
from os.path import join
from typing import List, Optional

import numpy as np
from astropy import constants as const

import eniric
from eniric.atmosphere import Atmosphere
from eniric.utilities import band_limits, doppler_shift_wav

atmmodel = "{0}.txt".format(eniric.atmmodel["base"])
choices = ["ALL"]
choices.extend(eniric.bands["all"])


def _parser():
    """Take care of all the argparse stuff.

    :returns: the args
    """

    parser = argparse.ArgumentParser(description="Band separate out atmospheric model.")

    parser.add_argument("-m", "--model", help="Model name", type=str, default=atmmodel)
    parser.add_argument(
        "-b",
        "--bands",
        type=str,
        default="ALL",
        nargs="+",
        choices=choices,
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
    parser.add_argument(
        "-c",
        "--cutoff_depth",
        default=2,
        type=float,
        help=r"Telluric line depth cutoff. Default = 2 precent.",
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
    model: str = atmmodel,
    bands: Optional[List[str]] = None,
    new_name: Optional[str] = None,
    data_dir: Optional[str] = None,
    rv_extend: float = 100,
    cutoff_depth: float = 2.0,
):
    """Split the large atmospheric model transmission spectra into the separate bands.

    Keeps wavelength of atmosphere model as nanometers.

    Parameters
    ----------
    model: str
        Telluric model file to load. It has columns wavelength, flux, std_flux, mask.
    bands: list[str]
        List bands to split model into separate files.
    new_name: str
        New file name base.
    data_dir: Optional[str]
        Directory for results. Can also be given in config.yaml "paths:atmmodel:"...
    rv_extend: float (positive) (default 100)
        Rv amount to extend wavelength range of telluric band. To later apply barycenter shifting.
    cutoff_depth: float
       Telluric line depth cutoff. Default = 2%.
    """
    if (bands is None) or ("ALL" in bands):
        bands_ = eniric.bands["all"]
    else:
        bands_ = bands

    if new_name is None:
        new_name = model.split(".")[0]
    if data_dir is None:
        data_dir_ = eniric.paths["atmmodel"]
    else:
        data_dir_ = str(data_dir)

    model_name = join(data_dir_, model)

    # If trying to obtain the provided model extract from and it doesn't yet exist
    # extract from tar.gz file. (Extracted it is 230 MB which is to large for Git)
    if "Average_TAPAS_2014.txt" == atmmodel:
        if not os.path.exists(model_name):
            print("Unpacking Average_TAPAS_2014.txt.tar.gz...")
            import tarfile

            with tarfile.open(str(model_name) + ".tar.gz", "r") as tar:
                tar.extractall(data_dir_)
            print("Unpacked")
    print("Loading from_file {0}".format(model_name))
    atm = Atmosphere.from_file(model_name)

    # Return value from saving each band
    write_status = np.empty_like(bands_, dtype=int)

    for i, band in enumerate(bands_):
        print("Starting {0}".format(band))
        filename_band = "{0}_{1}.txt".format(new_name, band)
        band_min, band_max = band_limits(band)

        # * 1000 to convert into km/s
        band_min = doppler_shift_wav(band_min, -rv_extend)
        band_max = doppler_shift_wav(band_max, rv_extend)

        split_atm = atm.wave_select(band_min, band_max)

        # Apply telluric line mask
        atm.mask_transmission(depth=cutoff_depth)

        # Save the result to file
        filename = join(data_dir_, filename_band)
        header = ["# atm_wav(nm)", "atm_flux", "atm_std_flux", "atm_mask"]
        print("Saving to_file {}".format(filename))
        write_status[i] = split_atm.to_file(filename, header=header, fmt="%11.8f")
    print("Done Splitting")

    return np.sum(write_status)  # If any extracts fail they will turn up here.


if __name__ == "__main__":
    args = vars(_parser())
    opts = {k: args[k] for k in args}
    exit_status = int(main(**opts))
    print("exit_status", exit_status)
    sys.exit(exit_status)

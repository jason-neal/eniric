#!/usr/bin/env python
"""Prepare_data.py

Code to take all phoenix-aces fits files and create .dat files with wavelength
and flux.
Adds them to the data directory of eniric for convolutions etc.

This wastes a lot of memory duplicating wavelength vector.

Jason Neal January 2017
"""
import argparse
import os
import re
import sys

import numpy as np
from astropy.io import fits

import eniric
import eniric.IOmodule as io


def _parser():
    """Take care of all the argparse stuff.

    :returns: the args
    """
    parser = argparse.ArgumentParser(
        description="Transform spectra to prepare for convolution."
    )
    parser.add_argument(
        "-s",
        "--startype",
        help='Spectral Type e.g "MO"',
        default=["M0"],
        type=str,
        nargs="+",
    )
    parser.add_argument(
        "-t",
        "--temp",
        help="Temperature of stars to prepare",
        type=float,
        nargs="+",
        default=[3900.0],
        choices=list(np.arange(2300, 7000, 100.0))
        + list(np.arange(7000, 12001, 200.0)),
    )
    parser.add_argument(
        "-l",
        "--logg",
        help="Logg for stellar models.",
        default=[4.50],
        type=float,
        nargs="+",
        choices=np.arange(0, 6.01, 0.5),
    )
    parser.add_argument(
        "-m",
        "--metallicity",
        type=float,
        default=[0.0],
        help="Metallicity values.",
        nargs="+",
    )
    # choices=[list(np.arange(-4.0, -2.0, 1))+list(np.arange(-2.0, 1.01, 0.5))]
    parser.add_argument(
        "-a",
        "--alpha",
        type=float,
        default=[0.0],
        choices=np.arange(-0.2, 1.201, 0.2),
        help="Alpha values. Default = [0.0]",
        nargs="+",
    )
    parser.add_argument(
        "-f",
        "--flux_type",
        type=str,
        default="photon",
        choices=["photon", "energy"],
        help="Type of flux to use. Default converts it to photons.",
    )
    parser.add_argument(
        "-d",
        "--data_dir",
        help="Data directory to save results.",
        type=str,
        default=None,
    )
    parser.add_argument(
        "-p",
        "--phoenix_dir",
        default=None,
        type=str,
        help="Phoenix directory to find fits files",
    )
    parser.add_argument(
        "-r",
        "--replace",
        action="store_true",
        help="Replace data files if already created.",
    )

    return parser.parse_args()


def main(
    startype,
    temp,
    logg,
    metallicity,
    alpha,
    flux_type="photon",
    data_dir=None,
    phoenix_dir=None,
    replace=False,
):
    """Prepare datafiles for phoenix models that match the input parameters.

    This add the wavelength information to each spectra and converts
    to microns/photons if the flux_type="photons" is given.
    We do realise that this is a waste of space and it would be more
    storage efficient to just read in the phoenix raw fits files and
    wavelength file when needed.

    """
    if data_dir is None:
        data_dir = eniric.paths["phoenix_dat"]
    os.makedirs(data_dir, exist_ok=True)

    if phoenix_dir is None:
        phoenix_dir = eniric.paths["phoenix_raw"]

    # Get Phoenix wavelength data
    wavelength_file = "WAVE_PHOENIX-ACES-AGSS-COND-2011.fits"
    wavelength = fits.getdata(os.path.join(phoenix_dir, wavelength_file))

    if flux_type == "photon":
        file_suffix = "_wave_photon.dat"  # For saving output
    else:
        file_suffix = "_wave.dat"

    # Convert Stellar_types into
    stellar_dict = {"M0": 3900.0, "M3": 3500.0, "M6": 2800.0, "M9": 2600.0}
    # Add temperature of stellar_type to temp list
    for star in startype:
        try:
            temp.append(stellar_dict[star])
        except KeyError:
            print(
                "Stellar type {0} is not implemented here (yet), submit and issue.".format(
                    star
                )
            )

    # Get all phoenix fits files we want to convert
    for (path, dirs, files) in os.walk(phoenix_dir):

        phoenix_files = []
        for f in files:
            # Test if filename meets conditions
            end_cond = f.endswith("PHOENIX-ACES-AGSS-COND-2011-HiRes.fits")

            try:
                if "Alpha=" in f:
                    (match_temp, match_logg, match_feh, match_alpha) = re.search(
                        r"(\d{5})-(\d\.\d\d)([+\-]\d\.\d)\.Alpha=([+\-]\d\.\d\d)\.", f
                    ).groups()
                    alpha_cond = float(match_alpha) in alpha
                else:
                    (match_temp, match_logg, match_feh) = re.search(
                        r"(\d{5})-(\d\.\d\d)([+\-]\d\.\d)", f
                    ).groups()
                    alpha_cond = True  # To make work
            except AttributeError:
                # This file doesn't match what we want so continue with next loop
                continue

            temp_cond = float(match_temp) in temp
            feh_cond = float(match_feh) in metallicity
            logg_cond = float(match_logg) in logg

            if np.all(
                [end_cond, temp_cond, feh_cond, logg_cond, alpha_cond]
            ):  # All conditions met
                # Matching file found
                phoenix_files.append(f)
            else:
                pass

        for phoenix_file in phoenix_files:
            z_folder = path.split(os.sep)[-1]
            os.makedirs(
                os.path.join(data_dir, z_folder), exist_ok=True
            )
            output_filename = os.path.join(
                data_dir, z_folder, phoenix_file[:-5] + file_suffix
            )  # Name of .dat file
            if os.path.exists(output_filename) and not replace:
                print(
                    "Skipping as {0} already exists (use -r to replace)".format(
                        output_filename
                    )
                )
                continue
            spectra = fits.getdata(os.path.join(path, phoenix_file))

            # Need to add conversions pedro preformed to flux!
            """The energy units of Phoenix fits files is erg/s/cm**2/cm
            We transform the flux into photons in the read_spectrum()
            function by multiplying the flux result by the wavelength (lambda)

                Flux_photon = Flux_energy/Energy_photon
            with
                Energy_photon = h*c/lambda
            Flux_photon = Flux_energy * lambda / (h * c)

            Here we convert the flux into erg/s/cm**2/\mum by multiplying by 10**-4 cm/\mum
            Flux_e(erg/s/cm**2/\mum)  = Flux_e(erg/s/cm**2/cm) * (1 cm) / (10000 \mum)
            """

            spectra_micron = spectra * 10 ** -4  # Convert   /cm    to  /micron

            if flux_type == "photon":
                wavelength_micron = (
                    wavelength * 10 ** -4
                )  # Convert Angstrom to   micron

                spectra_photon = (
                    spectra_micron * wavelength_micron
                )  # Ignoring constants h*c in photon energy equation

                result = io.pdwrite_cols(
                    output_filename,
                    wavelength_micron,
                    spectra_photon,
                    header=["# Wavelength (micron)", r"Flux (photon/s/cm^2)"],
                    float_format="%.7f",
                )

            else:
                result = io.pdwrite_cols(
                    output_filename,
                    wavelength,
                    spectra_micron,
                    header=["# Wavelength (Angstrom)", r"Flux (erg/s/cm^2/micron)"],
                    float_format=None,
                )

            if not result:
                print("Successfully wrote to ", output_filename)
            else:
                print("Failed to write to ", output_filename)

    print("Done")
    return 0


if __name__ == "__main__":
    args = vars(_parser())
    opts = {k: args[k] for k in args}
    sys.exit(main(**opts))

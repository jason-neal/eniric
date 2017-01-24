
""" Prepare_data.py

Code to take all phoenix-aces fits files and create .dat files with wavelength
and flux.
Adds them to the data directory of eniric for convolutions etc.

This wastes alot of memory duplicating wavelemgth vector.

Jason Neal Janurary 2017
"""
from __future__ import division, print_function
import re
import os
import sys
import argparse
import pandas as pd
from astropy.io import fits
from eniric.IOmodule import pdwrite_2col
from eniric.utilities import wav_selector


def _parser():
    """Take care of all the argparse stuff.

    :returns: the args
    """
    parser = argparse.ArgumentParser(description='Helpful discription')
    parser.add_argument("-s", '--startype', help='Spectral Type e.g "MO"', default=["M0"],
                        type=str, nargs="+")
    parser.add_argument("-t", "--temp", help="Temperature of stars to prepare",
                        type=float, nargs="+", default=[3900.0], choices=list(np.arange(2300, 7000, 100.0))+list(np.arange(7000, 12001, 200.0)))
    parser.add_argument("-l", "--logg", help="Logg for stellar models.", default=[4.50],
                        type=float, nargs="+", choices=np.arange(0, 6.01, 0.5))
    parser.add_argument("-m", "--metalicity", type=float, default=[0.0],
                        help="Metalicity values.", nargs="+")
                        # choices=[list(np.arange(-4.0, -2.0, 1))+list(np.arange(-2.0, 1.01, 0.5))]
    parser.add_argument("-a", "--alpha", type=float, default=[0.0],
                        choices=np.arange(-0.2,1.201, 0.2),
                        help="Metalicity values.", nargs="+")
    parser.add_argument("-f", "--flux_type", type=str, default="photon",
                        choices=["photon", "energy"],
                        help="Type of flux to use. Default converts it to photons.")
    parser.add_argument('-d', '--data_dir', help='Data directory to save results.',
                        type=str, default=None)
    parser.add_argument('-p', '--phoenix_dir', default=None, type=str,
                        help='Phoenix directory to find fits files')

    args = parser.parse_args()
    return args


def main(startype, temp, logg, metalicity, alpha, flux_type="photon", data_dir=None, phoenix_dir=None):
    if data_dir is None:
        data_dir = "../data/PHOENIX-ACES_spectra/"
    if phoenix_dir is None:
        phoenix_dir = "../../../data/fullphoenix/"

    # Get Phoenix wavelength data
    wavelength_file = "WAVE_PHOENIX-ACES-AGSS-COND-2011.fits"
    wavelength = fits.getdata(data_dir + wavelength_file)

    if flux_type == "photon":
        file_suffix = "_wave_photon.dat"
    else:
        file_suffix = "_wave.dat"

    # Convert Stellar_types into
    stellar_dict = {"M0": 3900.0, "M3": 3500.0, "M6": 2800.0, "M9": 2600.0}
    # Add temperature of stellar_type to temp list
    for star in startype:
        print(star)
        try:
            temp.append(stellar_dict[star])
        except KeyError:
            print("Stellar type {0} is not implemented here (yet), submit and issue.".format(star))

    # get all phoenix fits files we want to convert
    for (path, dirs, files) in os.walk(phoenix_dir):
        # print(path)
        # print(dirs)
        phoenix_files = [f for f in files if (
                         f.endswith("PHOENIX-ACES-AGSS-COND-2011-HiRes.fits") and (re.search("03900-4.50-0.0", f) is not None))]

        for phoenix_file in phoenix_files:
            if int(phoenix_file[3:8]) > 4000:
                # Not doing temperatures over 4000 K.
                pass
            else:
                Z_folder = path.split("/")[-1]
                os.makedirs(os.path.join(data_dir, Z_folder), exist_ok=True)  # make folder if doesn't exit
                output_filename = os.path.join(data_dir, Z_folder,
                                               phoenix_file[:-5] + file_suffix)  # Name of .dat file

                spectra = fits.getdata(os.path.join(path, phoenix_file))

                # Need to add conversions pedro preformed to flux!
                """ The energy units of Phoenix fits files is erg/s/cm**2/cm
                We transform the flux into photons in the read_spectrum()
                function by multiplying the flux result by the wavelength (lambda)

                    Flux_photon = Flux_energy/Energy_photon
                with
                    Energy_photon = h*c/lambda
                Flux_photon = Flux_energy * lambda / (h * c)

                Here we convert the flux into erg/s/cm**2/\mum by multiplying by 10**-4 cm/\mum
                Flux_e(erg/s/cm**2/\mum)  = Flux_e(erg/s/cm**2/cm) * (1 cm) / (10000 \mum)
                """

                spectra_micron = spectra * 10**-4              # Convert   /cm    to  /micron

                if flux_type == "photon":
                    wavelength_micron = wavelength * 10**-4    # Convert Angstrom to   micron

                    spectra_photon = spectra_micron * wavelength_micron  # Ignoring constants h*c in photon energy equation

                    result = pdwrite_2col(output_filename, wavelength_micron, spectra_photon,
                                          header=["# Wavelength (micron)", r"Flux (photon/s/cm^2)"])

                else:
                    result = pdwrite_2col(output_filename, wavelength, spectra_micron,
                                          header=["# Wavelength (Angstom)", r"Flux (erg/s/cm^2/micron)"])

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

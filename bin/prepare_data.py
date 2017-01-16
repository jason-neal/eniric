
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

data_dir = "../data/PHOENIX-ACES_spectra/"
phoenix_dir = "../../../data/fullphoenix/"

wavelength_file = "WAVE_PHOENIX-ACES-AGSS-COND-2011.fits"


def _parser():
    """Take care of all the argparse stuff.

    :returns: the args
    """
    parser = argparse.ArgumentParser(description='Helpful discription')
    parser.add_argument("-s", '--startype', help='Spectral Type e.g "MO"', type=str, nargs="+")
    parser.add_argument("-v", "--vsini", help="Rotational velocity of source",
                        type=float, nargs="+")
    parser.add_argument("-R", "--resolution", help="Observational resolution",
                        type=float, nargs="+")
    parser.add_argument("-b", "--band", type=str, default=["ALL"],
                        choices=["ALL", "VIS", "GAP", "Z", "Y", "J", "H", "K"],
                        help="Wavelength band to select", nargs="+")
    parser.add_argument('-d', '--data_dir', help='Data directory', type=str, default=None)
    parser.add_argument('--sample_rate', default=[3.0], type=float, nargs="+",
                        help="Resample rate, pixels per FWHM. Default=3.0")
    parser.add_argument('--results', default=None, type=str,
                        help='Result directory Default=data_dir+"/results/"')
    parser.add_argument('--resamples', default=None, type=str,
                        help='Resample directory. Default=data_dir+"/resampled/"')
    parser.add_argument('--noresample', help='Resample output', default=False,
                        action="store_true")
    parser.add_argument('--normalize', help='Normalize for wavelength step', default=True,
                        action="store_false")
    parser.add_argument('--org', help='Only use original .dat files, (temporary option)',
                        default=False, action="store_true")
    args = parser.parse_args()
    return args


def main(startype, vsini, resolution, band, data_dir=None, results=None,
         resamples=None, sample_rate=3.0, noresample=False, normalize=True,
         org=False):
    wavelength = fits.getdata(data_dir + wavelength_file)  # Phoenix wavelength

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
                                               phoenix_file[:-5] + "_wave.dat")  # Name of .dat file

                spectra = fits.getdata(os.path.join(path, phoenix_file))

                # Need to add conversions pedro preformed to flux!

                if not pdwrite_2col(output_filename, wavelength, spectra,
                                header=["# Wavelength", "Flux"]):
                    print("Successfully wrote to ", output_filename)
                else:
                    print("Failed to write to ", output_filename)

    print("Done")
    return 0


if __name__ == "__main__":
    args = vars(_parser())
    opts = {k: args[k] for k in args}
    sys.exit(main(**opts))


""" Prepare_data.py

Code to take all phoenix-aces fits files and create .dat files with wavelength
and flux.
Adds them to the data directory of eniric for convolutions etc.

This wastes alot of memory duplicating wavelemgth vector.

Jason Neal Janurary 2017
"""
from __future__ import division, print_function
import os
import sys
import pandas as pd
from astropy.io import fits
from eniric.IOmodule import pdwrite_2col

data_dir = "../data/PHOENIX-ACES_spectra/"
phoenix_dir = "../../../data/fullphoenix/"

wavelength_file = "WAVE_PHOENIX-ACES-AGSS-COND-2011.fits"


def main(flux_type="photon"):
    wavelength = fits.getdata(data_dir + wavelength_file)  # Phoenix wavelength

    if flux_type == "photon":
        file_suffix = "_wave_photon.dat"
    else:
        file_suffix = "_wave.dat"

    # get all phoenix fits files we want to convert
    for (path, dirs, files) in os.walk(phoenix_dir):
        # print(path)
        # print(dirs)
        phoenix_files = [f for f in files if
                         f.endswith("PHOENIX-ACES-AGSS-COND-2011-HiRes.fits")]

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
    sys.exit(main())

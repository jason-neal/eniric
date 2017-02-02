"""
Functions for file resampling.

"""

import re
import numpy as np
import matplotlib.pyplot as plt
from os import listdir
from os.path import isfile, join
from eniric.IOmodule import pdread_3col, write_e_2col

results_dir = "../data/results/"
resampled_dir = "../data/resampled/"


def resample_allfiles(results_dir=results_dir, resampled_dir=resampled_dir):
    """Resample all files inside folder.

    Parameters
    ----------
    results_dir: str
        Directory containing results to resample.
    """
    # getting a list of all the files
    onlyfiles = [f for f in listdir(results_dir) if isfile(join(results_dir, f))]

    [resampler(spectrum_file, results_dir=results_dir,
               resampled_dir=resampled_dir) for spectrum_file in onlyfiles
     if spectrum_file.endswith(".txt")]

    return 0

def resampler(spectrum_name="Spectrum_M0-PHOENIX-ACES_Yband_vsini1.0_R60k.txt",
              results_dir=results_dir, resampled_dir=resampled_dir,
              sampling=3.0, plottest=False):
    """
    resamples a spectrum by interpolation onto a grid with a
    sampling of 3 pixels per resolution element.
    """
    # wavelength, theoretical_spectrum, spectrum = read_3col(spectrum_name)
    read_name = results_dir + spectrum_name
    # theoretical_spectrum = data["model"].values
    wavelength, __, spectrum = pdread_3col(read_name, noheader=True)

    wavelength_start = wavelength[1]  # because of border effects
    wavelength_end = wavelength[-2]   # because of border effects

    # match = re.match("Spectrum_(M\d)-PHOENIX-ACES_([A-Z]{1,4})band_vsini
    # (\d{1,2}.?\d?)_R(\d{2,3})k(_conv_normalize)?.txt", spectrum_name)
    match = re.search("_R(\d{2,3})k", spectrum_name)
    resolution = int(match.group(1)) * 1000
    # wav_grid = [wavelength_start]
    # while(wav_grid[-1] < wavelength_end):
    #     wav_grid.append(wav_grid[-1]*(1.0+1.0/(sampling*resolution)))
    # wav_grid = np.array(wav_grid)

    # Create grid using logarithms with base of (1.0+1.0/(sampling*resolution))
    base = 1.0 + 1.0 / (sampling * resolution)
    n = np.log(wavelength_end / wavelength_start) / np.log(base)
    powers = np.arange(np.ceil(n))
    wav_grid = wavelength_start * base ** powers

    interpolated_flux = np.interp(wav_grid, wavelength, spectrum)
    filetowrite = "{0}{1}_res{2}.txt".format(resampled_dir, spectrum_name[:-4],
                                             int(sampling))
    write_e_2col(filetowrite, wav_grid, interpolated_flux)

    if(plottest):
        plt.figure(1)
        plt.xlabel(r"wavelength [$\mu$m])")
        plt.ylabel(r"flux [counts] ")
        plt.plot(wavelength, spectrum, color='k', linestyle="-",
                 label="Original spectrum")
        plt.plot(wav_grid, interpolated_flux, color='b', linestyle="-",
                 label="Interpolated spectrum")
        plt.legend(loc='best')
        plt.show()

        plt.close()

    return 0

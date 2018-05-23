"""
Functions for file resampling.

"""

import os
import re
from os.path import isfile, join
from typing import Union

import matplotlib.pyplot as plt
import numpy as np

import eniric
import eniric.IOmodule as io

results_dir = eniric.paths["results"]
resampled_dir = eniric.paths["resampled"]


def resample_allfiles(results_dir: str = None, resampled_dir: str = None) -> int:
    """Resample all files inside results_dir folder.

    Input
    -----
    results_dir: str
        Directory containing results to resample.
    resampled_dir: str
        Directory to save resampled results.
    """
    if results_dir is None:
        results_dir = eniric.paths["results"]
    if resampled_dir is None:
        resampled_dir = eniric.paths["resampled"]
    # Getting a list of all the files
    onlyfiles = [f for f in os.listdir(results_dir) if isfile(join(results_dir, f))]

    [resampler(spectrum_file, results_dir=results_dir,
               resampled_dir=resampled_dir) for spectrum_file in onlyfiles
     if spectrum_file.endswith(".txt")]

    return 0


def resampler(spectrum_name: str = "Spectrum_M0-PHOENIX-ACES_Yband_vsini1.0_R60k.txt",
              results_dir: str = results_dir, resampled_dir: str = resampled_dir,
              sampling: Union[int, float] = 3.0, plottest: bool = False) -> int:
    """Resamples a spectrum file by interpolation onto a grid with a
    sampling of 3 pixels per resolution element.

    Inputs
    ------
    spectrum_name: str
        Name of spectrum.
    results_dir: str
        Directory to find the spectrum to load.
    resample_dir: str
        Directory to save the results.
    sampling: float (default=3.0)
        Sampling per pixel.
    plottest: bool
        Plot a test of resampling.
    """
    os.makedirs(resampled_dir, exist_ok=True)
    read_name = os.path.join(results_dir, spectrum_name)

    match = re.search("_R(\d{2,3})k", spectrum_name)
    resolution = int(match.group(1)) * 1000

    wavelength, __, spectrum = io.pdread_3col(read_name, noheader=True)
    wav_grid = log_resample(wavelength, sampling, resolution)

    interpolated_flux = np.interp(wav_grid, wavelength, spectrum)

    output_path = [resampled_dir, "{0}_res{1}.txt".format(spectrum_name[:-4], int(sampling))]
    filetowrite = os.path.join(*output_path)
    io.write_e_2col(filetowrite, wav_grid[1:-2], interpolated_flux[1:-2])  # [1:-2] for border effects

    if plottest:
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


def old_resample(wavelength, sampling: Union[int, float], resolution: Union[int, float]) -> np.ndarray:
    """Re-sample spectrum with a given sampling per resolution element.

    Inputs
    ------
    wavelength: np.ndarray
        Wavelength array.
    sampling: int, float
        Points to sample per resolution element
    resolution: int, float
        Instrumental resolution

    """
    wavelength_start = wavelength[0]  # because of border effects
    wavelength_end = wavelength[-1]  # because of border effects

    wav_grid = [wavelength_start]
    while (wav_grid[-1] < wavelength_end):
        wav_grid.append(wav_grid[-1] * (1.0 + 1.0 / (sampling * resolution)))
    wav_grid = np.array(wav_grid)
    return wav_grid


def log_resample(wavelength, sampling: Union[int, float], resolution: Union[int, float]) -> np.ndarray:
    """Re-sample spectrum with a given sampling per resolution element.

    Inputs
    ------
    wavelength: np.ndarray
        Wavelength array.
    sampling: int, float
        Points to sample per resolution element
    resolution: int, float
        Instrumental resolution

    Uses faster method using log and powers of a base.
    The base is (1.0 + 1.0/(sampling*resolution).
    """
    wavelength_start = wavelength[0]  # because of border effects
    wavelength_end = wavelength[-1]  # because of border effects

    # Create grid using logarithms with base of (1.0 + 1.0/(sampling*resolution))
    base = 1.0 + 1.0 / (sampling * resolution)
    n = np.log(wavelength_end / wavelength_start) / np.log(base)
    powers = np.arange(np.ceil(n + 1))
    wav_grid = wavelength_start * base ** powers
    return wav_grid

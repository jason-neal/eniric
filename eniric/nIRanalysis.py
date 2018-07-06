#!/usr/bin/env python
"""
Created on Sun Dec 14 15:43:13 2014

@author: pfigueira

Adapted December 2016 by Jason Neal
"""
import os
import sys
from typing import Optional

import matplotlib.pyplot as plt
import multiprocess as mprocess
import numpy as np
from joblib import Memory
from matplotlib import rc
from numpy import ndarray
from tqdm import tqdm

import eniric
import eniric.IOmodule as io
from eniric.utilities import (band_selector, mask_between, read_spectrum,
                              rotation_kernel, unitary_gaussian, wav_selector)

# Cache convolution results.
cachedir = os.path.join(os.path.expanduser("~"), ".joblib")
memory = Memory(cachedir=cachedir, verbose=0)

# set stuff for latex usage
rc('text', usetex=True)

data_rep = eniric.paths["phoenix_dat"]
results_dir = eniric.paths["results"]
resampled_dir = eniric.paths["resampled"]

# models form PHOENIX-ACES
M0_ACES = data_rep + "lte03900-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes_wave.dat"
M3_ACES = data_rep + "lte03500-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes_wave.dat"
M6_ACES = data_rep + "lte02800-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes_wave.dat"
M9_ACES = data_rep + "lte02600-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes_wave.dat"


def run_convolutions(spectrum_string: str, band: str) -> None:
    """
    Runs the convolutions for a set of spectra in batch
    """
    vsini = [1.0, 5.0, 10.0]
    R = [60000, 80000, 100000]

    spectrum = spectrum_string  # removed exec usage
    print(spectrum)
    print("Running the convolutions for spectra of {0:s} in band {1:s}\n.".format(spectrum, band))
    for vel in vsini:
        for res in R:
            convolve_spectra(spectrum, band, vel, res, plot=False)


def save_convolution_results(filename: str, wavelength: ndarray, flux: ndarray,
                             convolved_flux: ndarray) -> int:
    """Saves convolution results to a file.

    Parameters
    ----------
    filename: str
    wavelength: array-like
     flux, convolved_flux
    """
    print("Saving results...")

    # Note: difference in sampling at 1.0 and 1.5 microns makes jumps
    # in the beginning of Y and H bands
    io.write_e_3col(filename, wavelength, flux, convolved_flux)
    print("Done.")
    return 0


def convolve_spectra(spectrum, band, vsini, R, epsilon: float = 0.6, fwhm_lim: float = 5.0,
                     plot: bool = True, num_procs: Optional[int] = None,
                     results_dir: str = results_dir,
                     normalize: bool = True, output_name: Optional[str] = None) -> int:
    """Load Spectrum, apply convolution and then save results.

    """
    print("Reading the data...")
    wav, flux = read_spectrum(spectrum)  # In microns and  photon flux.
    print("Done.")

    wav_band, flux_band, convolved_flux = convolution(wav, flux, vsini, R, band,
                                                      epsilon=epsilon, fwhm_lim=fwhm_lim,
                                                      num_procs=num_procs, normalize=normalize)
    if not normalize:
        norm_ = "_unnormalized"
    else:
        norm_ = ""

    if output_name is None:
        name_model = name_assignment(spectrum)

        filename = "{0}Spectrum_{1}_{2}band_vsini{3:3.1f}_R{4:d}k{5}.txt".format(
            results_dir, name_model, band, vsini, R / 1000, norm_)
    else:
        filename = os.path.join(results_dir, output_name)

    save_convolution_results(filename, wav_band, flux_band, convolved_flux)

    if plot:
        fig = plt.figure(1)
        plt.xlabel(r"wavelength [$\mu$m])")
        plt.ylabel(r"flux [counts] ")
        plt.plot(wav_band, flux_band / np.max(flux_band), color='k', linestyle="-",
                 label="Original spectra")
        plt.plot(wav_band, convolved_flux / np.max(convolved_flux), color='b', linestyle="-",
                 label="{0:s} spectrum observed at vsini={1:.2f} and R={2:d} .".format(name_model,
                                                                                       vsini, R))
        plt.legend(loc='best')
        plt.show()

        fig.savefig(filename[:-3] + "pdf", facecolor='w', format='pdf', bbox_inches='tight')
        plt.close()

    return 0


@memory.cache
def convolution(wav, flux, vsini, R, band: str = "All", epsilon: float = 0.6, fwhm_lim: float = 5.0,
                num_procs: Optional[int] = None, normalize: bool = True):
    """Perform convolution of spectrum.

    Rotational convolution followed by a Gaussian of a specified resolution R.

    Parameters
    ----------
    wav: ndarray
        Wavelength in microns
    flux: ndarray
        Photon flux
    vsini: float
        Rotational velocity in km/s.
    R: int
        Resolution of instrumental profile.
    band: str
        Wavelength band to choose, default="All"
    num_procs: int, None
        Number of processes to use with multiprocess. If None it is assigned to 1 less then total number of cores.
        If num_procs = 0, then multiprocess is not used.

    Returns
    -------
    wav_band: ndarray
        Wavelength for the selected band.
    flux_band: ndarray
        Original flux for the selected band.
    flux_conv: ndarray
        Convolved flux for the selected band.
    """

    wav_band, flux_band = band_selector(wav, flux, band)

    # We need to calculate the fwhm at this value in order to set the starting point for the convolution
    fwhm_min = wav_band[0] / R  # fwhm at the extremes of vector
    fwhm_max = wav_band[-1] / R

    # performing convolution with rotation kernel
    print("Starting the Rotation convolution for vsini={0:.2f}...".format(vsini))

    delta_lambda_min = wav_band[0] * vsini / 3.0e5
    delta_lambda_max = wav_band[-1] * vsini / 3.0e5

    # widest wavelength bin for the rotation convolution
    lower_lim = wav_band[0] - delta_lambda_min - fwhm_lim * fwhm_min
    upper_lim = wav_band[-1] + delta_lambda_max + fwhm_lim * fwhm_max
    wav_ext_rotation, flux_ext_rotation = wav_selector(wav, flux, lower_lim, upper_lim)

    # wide wavelength bin for the resolution_convolution
    lower_lim = wav_band[0] - fwhm_lim * fwhm_min
    upper_lim = wav_band[-1] + fwhm_lim * fwhm_max
    wav_extended, flux_extended = wav_selector(wav, flux, lower_lim, upper_lim)

    # rotational convolution
    flux_conv_rot = rotational_convolution(wav_extended, wav_ext_rotation,
                                           flux_ext_rotation, vsini, epsilon,
                                           num_procs=num_procs,
                                           normalize=normalize)

    print("Starting the Resolution convolution...")

    flux_conv_res = resolution_convolution(wav_band, wav_extended,
                                           flux_conv_rot, R, fwhm_lim,
                                           num_procs=num_procs,
                                           normalize=normalize)

    return wav_band, flux_band, flux_conv_res


@memory.cache
def rotational_convolution(wav_extended, wav_ext_rotation, flux_ext_rotation,
                           vsini, epsilon, num_procs=None, normalize: bool = True):
    """Perform Rotational convolution part of convolution.
    """

    def wrapper_rot_parallel_convolution(args):
        """Wrapper for rot_parallel_convolution needed to unpack the arguments for
        fast_convolve as multiprocess.Pool.map does not accept multiple
        arguments
        """
        return element_rot_convolution(*args)

    def element_rot_convolution(wav, wav_extended, wav_ext_rotation,
                                flux_ext_rotation, vsini: float, epsilon: float,
                                normalize: bool):
        """Embarrassingly parallel part of rotational convolution"""
        # select all values such that they are within the fwhm limits
        delta_lambda_l = wav * vsini / 3.0e5

        index_mask = mask_between(wav_ext_rotation, wav - delta_lambda_l, wav + delta_lambda_l)

        flux_2convolve = flux_ext_rotation[index_mask]
        rotation_profile = rotation_kernel(wav_ext_rotation[index_mask] - wav, delta_lambda_l,
                                           vsini, epsilon)

        sum_val = np.sum(rotation_profile * flux_2convolve)

        if normalize:
            # Correct for the effect of non-equidistant sampling
            unitary_rot_val = np.sum(rotation_profile)  # Affects precision
            return sum_val / unitary_rot_val
        else:
            return sum_val

    if num_procs != 0:
        if num_procs is None:
            num_procs = mprocess.cpu_count() - 1

        mproc_pool = mprocess.Pool(processes=num_procs)

        args_generator = tqdm([[wav, wav_extended, wav_ext_rotation,
                                flux_ext_rotation, vsini, epsilon, normalize]
                               for wav in wav_extended])

        flux_conv_rot = np.array(mproc_pool.map(wrapper_rot_parallel_convolution,
                                                args_generator))

        mproc_pool.close()

    else:  # num_procs == 0
        flux_conv_rot = np.empty_like(wav_extended)  # Memory assignment
        for ii, wav in enumerate(tqdm(wav_extended)):
            flux_conv_rot[ii] = element_rot_convolution(wav, wav_extended,
                                                        wav_ext_rotation,
                                                        flux_ext_rotation,
                                                        vsini, epsilon,
                                                        normalize=normalize)
        print("Done.\n")
    return flux_conv_rot


@memory.cache
def resolution_convolution(wav_band, wav_extended, flux_conv_rot, R, fwhm_lim,
                           num_procs: int = 1, normalize: bool = True):
    """Perform Resolution convolution part of convolution."""

    # Define inner convolution functions
    def element_res_convolution(wav, R, wav_extended, flux_conv_rot, fwhm_lim,
                                normalize: bool = True):
        """Embarrassingly parallel component of resolution convolution"""
        fwhm = wav / R
        # Mask of wavelength range within fwhm_lim* fwhm of wav
        fwhm_space = fwhm_lim * fwhm
        index_mask = mask_between(wav_extended, wav - fwhm_space, wav + fwhm_space)

        flux_2convolve = flux_conv_rot[index_mask]
        # Gaussian Instrument Profile for given resolution and wavelength
        IP = unitary_gaussian(wav_extended[index_mask], wav, fwhm)

        sum_val = np.sum(IP * flux_2convolve)
        if normalize:
            # Correct for the effect of convolution with non-equidistant positions
            unitary_val = np.sum(IP)  # Affects precision
            return sum_val / unitary_val
        else:
            return sum_val

    def wrapper_res_parallel_convolution(args):
        """Wrapper for res_parallel_convolution needed to unpack the arguments
        for fast_convolve as multiprocess.Pool.map does not accept multiple
        arguments
        """
        return element_res_convolution(*args)

    if num_procs != 0:
        if num_procs is None:
            num_procs = mprocess.cpu_count() - 1

        mproc_pool = mprocess.Pool(processes=num_procs)
        # Need to update the values here
        args_generator = tqdm([[wav, R, wav_extended, flux_conv_rot, fwhm_lim,
                                normalize] for wav in wav_band])
        flux_conv_res = np.array(mproc_pool.map(wrapper_res_parallel_convolution,
                                                args_generator))
        mproc_pool.close()

    else:  # num_procs == 0
        flux_conv_res = np.empty_like(wav_band)  # Memory assignment
        for jj, wav in enumerate(tqdm(wav_band)):
            flux_conv_res[jj] = element_res_convolution(wav, R, wav_extended,
                                                        flux_conv_rot, fwhm_lim,
                                                        normalize=normalize)
        print("Done.\n")
    return flux_conv_res


###############################################################################
def name_assignment(spectrum: str):
    """
    assigns a name to the filename in which the spectrum is going to be saved
    """
    # Simplified to temperature and base in spectrum name.
    m0_aces = "lte03900"
    m3_aces = "lte03500"
    m6_aces = "lte02800"
    m9_aces = "lte02600"
    base = "PHOENIX-ACES-AGSS-COND-2011-HiRes_wave.dat"
    if (m0_aces in spectrum) and (base in spectrum):
        name = "M0-PHOENIX-ACES"
    elif (m3_aces in spectrum) and (base in spectrum):
        name = "M3-PHOENIX-ACES"
    elif (m6_aces in spectrum) and (base in spectrum):
        name = "M6-PHOENIX-ACES"
    elif (m9_aces in spectrum) and (base in spectrum):
        name = "M9-PHOENIX-ACES"
    else:
        raise ValueError("Name {0} not found!".format(spectrum))
    return name


if __name__ == "__main__":
    if len(sys.argv) == 3:
        run_convolutions(sys.argv[1], sys.argv[2])
    else:
        print("Arguments not compatible with called function.")

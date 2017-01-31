"""
Created on Sun Dec 14 15:43:13 2014

@author: pfigueira

Adapted December 2016 by Jason Neal
"""
from __future__ import division, print_function
import sys
import numpy as np
from tqdm import tqdm
import multiprocess as mprocess

# from eniric.IOmodule import read_2col, read_3col
from eniric.IOmodule import pdread_2col,  write_e_3col
# from eniric.Qcalculator import RVprec_calc, SqrtSumWis
from eniric.utilities import wav_selector, unitary_Gauss, rotation_kernel
from eniric.utilities import band_selector, read_spectrum
import matplotlib.pyplot as plt
from matplotlib import rc
# set stuff for latex usage
rc('text', usetex=True)

data_rep = "../data/PHOENIX_ACES_spectra/"
results_dir = "../data/results/"
resampled_dir = "../data/resampled/"

# models form PHOENIX-ACES
M0_ACES = data_rep+"lte03900-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes_wave.dat"
M3_ACES = data_rep+"lte03500-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes_wave.dat"
M6_ACES = data_rep+"lte02800-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes_wave.dat"
M9_ACES = data_rep+"lte02600-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes_wave.dat"


def run_convolutions(spectrum_string, band):
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


def save_convolution_results(filename, wavelength, flux, convolved_flux):
    """Saves covolution results to a file.

    Parameters
    ----------
    filename: str
    wavelength: array-like
     flux, convolved_flux
    """
    print("Saving results...")

    # Note: difference in sampling at 1.0 and 1.5 microns makes jumps
    # in the beginning of Y and H bands
    write_e_3col(filename, wavelength, flux, convolved_flux)
    print("Done.")
    return 0


def convolve_spectra(spectrum, band, vsini, R, epsilon=0.6, fwhm_lim=5.0,
                     plot=True, num_procs=None, results_dir=results_dir,
                     data_rep=data_rep, normalize=True, output_name=None):
    """Load Spectrum, apply convolution and then save results.

    """
    print("Reading the data...")
    wav, flux = read_spectrum(spectrum)    # In microns and  photon flux.
    print("Done.")

    wav_band, flux_band, convolved_flux = convolution(wav, flux, vsini, R, band,
                                                      epsilon=epsilon, fwhm_lim=fwhm_lim,
                                                      num_procs=num_procs, normalize=normalize)

    if output_name is None:
        name_model = name_assignment(spectrum)
        if normalize:
            filename = ("{0}Spectrum_{1}_{2}band_vsini{3:3.1f}_R{4:d}k.txt"
                        "").format(results_dir, name_model, band, vsini, R/1000)
        else:
            filename = ("{0}Spectrum_{1}_{2}band_vsini{3:3.1f}_R{4:d}k_unnormalized.txt"
                        "").format(results_dir, name_model, band, vsini, R/1000)
    else:
        filename = results_dir + output_name

    save_convolution_results(filename, wav_band, flux_band, convolved_flux)

    if plot:
        fig = plt.figure(1)
        plt.xlabel(r"wavelength [$\mu$m])")
        plt.ylabel(r"flux [counts] ")
        plt.plot(wav_band, flux_band/np.max(flux_band), color='k', linestyle="-", label="Original spectra")
        plt.plot(wav_band, convolved_flux/np.max(convolved_flux), color='b', linestyle="-", label="{0:s} spectrum observed at vsini={1:.2f} and R={2:d} .".format(name_model, vsini, R))
        plt.legend(loc='best')
        plt.show()

        fig.savefig(filename[:-3]+"pdf", facecolor='w', format='pdf', bbox_inches='tight')
        plt.close()

    return 0


def convolution(wav, flux, vsini, R, band="All", epsilon=0.6, fwhm_lim=5.0,
                num_procs=None, normalize=True, output_name=None):
    """Perform convolution of spectrum.

    Rotational convolution followed by a guassian a a specified resolution R.

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
        Number of processes to use with multiprocess. If None it is asigned to 1 less then total number of cores.
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
    fwhm_min = wav_band[0] / R    # fwhm at the extremes of vector
    fwhm_max = wav_band[-1] / R

    # performing convolution with rotation kernel
    print("Starting the Rotation convolution for vsini={0:.2f}...".format(vsini))

    delta_lambda_min = wav_band[0] * vsini / 3.0e5
    delta_lambda_max = wav_band[-1] * vsini / 3.0e5

    # widest wavelength bin for the rotation convolution
    wav_ext_rotation, flux_ext_rotation = wav_selector(wav, flux, wav_band[0]-delta_lambda_min-fwhm_lim*fwhm_min, wav_band[-1]+delta_lambda_max+fwhm_lim*fwhm_max)

    # wide wavelength bin for the resolution_convolution
    wav_extended, flux_extended = wav_selector(wav, flux, wav_band[0]-fwhm_lim*fwhm_min, wav_band[-1]+fwhm_lim*fwhm_max)

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


def rotational_convolution(wav_extended, wav_ext_rotation, flux_ext_rotation,
                           vsini, epsilon, num_procs=None, normalize=True):
    """Perform Rotational convolution part of convolution.
    """

    def wrapper_rot_parallel_convolution(args):
        """Wrapper for rot_parallel_convolution needed to unpack the arguments for
        fast_convolve as multiprocess.Pool.map does not accept multiple
        arguments
        """
        return element_rot_convolution(*args)

    def element_rot_convolution(wav, wav_extended, wav_ext_rotation,
                                flux_ext_rotation, vsini, epsilon,
                                normalize):
        """Embarisingly parallel part of rotational convolution"""
        # select all values such that they are within the fwhm limits
        delta_lambda_l = wav * vsini / 3.0e5

        index_mask = ((wav_ext_rotation > (wav - delta_lambda_l)) &
                      (wav_ext_rotation < (wav + delta_lambda_l)))
        flux_2convolve = flux_ext_rotation[index_mask]
        rotation_profile = rotation_kernel(wav_ext_rotation[index_mask] - wav, delta_lambda_l, vsini, epsilon)

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
                                                        normalize)
        print("Done.\n")
    return flux_conv_rot


def resolution_convolution(wav_band, wav_extended, flux_conv_rot, R, fwhm_lim,
                           num_procs=1, normalize=True):
    """Perform Resolution convolution part of convolution.
    """

    # Define inner convolution functions
    def element_res_convolution(wav, R, wav_extended, flux_conv_rot, fwhm_lim,
                                normalize):
        """Embarisingly parallel component of resolution convolution"""
        fwhm = wav / R
        # Mask of wavelength range within 5 fwhm of wav
        index_mask = ((wav_extended > (wav - fwhm_lim*fwhm)) &
                      (wav_extended < (wav + fwhm_lim*fwhm)))

        flux_2convolve = flux_conv_rot[index_mask]
        # Gausian Instrument Profile for given resolution and wavelength
        IP = unitary_Gauss(wav_extended[index_mask], wav, fwhm)

        sum_val = np.sum(IP * flux_2convolve)
        if normalize:
            # Correct for the effect of convolution with non-equidistant postions
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
        flux_conv_res = np.empty_like(wav_band)   # Memory assignment
        for jj, wav in enumerate(tqdm(wav_band)):
            flux_conv_res[jj] = element_res_convolution(wav, R, wav_extended,
                                                        flux_conv_rot, fwhm_lim)
        print("Done.\n")
    return flux_conv_res


###############################################################################
def name_assignment(spectrum):
    """
    assigns a name to the filename in which the spectrum is going to be saved
    """
    # Simplified to temperature and base in spectrum name.
    M0_ACES = "lte03900"
    M3_ACES = "lte03500"
    M6_ACES = "lte02800"
    M9_ACES = "lte02600"
    base = "PHOENIX-ACES-AGSS-COND-2011-HiRes_wave.dat"
    if (M0_ACES in spectrum) and (base in spectrum):
        name = "M0-PHOENIX-ACES"
    elif(M3_ACES in spectrum) and (base in spectrum):
        name = "M3-PHOENIX-ACES"
    elif(M6_ACES in spectrum) and (base in spectrum):
        name = "M6-PHOENIX-ACES"
    elif(M9_ACES in spectrum) and (base in spectrum):
        name = "M9-PHOENIX-ACES"
    else:
        print("Name {0} not found!".format(spectrum))
        exit(1)
    return name


if __name__ == "__main__":
    if len(sys.argv) == 3:
        run_convolutions(sys.argv[1], sys.argv[2])
    else:
        print("Arguments not compatible with called functtion.")

"""
Created on Sun Dec 14 15:43:13 2014

@author: pfigueira

Adapted December 2016 by Jason Neal
"""
from __future__ import division, print_function
import sys
import numpy as np
from tqdm import tqdm
import pandas as pd
import multiprocess as mprocess
from os import listdir
from os.path import isfile, join

# from eniric.IOmodule import read_2col, read_3col
# from eniric.Qcalculator import RVprec_calc, SqrtSumWis

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


def read_spectrum(spec_name):
    """
    function that reads spectra from the database
    """
    data = pd.read_table(spec_name, comment='#', names=["wavelength", "flux"],
                         dtype=np.float64, delim_whitespace=True)
    wav, flux = data["wavelength"].values, data["flux"].values
    wav *= 1.0e-4  # conversion to microns

    flux_photons = flux * wav

    return wav, flux_photons


def band_selector(wav, flux, band):
    """ Select a specific wavelength band.

    Parameters
    ----------
    wav: array-like
        Wavelength values.
    flux: array-like
        Flux values.
    band: str
        Band letter to select, upper or lower case is supported. Options
        are ("ALL" or ""), "VIS", "GAP", "Z", "Y", "J", "H", "K".
    """
    band = band.upper()

    bands = {"VIS": (0.38, 0.78), "GAP": (0.78, 0.83), "Z": (0.83, 0.93),
             "Y": (1.0, 1.1), "J": (1.17, 1.33), "H": (1.5, 1.75),
             "K": (2.07, 2.35)}
    if(band in ["ALL", ""]):
        return wav, flux
    elif band in bands:
        # select values form the band
        bandmin = bands[band][0]
        bandmax = bands[band][1]

        return wav_selector(wav, flux, bandmin, bandmax)
    else:
        print("Unrecognized band tag.")
        exit(1)


def plotter(spectrum, band, vsini=0, R=0):
    """
    Reads and plots the selected spectrum in a given band
    """
    wav, flux = read_spectrum(spectrum)
    wav_band, flux_band = band_selector(wav, flux, band)

    plt.figure(1)
    plt.xlabel(r"wavelength [$\mu$m])")
    plt.ylabel(r"flux [counts] ")
    plt.plot(wav_band, flux_band, color='k', marker="o", linestyle="-")
    plt.show()
    plt.close()

def run_convolutions(spectrum_string, band):
    """
    Runs the convolutions for a set of spectra in batch
    """
    vsini = [1.0, 5.0, 10.0]
    R = [60000, 80000, 100000]

    exec('spectrum = ' + spectrum_string)        #note: warnings to be dismissed, due to exec usage
    print(spectrum)
    print("Running the convolutions for spectra of %s in band %s\n." % (spectrum, band))
    for vel in vsini:
        for res in R:
            convolution(spectrum, band, vel, res, plot=False)

def convolution(spectrum, band, vsini, R, epsilon=0.6, FWHM_lim=5.0, plot=True, numProcs=None, data_rep=data_rep, results_dir=results_dir):

    """
    function that convolves a given spectra to a resolution of R
    R = 60 000 , R = 80 000, R = 100 000

    Parameters
    ----------
    numProcs: int, None
        Number of processes to use with multiprocess. If None it is asigned to 1 less then total number of cores.
        If numProcs = 0, then multiprocess is not used.

    Returns
    -------
    Saves figures to results_dir.

    """

    print("Reading the data...")
    wav, flux = read_spectrum(spectrum)
    wav_band, flux_band = band_selector(wav, flux, band)
    print("Done.")

    # we need to calculate the FWHM at this value in order to set the starting point for the convolution
    FWHM_min = wav_band[0] / R    # FWHM at the extremes of vector
    FWHM_max = wav_band[-1] / R

    # performing convolution with rotation kernel
    print("Starting the Rotation convolution for vsini=%.2f..." % (vsini))

    delta_lambda_min = wav_band[0]*vsini/3.0e5
    delta_lambda_max = wav_band[-1]*vsini/3.0e5
    # widest wavelength bin for the rotation convolution
    wav_ext_rotation, flux_ext_rotation = wav_selector(wav, flux, wav_band[0]-delta_lambda_min-FWHM_lim*FWHM_min, wav_band[-1]+delta_lambda_max+FWHM_lim*FWHM_max)
    wav_ext_rotation = np.array(wav_ext_rotation, dtype="float64")
    flux_ext_rotation = np.array(flux_ext_rotation, dtype="float64")

    # wide wavelength bin for the resolution_convolution
    wav_extended, flux_extended = wav_selector(wav, flux, wav_band[0]-FWHM_lim*FWHM_min, wav_band[-1]+FWHM_lim*FWHM_max)
    wav_extended = np.array(wav_extended, dtype="float64")
    flux_extended = np.array(flux_extended, dtype="float64")

    # rotational convolution
    flux_conv_rot = rotational_convolution(wav_extended, wav_ext_rotation, flux_ext_rotation, vsini, epsilon)

    print("Starting the Resolution convolution...")

    flux_conv_res = resolution_convolution(wav_band, wav_extended, flux_conv_rot, R, FWHM_lim, numProcs=numProcs)

    print("Saving results...")

    # Note: difference in sampling at 1.0 and 1.5 microns makes jumps in the beginning of Y and H bands

    name_model = name_assignment(spectrum)
    filename = results_dir+"Spectrum_"+name_model+"_"+band+"band_vsini"+str(vsini)+"_R"+str(int(R/1000))+"k.txt"
    write_3col(filename, wav_band, flux_band, flux_conv_res)
    print("Done.")

    if(plot):
        fig = plt.figure(1)
        plt.xlabel(r"wavelength [$\mu$m])")
        plt.ylabel(r"flux [counts] ")
        plt.plot(wav_band, flux_band/np.max(flux_band), color='k', linestyle="-", label="Original spectra")
        plt.plot(wav_band, flux_conv_res/np.max(flux_conv_res), color='b', linestyle="-", label="%s spectrum observed at vsini=%.2f and R=%d ." % (name_model, vsini, R))
        plt.legend(loc='best')
        plt.show()

        fig.savefig(filename[:-3]+"pdf", facecolor='w', format='pdf', bbox_inches='tight')
        plt.close()

    return wav_band, flux_conv_res


def rotational_convolution(wav_extended, wav_ext_rotation, flux_ext_rotation, vsini, epsilon, numProcs=None):
    """ Perform Rotational convolution part of convolution.
    """

    def wrapper_rot_parallel_convolution(args):
        """ Wrapper for rot_parallel_convolution needed to unpack the arguments for
        fast_convolve as multiprocess.Pool.map does not accept multiple
        arguments
        """
        return element_rot_convolution(*args)

    def element_rot_convolution(wav, wav_extended, wav_ext_rotation, flux_ext_rotation, vsini, epsilon):
        """Embarisingly parallel part of rotational convolution"""
        # select all values such that they are within the FWHM limits
        delta_lambda_L = wav * vsini / 3.0e5

        index_mask = ((wav_ext_rotation > (wav - delta_lambda_L)) &
                      (wav_ext_rotation < (wav + delta_lambda_L)))
        flux_2convolve = flux_ext_rotation[index_mask]
        rotation_profile = rotation_kernel(wav_ext_rotation[index_mask] - wav, delta_lambda_L, vsini, epsilon)
        return np.sum(rotation_profile * flux_2convolve)

    if numProcs != 0:
        if numProcs is None:
            numProcs = mprocess.cpu_count() - 1

        mprocPool = mprocess.Pool(processes=numProcs)

        args_generator = tqdm([[wav, wav_extended, wav_ext_rotation,
                                flux_ext_rotation, vsini, epsilon]
                              for wav in wav_extended])

        flux_conv_rot = np.array(mprocPool.map(wrapper_rot_parallel_convolution,
                                 args_generator))

        mprocPool.close()

    else:  # numProcs == 0
        flux_conv_rot = np.empty_like(wav_extended)  # Memory assignment
        for ii, wav in enumerate(tqdm(wav_extended)):
            flux_conv_rot[ii] = element_rot_convolution(wav, wav_extended,
                                                        wav_ext_rotation,
                                                        flux_ext_rotation,
                                                        vsini, epsilon)
        print("Done.\n")
    return flux_conv_rot



def resolution_convolution(wav_band, wav_extended, flux_conv_rot, R, FWHM_lim, numProcs=1):
    """ Perform Resolution convolution part of convolution.
    """

    # Define inner convolution functions
    def element_res_convolution(wav, R, wav_extended, flux_conv_rot, FWHM_lim):
        """ Embarisingly parallel component of resolution convolution"""
        FWHM = wav / R
        # Mask of wavelength range within 5 FWHM of wav
        index_mask = ((wav_extended > (wav - FWHM_lim*FWHM)) &
              (wav_extended < (wav + FWHM_lim*FWHM)))

        flux_2convolve = flux_conv_rot[index_mask]
        # Gausian Instrument Profile for given resolution and wavelength
        IP = unitary_Gauss(wav_extended[index_mask], wav, FWHM)

        sum_val = np.sum(IP * flux_2convolve)
        # Correct for the effect of convolution with non-equidistant postions
        # unitary_val = np.sum(IP * np.ones_like(flux_2convolve))  # Affects precision
        return sum_val # / unitary_val


    def wrapper_res_parallel_convolution(args):
        """ Wrapper for res_parallel_convolution needed to unpack the arguments for
        fast_convolve as multiprocess.Pool.map does not accept multiple
        arguments
        """
        return element_res_convolution(*args)

    if numProcs != 0:
        if numProcs is None:
            numProcs = mprocess.cpu_count() - 1

        mprocPool = mprocess.Pool(processes=numProcs)
        # Need to update the values here
        args_generator = tqdm([[wav, R, wav_extended, flux_conv_rot, FWHM_lim]
                              for wav in wav_band])
        flux_conv_res = np.array(mprocPool.map(wrapper_res_parallel_convolution,
                                               args_generator))
        mprocPool.close()

    else:  # numProcs == 0
        flux_conv_res = np.empty_like(wav_band)   # Memory assignment
        for jj, wav in enumerate(tqdm(wav_band)):
            flux_conv_res[jj] = element_res_convolution(wav, R, wav_extended,
                                                        flux_conv_rot, FWHM_lim)
        print("Done.\n")
    return flux_conv_res


###############################################################################
#
#   Auxiliary functions
#
###############################################################################

def wav_selector(wav, flux, wav_min, wav_max):
    """
    function that returns wavelength and flux withn a giving range

    Parameters
    ----------
    wav: array-like
        Wavelength array.
    flux: array-like
        Flux array.
    wav_min: float
        Lower bound wavelength value.
    wav_max: float
        Upper bound wavelength value.

    Returns
    -------
    wav_sel: array
        New wavelength array within bounds wav_min, wav_max
    flux_sel: array
        New wavelength array within bounds wav_min, wav_max
        """
    wav = np.asarray(wav, dtype="float64")
    flux = np.asarray(flux, dtype="float64")

    mask = (wav > wav_min) & (wav < wav_max)
    flux_sel = flux[mask]
    wav_sel = wav[mask]

    return wav_sel, flux_sel


def unitary_Gauss(x, center, FWHM):
    """ Gaussian function of area = 1.

    Parameters
    ----------
    x: array-like
        Position array
    center: float
        Central postion of Gaussian
    FHWM: float
        Full Width at Half Maximum

    Returns
    -------
    result: array-like
        Result of gaussian function sampled at x values.
    """

    sigma = np.abs(FWHM) / (2 * np.sqrt(2 * np.log(2)))
    Amp = 1.0 / (sigma * np.sqrt(2 * np.pi))
    tau = -((x - center)**2) / (2 * (sigma**2))
    result = Amp * np.exp(tau)

    return result


def rotation_kernel(delta_lambdas, delta_lambda_L, vsini, epsilon):
    """ Calculate the rotation kernel for a given wavelength

    Parameters
    ----------
    delta_lambdas: array
        Wavelength values selected within delta_lambda_L around central value. (check)
    delta_lambda_L: float
        FWHM of rotational broading. (check)
    vsini: float
        Projected rotational velocity [km/s]
    epsilon: float
        Linear limb-darkening coefficient (0-1).

    Returns
    -------
        Rotational kernal

    Notes:
    Equations * from .... book.

    """
    denominator = (np.pi * vsini * (1.0 - epsilon / 3.0))
    lambda_ratio_sqr = (delta_lambdas / delta_lambda_L)**2.0

    c1 = 2.0 * (1.0 - epsilon) / denominator
    c2 = 0.5 * np.pi * epsilon / denominator

    return (c1 * np.sqrt(1.0-lambda_ratio_sqr) + c2 * (1.0-lambda_ratio_sqr))


###############################################################################
def resample_allfiles(results_dir=results_dir, resampled_dir=resampled_dir):
    """
    reample all files inside folder
    Parameters
    ----------
    results_dir: str
        Directory containing results to resample.
    """
    # getting a list of all the files
    onlyfiles = [f for f in listdir(results_dir) if isfile(join(results_dir, f))]

    [resampler(spectrum_file, results_dir=results_dir,
               resampled_dir=resampled_dir) for spectrum_file in onlyfiles
     if spectrum_file[-4:] == ".txt"]


def resampler(spectrum_name="Spectrum_M0-PHOENIX-ACES_Yband_vsini1.0_R60k.txt",
              results_dir=results_dir, resampled_dir=resampled_dir,
              sampling=3.0, plottest=False):
    """
    resamples a spectrum by interpolation onto a grid with a sampling of 3 pixels per resolution element
    """
    # wavelength, theoretical_spectrum, spectrum = read_3col(spectrum_name)
    read_name = results_dir + spectrum_name
    data = pd.read_table(read_name, header=None, dtype=np.float64,
                         names=["wavelength", "model", "spectrum"],
                         delim_whitespace=True)
    wavelength = data["wavelength"].values
    # theoretical_spectrum = data["model"].values
    spectrum = data["spectrum"].values

    wavelength_start = wavelength[1]  # because of border effects
    wavelength_end = wavelength[-2]   # because of border effects
    resolution_string = spectrum_name[-8:-5]

    if(resolution_string[0] == "R"):
        resolution = int(resolution_string[1:])*1000
    else:
        resolution = int(resolution_string)*1000

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
    write_2col(filetowrite, wav_grid, interpolated_flux)

    if(plottest):
        plt.figure(1)
        plt.xlabel(r"wavelength [$\mu$m])")
        plt.ylabel(r"flux [counts] ")
        plt.plot(wavelength, spectrum, color='k', linestyle="-", label="Original spectrum")
        plt.plot(wav_grid, interpolated_flux, color='b', linestyle="-", label="Interpolated spectrum")
        plt.legend(loc='best')
        plt.show()

        plt.close()


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
        print("Name {} not found!".format(spectrum))
        exit(1)
    return name


def write_2col(filename, data1, data2):
    # Writes data in 2 columns separated by tabs in a "filename" file.

    f = open(filename, "w")

    for i in range(len(data1)):
        # f.write("\t"+str(data1[i])+"\t\t"+str(data2[i])+"\t\t"+str(data3[i])+"\n")
        f.write("\t%e\t\t%e\n" % (data1[i], data2[i]))

    f.close()


def write_3col(filename, data1, data2, data3):
    # Writes data in 3 columns separated by tabs in a "filename" file.

    f = open(filename, "w")

    for i in range(len(data1)):
        # f.write("\t"+str(data1[i])+"\t\t"+str(data2[i])+"\t\t"+str(data3[i])+"\n")
        f.write("\t%e\t\t%e\t\t%e\n" % (data1[i], data2[i], data3[i]))

    f.close()
###############################################################################


def list_creator(spectrum, band):
    """
    creates a list of potential lines from a brute-force analysis of the band
    """
    wav, flux = read_spectrum(spectrum)
    wav_band, flux_band = band_selector(wav, flux, band)

    print(band + " band list:")
    short_flux = flux_band[2:-2]
    left_mask = ((short_flux < flux_band[:-4]) &
                 (short_flux < flux_band[1:-3]) &
                 (flux_band[:-4] > flux_band[1:-3]))

    right_mask = mask_right = ((short_flux < flux_band[3:-1]) &
                               (short_flux < flux_band[4:]) &
                               (flux_band[4:] > flux_band[3:-1]))

    line_centers = wav_band[2:-2][left_mask * right_mask]  # find peaks using masking
    print("Line centers", line_centers * 1.0e4)
    print("In a spectrum with {} points".format(len(wav_band)),
          ", {} lines were found.".format(len(line_centers)))
    return line_centers


###############################################################################

if __name__ == "__main__":
    if len(sys.argv) == 3 :
        run_convolutions(sys.argv[1], sys.argv[2])
    else:
        print("Arguments not compatible with called functtion.")

"""
Created on Fri Feb  6 15:42:03 2015

@author: pfigueira

Updated for eniric/python3 - Janurary 2017
@author: Jason Neal
"""
import sys
import argparse
import itertools
import numpy as np
from sys import exit
import matplotlib.pyplot as plt

# to remove labels in one tick
from matplotlib.ticker import MaxNLocator

import eniric.IOmodule as IOmodule
import eniric.Qcalculator as Qcalculator
import eniric.utilities as utils
from eniric.utilities import band_selector

import eniric.plotting_functions as plt_functions

from matplotlib import rc
# set stuff for latex usage
rc('text', usetex=True)

def _parser():
    """Take care of all the argparse stuff.

    :returns: the args
    """
    parser = argparse.ArgumentParser(description='Calculate radial velocity precision of model spectra.')

    parser.add_argument("-b", "--bands", type=str, default="J",
                        choices=["ALL", "VIS", "GAP", "Z", "Y", "J", "H", "K", None],
                        help="Wavelength bands to select. Default=J.", nargs="+")
    parser.add_argument("--plot_bary", default=False, action="store_true",
                        help="Plot the barycentric shift and exit.")
    args = parser.parse_args()
    return args


def main(bands="J", plot_bary=False):
    """ Main function that calls calc_precision


    Parameters
    ----------

    bands: str or list of str or None, Default="J"
        Band letters to use. None does the bands Z through K.
    plot_bary: bool
        Flag to plot and test the barycentric masking then exit.
    """

    resampled_dir = "resampled_cont/"
    resampled_dir_OLD = "resampled/"

    spectral_types = ["M0", "M3", "M6", "M9"]
    if isinstance(bands, str):
        bands = [bands]
    elif (bands is None) or (bands is "None"):
        bands = ["Z", "Y", "J", "H", "K"]

    vsini = ["1.0", "5.0", "10.0"]
    R = ["60k", "80k", "100k"]
    sampling = ["3"]

    results = calculate_prec(bands, plot_bary=plot_bary, plot_atm=False, plot_ste=False, plot_flux=False, paper_plots=False, offset_RV=0.0)

def read_contfit():
    """
    function that reads continuum fitting
    """
    M0_contfit = "PHOENIX_ACES_spectra/M0_nIRcont.txt"
    M3_contfit = "PHOENIX_ACES_spectra/M3_nIRcont.txt"
    M6_contfit = "PHOENIX_ACES_spectra/M6_nIRcont.txt"
    M9_contfit = "PHOENIX_ACES_spectra/M9_nIRcont.txt"

    wav_M0, flux_M0 = IOmodule.pdread_2col(M0_contfit)
    wav_M3, flux_M3 = IOmodule.pdread_2col(M3_contfit)
    wav_M6, flux_M6 = IOmodule.pdread_2col(M6_contfit)
    wav_M9, flux_M9 = IOmodule.pdread_2col(M9_contfit)

    wav_M0 = np.array(wav_M0, dtype="float64")*1.0e-4  # conversion to microns
    flux_M0 = np.array(flux_M0, dtype="float64")*wav_M0
    wav_M3 = np.array(wav_M3, dtype="float64")*1.0e-4  # conversion to microns
    flux_M3 = np.array(flux_M3, dtype="float64")*wav_M3
    wav_M6 = np.array(wav_M6, dtype="float64")*1.0e-4  # conversion to microns
    flux_M6 = np.array(flux_M6, dtype="float64")*wav_M6
    wav_M9 = np.array(wav_M9, dtype="float64")*1.0e-4  # conversion to microns
    flux_M9 = np.array(flux_M9, dtype="float64")*wav_M9

    return [wav_M0, flux_M0, wav_M3, flux_M3, wav_M6, flux_M6, wav_M9, flux_M9]


def read_nIRspectra():
    """
    function that reads nIR spectra
    """
    M0_contfit = "PHOENIX_ACES_spectra/lte03900-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes_wave_CUT_nIR.dat"
    M3_contfit = "PHOENIX_ACES_spectra/lte03500-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes_wave_CUT_nIR.dat"
    M6_contfit = "PHOENIX_ACES_spectra/lte02800-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes_wave_CUT_nIR.dat"
    M9_contfit = "PHOENIX_ACES_spectra/lte02600-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes_wave_CUT_nIR.dat"

    print("Reading PHOENIX original spectrum...")
    wav_M0, flux_M0 = IOmodule.pdread_2col(M0_contfit)
    wav_M3, flux_M3 = IOmodule.pdread_2col(M3_contfit)
    wav_M6, flux_M6 = IOmodule.pdread_2col(M6_contfit)
    wav_M9, flux_M9 = IOmodule.pdread_2col(M9_contfit)
    print("Done.")

    wav_M0 = np.array(wav_M0, dtype="float64")*1.0e-4  # conversion to microns
    flux_M0 = np.array(flux_M0, dtype="float64")*wav_M0
    wav_M3 = np.array(wav_M3, dtype="float64")*1.0e-4  # conversion to microns
    flux_M3 = np.array(flux_M3, dtype="float64")*wav_M3
    wav_M6 = np.array(wav_M6, dtype="float64")*1.0e-4  # conversion to microns
    flux_M6 = np.array(flux_M6, dtype="float64")*wav_M6
    wav_M9 = np.array(wav_M9, dtype="float64")*1.0e-4  # conversion to microns
    flux_M9 = np.array(flux_M9, dtype="float64")*wav_M9

    wav_M0 = wav_M0[1000:-1000]
    wav_M3 = wav_M3[1000:-1000]
    wav_M6 = wav_M6[1000:-1000]
    wav_M9 = wav_M9[1000:-1000]

    flux_M0 = moving_average(flux_M0, 2000)[1000:-1000]
    flux_M3 = moving_average(flux_M3, 2000)[1000:-1000]
    flux_M6 = moving_average(flux_M6, 2000)[1000:-1000]
    flux_M9 = moving_average(flux_M9, 2000)[1000:-1000]

    return [wav_M0, flux_M0, wav_M3, flux_M3, wav_M6, flux_M6, wav_M9, flux_M9]


def ratios_calc(wav_bin):
    """
    funtion that calculates ratios as set in the continuum to apply to the spectrum
    """
    wav_M0, flux_M0, wav_M3, flux_M3, wav_M6, flux_M6, wav_M9, flux_M9 = read_contfit()

    index_M0_l = np.searchsorted(wav_M0, [0.83, 1.0, 1.17, 1.5, 2.07])
    index_M3_l = np.searchsorted(wav_M3, [0.83, 1.0, 1.17, 1.5, 2.07])
    index_M6_l = np.searchsorted(wav_M6, [0.83, 1.0, 1.17, 1.5, 2.07])
    index_M9_l = np.searchsorted(wav_M9, [0.83, 1.0, 1.17, 1.5, 2.07])

    index_M0_r = np.searchsorted(wav_M0, [0.83+wav_bin, 1.0+wav_bin, 1.17+wav_bin, 1.5+wav_bin, 2.07+wav_bin])
    index_M3_r = np.searchsorted(wav_M3, [0.83+wav_bin, 1.0+wav_bin, 1.17+wav_bin, 1.5+wav_bin, 2.07+wav_bin])
    index_M6_r = np.searchsorted(wav_M6, [0.83+wav_bin, 1.0+wav_bin, 1.17+wav_bin, 1.5+wav_bin, 2.07+wav_bin])
    index_M9_r = np.searchsorted(wav_M9, [0.83+wav_bin, 1.0+wav_bin, 1.17+wav_bin, 1.5+wav_bin, 2.07+wav_bin])

    return [flux_bin(flux_M0, index_M0_l, index_M0_r), flux_bin(flux_M3, index_M3_l, index_M3_r), flux_bin(flux_M6, index_M6_l, index_M6_r), flux_bin(flux_M9, index_M9_l, index_M9_r)]


def flux_bin(flux, index_left, index_right):
    fluxes = []
    for ind_l, ind_r in zip(index_left, index_right):
        fluxes.append(np.average(flux[ind_l: ind_r]))
    return fluxes


def prepare_atmopshere(atmmodel):
    """ Read in atmopheric model and prepare. """
    wav_atm, flux_atm, std_flux_atm, mask_atm = IOmodule.pdread_4col(atmmodel)
    # pandas lready returns numpy arrays
    wav_atm = wav_atm / 1000.0  # conversion from nanometers to micrometers
    mask_atm = np.array(mask_atm, dtype=bool)
    return wav_atm, flux_atm, std_flux_atm, mask_atm


def old_barycenter_shift(wav_atm, mask_atm, offset_RV=0.0):
    """ Old version Calculating impact of Barycentric movement on mask...

    Extends the masked region to +-30 km/s due to the barycentic motion of the earth.
    """
    pixels_total = len(mask_atm)
    masked_start = pixels_total - np.sum(mask_atm)

    mask_atm_30kms = []
    for value in zip(wav_atm, mask_atm):
        if (value[1] == False) and (offset_RV == 666.0):    # if the mask is false and the offset is equal to zero
            mask_atm_30kms.append(value[1])

        else:

            delta_lambda = value[0] * 3.0e4/Qcalculator.c.value
            starting_lambda = value[0] * offset_RV*1.0e3/Qcalculator.c.value
            indexes_30kmslice = np.searchsorted(wav_atm, [starting_lambda+value[0]-delta_lambda,
                                                          starting_lambda+value[0]+delta_lambda])
            indexes_30kmslice = [index if(index < len(wav_atm)) else len(wav_atm)-1 for index in indexes_30kmslice]

            mask_atm_30kmslice = np.array(mask_atm[indexes_30kmslice[0]:indexes_30kmslice[1]], dtype=bool)    # selecting only the slice in question

            # if(False in mask_atm_30kmslice):
            #    mask_atm_30kms.append(False)
            # else:
            #    mask_atm_30kms.append(True)

            mask_atm_30kmslice_reversed = [not i for i in mask_atm_30kmslice]

            clump = np.array_split(mask_atm_30kmslice, np.where(np.diff(mask_atm_30kmslice_reversed))[0]+1)[::2]

            tester = True
            for block in clump:
                if len(clump) >= 3:
                    tester = False
                    break

            mask_atm_30kms.append(tester)

    mask_atm = np.array(mask_atm_30kms, dtype=bool)
    masked_end = pixels_total - np.sum(mask_atm)
    print(("Old Barycentric impact affects number of masked pixels by {0:04.1%} due to the atmospheric"
          " spectrum").format((masked_end-masked_start)/pixels_total))
    print(("Pedros Pixels start = {1}, Pixel_end = {0}, Total = {2}").format(masked_end, masked_start, pixels_total))
    return mask_atm


def barycenter_shift(wav_atm, mask_atm, offset_RV=0.0):
    """ Calculating impact of Barycentric movement on mask...

    Extends the masked region to +-30 km/s due to the barycentic motion of the earth.
    """
    # Mask values to the left and right side of mask_atm. To avoid indexing errors have padded with first and last values.
    mask_iminus1 = np.concatenate(mask_atm[0], mask_atm[0], mask_atm[:-2])  # padding with first value
    mask_iplus1 = np.concatenate(mask_atm[2:], mask_atm[-1], mask_atm[-1])  # padding with last value

    pixels_total = len(mask_atm)
    masked_start = pixels_total - np.sum(mask_atm)

    barycenter_rv =  30000           # 30 km/s in m/s
    offset_rv = offset_RV * 1.0e3    # Convert to m/s

    # Doppler shift  applied to the vectors
    delta_lambdas = wav_atm * barycenter_rv / Qcalculator.c.value
    offset_lambdas = wav_atm * offset_rv / Qcalculator.c.value   # offset lambda

    # Dopler shift limits of each pixel
    wav_lower_barys = wav_atm + offset_lambdas - delta_lambdas
    wav_upper_barys = wav_atm + offset_lambdas + delta_lambdas

    mask_atm_30kms = np.empty_like(mask_atm, dtype=bool)

    for i, (wav_value, mask_val) in enumerate(zip(wav_atm, mask_atm)):
        """ If there are 3 consecutive zeros within +/-30km/s then make the value 0."""

        # Offset_RV is the offset applied for the star RV.
        if (mask_val is False) and (mask_iminus1[i] == False) and (mask_iplus1[i] == False) and (offset_RV == 0):    # if the mask is false and the offset is equal to zero
        """ If the value and its friends are already zero don't do the barycenter shifts"""
            mask_atm_30kms[i] = False
        else:
            # np.searchsorted is faster then the boolean masking wavlength range
            slice_limits = np.searchsorted(wav_atm, [wav_lower_barys[i], wav_upper_barys[i]]) # returns index to place the two shifted values
            slice_limits = [index if(index < len(wav_atm)) else len(wav_atm)-1 for index in slice_limits]  # replace index if above lenght of array
            mask_atm_slice = mask_atm[slice_limits[0]:slice_limits[1]]    # selecting only the slice in question

            mask_atm_slice = np.asarray(mask_atm_slice, dtype=bool)    # Assuring type bool

            # Make mask value false if there are 3 or more consecutive zeros in slice.
            len_consec_zeros = consecutive_truths(mask_atm_slice == False)
            if np.max(len_consec_zeros) >= 3:  # Invert mask to make zeros true with ~
                mask_atm_30kms[i] = False
            else:
                mask_atm_30kms[i] = True
                if np.sum(~mask_atm_slice) > 3:
                    print("There were {0} zeros in this barycentric shift but none were 3 consecutive!".format(np.sum(~mask_atm_slice)))

    masked_end = pixels_total - np.sum(mask_atm_30kms)
    print(("New Barycentric impact affects the number of masked pixels by {0:04.1%} due to the atmospheric"
          " spectrum").format((masked_end-masked_start)/pixels_total))
    print(("Masked Pixels start = {1}, masked_pixel_end = {0}, Total = {2}").format(masked_end, masked_start, pixels_total))
    return mask_atm_30kms


def consecutive_truths(condition):
    """ Length of consecutive true values in an bool array.

    Parameters
    ----------
    condition: ndarray of bool
        True False array of a condition.

    Returns
    -------
    len_consecutive: ndarray of ints
        Array of lengths of consecutive true values of the condition.

    Notes
    -----
    Solution found at http://stackoverflow.com/questions/24342047/count-consecutive-occurences-of-values-varying-in-length-in-a-numpy-array
    """
    if not np.any(condition):  # No match to condition
        return np.array([0])
    else:
        unequal_consec = np.concatenate(([condition[0]], condition[:-1] != condition[1:], [True]))
        where_changes = np.where(unequal_consec)[0]         # indices where condition changes
        len_consecutive = np.diff(where_changes)[::2]       # step through every second to get the "True" lenghts.
    return len_consecutive

def normalize_flux(flux_stellar, id_string):
    """Normalize flux to have SNR of 100 in middle of J band."""

    if("M0" in id_string):
        norm_constant = 1607

    elif("M3" in id_string):
        norm_constant = 1373

    elif("M6" in id_string):
        if("1.0" in id_string):
            norm_constant = 933
        elif("5.0" in id_string):
            norm_constant = 967
        else:
            norm_constant = 989

    elif("M9" in id_string):
        if("1.0" in id_string):
            norm_constant = 810
        elif("5.0" in id_string):
            norm_constant = 853
        else:
            norm_constant = 879
    else:
        print("Constant not defined. Aborting...")
        exit()

    return flux_stellar / ((norm_constant / 100.0)**2.0)


def calculate_prec(bands, plot_bary=False, plot_atm=False, plot_ste=False, plot_flux=True, paper_plots=True, offset_RV=0.0):

    for band in bands:

        atmmodel = "../data/atmmodel/Average_TAPAS_2014_{}.txt".format(band)
        print("Reading atmospheric model...")
        wav_atm, flux_atm, std_flux_atm, mask_atm = prepare_atmopshere(atmmodel)
        print(("There were {0:d} unmasked pixels out of {1:d}., or {2:.1%}."
              "").format(np.sum(mask_atm), len(mask_atm), np.sum(mask_atm) / len(mask_atm)))
        print("The model ranges from {0:4.2f} to {1:4.2f} micron.".format(wav_atm[0], wav_atm[-1]))
        print("Done.")

        print("Calculating impact of Barycentric movement on mask...")

        if plot_bary: # Ploting the two masks alongside the flux
            # Shorten arrays to make quicker
            save_results = True
            if not save_results:
                __ , flux_atm = utils.wav_selector(wav_atm, flux_atm, 2.135, 2.137)
                wav_atm, mask_atm = utils.wav_selector(wav_atm, mask_atm, 2.135, 2.137)

            new_mask_atm = barycenter_shift(wav_atm, mask_atm, offset_RV=offset_RV)
            old_mask_atm = old_barycenter_shift(wav_atm, mask_atm, offset_RV=offset_RV)  # Extend masked regions

            plt.plot(wav_atm, new_mask_atm + 0.01, "b.-", label="New Bary mask")
            plt.plot(wav_atm, old_mask_atm + 0.02, "ko-", label="Pedro Bary mask")
            plt.plot(wav_atm, mask_atm, "gs-", label="Orignal mask")
            neg30kms = wav_atm * (1 - 3e4/Qcalculator.c.value)  # doppler shift
            pos30kms = wav_atm * (1 - 3e4/Qcalculator.c.value)  # doppler shift
            plt.plot(neg30kms, mask_atm-0.02, "y", label="-30km/s")
            plt.plot(pos30kms, mask_atm-0.01, "m", label="+30km/s")
            plt.plot(wav_atm, flux_atm/np.max(flux_atm),"r--", label="Flux atm")
            plt.ylim([0.9,1.05])
            plt.legend()
            plt.show()

            if save_results:
                IOmodule.pdwrite_cols("../data/Barycenter_masking_tests.txt", wav_atm, flux_atm, std_flux_atm, mask_atm, old_mask_atm, new_mask_atm, neg30kms, pos30kms, header=["#wav_atm", "flux_atm", "std_flux_atm", "mask_atm", "Pedro mask", "Jason mask", "-30kms_wav_atm", "+30kms_wav_atm"])

            sys.exit(0)
        else:
            nmask_atm = barycenter_shift(wav_atm, mask_atm, offset_RV=offset_RV)
        print(("There were {0:d} unmasked pixels out of {1:d}, or {2:.1%}."
               "").format(np.sum(mask_atm), len(mask_atm), np.sum(mask_atm) /
                          len(mask_atm)))


        # calculating the number of pixels inside the mask
        wav_Z, mask_Z = band_selector(wav_atm, mask_atm, "Z")
        wav_Y, mask_Y = band_selector(wav_atm, mask_atm, "Y")
        wav_J, mask_J = band_selector(wav_atm, mask_atm, "J")
        wav_H, mask_H = band_selector(wav_atm, mask_atm, "H")
        wav_K, mask_K = band_selector(wav_atm, mask_atm, "K")

        bands_masked = np.concatenate((mask_Z, mask_Y, mask_J, mask_H, mask_K))

        print(("Inside the bands, there were {0:.0f} unmasked pixels out of {1:d}"
               ", or {2:.1%}.").format(np.sum(bands_masked), len(bands_masked),
                np.sum(bands_masked) / len(bands_masked)))

        if plot_atm:
            # moved ploting code to separate code, eniric.plotting_functions.py
            plt_functions.plot_atmopshere_model(wav_atm, flux_atm, mask_atm)

        # theoretical ratios calculation
        wav_M0, flux_M0, wav_M3, flux_M3, wav_M6, flux_M6, wav_M9, flux_M9 = read_nIRspectra()

        results = {}    # creating empty dictionary for the results
        wav_plot_M0 = []   # creating empty lists for the plots
        flux_plot_M0 = []
        wav_plot_M3 = []
        flux_plot_M3 = []
        wav_plot_M6 = []
        flux_plot_M6 = []
        wav_plot_M9 = []
        flux_plot_M9 = []
        #iterations = itertools.product(spectral_types, vsini, R, sampling)
        for star in spectral_types:
            for vel in vsini:
                for resolution in R:
                    for smpl in sampling:
                        file_to_read = ("Spectrum_{0}-PHOENIX-ACES_{1}band_vsini"
                                        "{2}_R{3}_res{4}.txt").format(star, band,
                                                                   vel,
                                                                   resolution,
                                                                   smpl)
                        # print("Working on "+file_to_read+".")
                        wav_stellar, flux_stellar = IOmodule.pdread_2col(resampled_dir + file_to_read)
                        # removing boundary effects
                        wav_stellar = wav_stellar[2:-2]
                        flux_stellar = flux_stellar[2:-2]

                        id_string = "{0}-{1}-{2}-{3}".format(star, band, vel,
                                                         resolution)   # sample was left aside because only one value existed

                        # Getting the wav, flux and mask values from the atm model
                        # that are the closest to the stellar wav values, see
                        # https://stackoverflow.com/questions/2566412/find-nearest-value-in-numpy-array
                        index_atm = np.searchsorted(wav_atm, wav_stellar)
                        # replace indexes outside the array, at the very end, by the value at the very end
                        # index_atm = [index if(index < len(wav_atm)) else len(wav_atm)-1 for index in index_atm]
                        indx_mask = (index_atm >= len(wav_atm))  # find broken indexs
                        index_atm[indx_mask] = len(wav_atm) - 1  # replace with index of end.

                        wav_atm_selected = wav_atm[index_atm]
                        flux_atm_selected = flux_atm[index_atm]
                        mask_atm_selected = mask_atm[index_atm]

                        # Normaize to SNR 100 in middle of J band 1.25 micron!
                        flux_stellar = normalize_flux(flux_stellar, id_string)

                        if(id_string in ["M0-J-1.0-100k", "M3-J-1.0-100k", "M6-J-1.0-100k", "M9-J-1.0-100k"]):
                            index_reference = np.searchsorted(wav_stellar, 1.25)    # searching for the index closer to 1.25 micron
                            SN_estimate = np.sqrt(np.sum(flux_stellar[index_reference-1:index_reference+2]))
                            print("\tSanity Check: The S/N for the {0:s} reference model was of {1:4.2f}.".format(id_string, SN_estimate))
                        elif("J" in id_string):
                            index_reference = np.searchsorted(wav_stellar, 1.25)    # searching for the index closer to 1.25 micron
                            SN_estimate = np.sqrt(np.sum(flux_stellar[index_reference-1:index_reference+2]))
                            print("\tSanity Check: The S/N for the {0:s} non-reference model was of {1:4.2f}.".format(id_string, SN_estimate))

                        # Precision given by the first method:
                        print("Performing analysis for: ", id_string)
                        prec_1 = Qcalculator.RVprec_calc(wav_stellar, flux_stellar)

                        # precision as given by the second_method
                        """
                        Example Joao
                        a = np.array([1, 5, 6, 8, 16, 34, 5, 7, 10, 83, 12, 6, 17, 18])
                        b = np.array([1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1], dtype=bool)

                        # this will give you a list of numpy arrays
                        c = np.array_split(a, np.where(np.diff(b))[0]+1)[::2]

                        # this will give you a list of lists
                        d = [list(cc) for cc in c]
                        print(d)
                        >>> [[1, 5], [16, 34, 5], [83, 12], [17, 18]]
                        """

                        wav_stellar_chunks_unformated = np.array_split(wav_stellar, np.where(np.diff(mask_atm_selected))[0]+1)[::2]
                        wav_stellar_chunks = [list(chunk) for chunk in wav_stellar_chunks_unformated]

                        """
                        # test section
                        print("check that lengths are the same", len(wav_stellar), len(mask_atm_selected))
                        print("size of spectra {0:d} vs number of chunks {1:d}".format(len(wav_stellar), len(wav_stellar_chunks)))
                        print("number of true elements in all chunks: {0:d}".format(len(mask_atm_selected[mask_atm_selected])))
                        """

                        flux_stellar_chunks_unformated = np.array_split(flux_stellar, np.where(np.diff(mask_atm_selected))[0]+1)[::2]
                        flux_stellar_chunks = [list(chunk) for chunk in flux_stellar_chunks_unformated]

                        """
                        # histogram checking
                        lengths = [len(chunk) for chunk in flux_stellar_chunks_unformated]
                        n, bins, patches = plt.hist(lengths, 500, range=[0.5, 500.5], histtype='stepfilled')
                        plt.title(id_string)
                        plt.show()
                        """

                        prec_2 = Qcalculator.RVprec_calc_chunks(wav_stellar_chunks, flux_stellar_chunks)

                        # Precision as given by the third_method
                        prec_3 = Qcalculator.RV_prec_calc_Trans(wav_stellar, flux_stellar, flux_atm_selected)

                        # Adding Precision results to the dictionary
                        results[id_string] = [prec_1, prec_2, prec_3]

                        # Prepare/Do for the ploting.
                        if(plot_ste or plot_ste == id_string):
                            plt_functions.plot_stellar_spectum(wav_stellar, flux_stellar, wav_atm_selected, mask_atm_selected)

                        if(plot_flux and id_string in ["M0-Z-1.0-100k", "M0-Y-1.0-100k", "M0-J-1.0-100k", "M0-H-1.0-100k", "M0-K-1.0-100k"]):
                            wav_plot_M0.append(wav_stellar)
                            flux_plot_M0.append(flux_stellar)
                        if(plot_flux and id_string in ["M3-Z-1.0-100k", "M3-Y-1.0-100k", "M3-J-1.0-100k", "M3-H-1.0-100k", "M3-K-1.0-100k"]):
                            wav_plot_M3.append(wav_stellar)
                            flux_plot_M3.append(flux_stellar)
                        if(plot_flux and id_string in ["M6-Z-1.0-100k", "M6-Y-1.0-100k", "M6-J-1.0-100k", "M6-H-1.0-100k", "M6-K-1.0-100k"]):
                            wav_plot_M6.append(wav_stellar)
                            flux_plot_M6.append(flux_stellar)
                        if(plot_flux and id_string in ["M9-Z-1.0-100k", "M9-Y-1.0-100k", "M9-J-1.0-100k", "M9-H-1.0-100k", "M9-K-1.0-100k"]):
                            wav_plot_M9.append(wav_stellar)
                            flux_plot_M9.append(flux_stellar)

    if(plot_flux):
       plt_functions.plot_nIR_flux()

    if paper_plot:
        plt_functions.plot_paper_plots()


    else:
        return results


###############################################################################
def compare_runs():
    """
    Function that compares spectra as resampled in the two versions of the code
    """
    for star in spectral_types:
        for band in bands:
            for vel in vsini:
                for resolution in R:
                    for smpl in sampling:
                        file_to_read = "Spectrum_"+star+"-PHOENIX-ACES_"+band+"band_vsini"+vel+"_R"+resolution+"_res"+smpl+".txt"
                        # print "Working on "+file_to_read+"."
                        wav_stellar, flux_stellar = IOmodule.pdread_2col(resampled_dir+file_to_read)
                        wav_stellar = wav_stellar
                        flux_stellar = flux_stellar / ((1.634e4)**2.0)

                        wav_stellar_OLD, flux_stellar_OLD = IOmodule.pdread_2col(resampled_dir_OLD+file_to_read)
                        flux_stellar_OLD = np.array(flux_stellar_OLD) / ((1.634e4)**2.0)

                        plt.figure(1)
                        plt.xlabel(r"wavelength [$\mu$m])")
                        plt.ylabel(r"Flux_stellar difference [ ] ")
                        plt.plot(wav_stellar, flux_stellar-flux_stellar_OLD, color='k')
                        plt.show()
                        plt.close()


def compare_output():
    """
    function that compares a spectrum prior to convolution, after, and after resampling
    """

    pre_convolution = "PHOENIX_ACES_spectra/lte03900-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes_wave_CUT_nIR.dat"
    pre_wav, pre_flux = IOmodule.pdread_2col(pre_convolution)
    pre_wav = np.array(pre_wav, dtype="float64")*1.0e-4  # conversion to microns
    pre_flux = np.array(pre_flux, dtype="float64")*pre_wav

    convolved = "results_new/Spectrum_M6-PHOENIX-ACES_Jband_vsini1.0_R100k.txt"
    sampled = "resampled_new/Spectrum_M6-PHOENIX-ACES_Jband_vsini1.0_R100k_res3.txt"

    conv_wav, theor_flux, conv_flux = IOmodule.pdread_3col(convolved)
    sampled_wav, sampled_flux = IOmodule.pdread_2col(sampled)

    theor_flux = np.array(theor_flux)
    conv_flux = np.array(conv_flux)

    ratio_flux = moving_average(conv_flux, 300) / moving_average(theor_flux, 300)
    ratio_flux = ratio_flux/ratio_flux[0]

    plt.figure(1)
    plt.xlabel(r"wavelength [$\mu$m])")
    plt.ylabel(r"Flux[ ] ")
    plt.plot(conv_wav, np.array(theor_flux)/theor_flux[0], color='k')
    plt.plot(conv_wav, np.array(conv_flux)/conv_flux[0], color='b')
    plt.plot(conv_wav, ratio_flux, color='g', linestyle='--')
    plt.show()
    plt.close()

    conv_flux_corrected = conv_flux / ratio_flux

    plt.figure(1)
    plt.xlabel(r"wavelength [$\mu$m])")
    plt.ylabel(r"Flux corrected[ ] ")
    plt.plot(conv_wav, np.array(theor_flux)/theor_flux[0], color='k')
    plt.plot(conv_wav, np.array(conv_flux_corrected)/conv_flux_corrected[0], color='b')
    plt.show()
    plt.close()


def RV_cumulative(RV_vector):
    """
    funtion that calculates the cumulative RV vector weighted_error
    """

    return[weighted_error(RV_vector[:2]), weighted_error(RV_vector[:3]), weighted_error(RV_vector[:4]), weighted_error(RV_vector)]


def weighted_error(RV_vector):
    """
    function that calculates the average weighted error from a vector of errors
    """

    RV_vector = np.array(RV_vector)
    RV_value = 1.0/(np.sqrt(np.sum((1.0/RV_vector)**2.0)))

    return RV_value


def moving_average(x, window_size):
    """
    moving average
    """
    window = np.ones(int(window_size))/float(window_size)
    return np.convolve(x, window, 'same')


###############################################################################

if __name__ == "__main__":
    args = vars(_parser())
    opts = {k: args[k] for k in args}
    sys.exit(main(**opts))

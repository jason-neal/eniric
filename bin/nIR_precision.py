"""Near-Infrared radial velocity precision"""
import re
import sys
import argparse
import itertools
import numpy as np
import matplotlib.pyplot as plt

import eniric.IOmodule as IO
import eniric.Qcalculator as Qcalculator
import eniric.utilities as utils
import eniric.atmosphere as atm
import eniric.plotting_functions as plt_functions

from matplotlib import rc
rc('text', usetex=True)   # set stuff for latex usage


def _parser():
    """Take care of all the argparse stuff.

    :returns: the args
    """
    parser = argparse.ArgumentParser(description='Calculate radial velocity precision of model spectra.')

    parser.add_argument("-b", "--bands", type=str, default="J",
                        choices=["ALL", "VIS", "GAP", "Z", "Y", "J", "H", "K", None],
                        help="Wavelength bands to select. Default=J.", nargs="+")
    args = parser.parse_args()
    return args

file_error_to_catch = getattr(__builtins__, 'FileNotFoundError', IOError)

def main(bands="J", plot_bary=False):
    """ Main function that calls calc_precision.

    Parameters
    ----------

    bands: str or list of str or None, Default="J"
        Band letters to use. None does the bands Z through K.
    plot_bary: bool
        Flag to plot and test the barycentric masking then exit.
    """

    resampled_dir = "../data/resampled/"

    spectral_types = ["M0", "M3", "M6", "M9"]
    if isinstance(bands, str):
        bands = [bands]
    elif (bands is None) or (bands is "None"):
        bands = ["Z", "Y", "J", "H", "K"]

    vsini = ["1.0", "5.0", "10.0"]
    resolution = ["60k", "80k", "100k"]
    sampling = ["3"]

    results = calculate_prec(spectral_types, bands, vsini, resolution, sampling,
                             resampled_dir=resampled_dir, plot_bary=plot_bary,
                             plot_atm=False, plot_ste=False, plot_flux=False,
                             paper_plots=False, offset_RV=0.0)

    print("{Combination\t\tPrec_1\t\tPrec_2\t\tPrec_3")
    print("-"*20)
    for key in results:
        print("{0:s}\t\t{1:0.4f}\t{2:0.4f}\t{3:0.4f}".format(key, results[key][0], results[key][1], results[key][2]))
    # Save precision results

    # return results


def strip_result_quantities(results):
    """Remove the units from Quantity results."""
    for key in results:
        results[key] = [results[key][0].value, results[key][1].value, results[key][2].value]
    return results



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
        sys.exit(1)

    return flux_stellar / ((norm_constant / 100.0)**2.0)


def calculate_prec(spectral_types, bands, vsini, resolution, sampling,
                   resampled_dir, plot_bary=False, plot_atm=False,
                   plot_ste=False, plot_flux=True, paper_plots=True,
                   offset_RV=0.0):

    for band in bands:

        atmmodel = "../data/atmmodel/Average_TAPAS_2014_{}.txt".format(band)
        print("Reading atmospheric model...")
        wav_atm, flux_atm, std_flux_atm, mask_atm = prepare_atmopshere(atmmodel)
        print(("There were {0:d} unmasked pixels out of {1:d}., or {2:.1%}."
              "").format(np.sum(mask_atm), len(mask_atm), np.sum(mask_atm) / len(mask_atm)))
        print("The model ranges from {0:4.2f} to {1:4.2f} micron.".format(wav_atm[0], wav_atm[-1]))
        print("Done.")

        print("Calculating impact of Barycentric movement on mask...")

        else:
            mask_atm = barycenter_shift(wav_atm, mask_atm, offset_RV=offset_RV)
        print(("There were {0:d} unmasked pixels out of {1:d}, or {2:.1%}."
               "").format(np.sum(mask_atm), len(mask_atm), np.sum(mask_atm) /
                          len(mask_atm)))

        if plot_atm:
            # moved ploting code to separate code, eniric.plotting_functions.py
            plt_functions.plot_atmopshere_model(wav_atm, flux_atm, mask_atm)

        # theoretical ratios calculation
        # wav_M0, flux_M0, wav_M3, flux_M3, wav_M6, flux_M6, wav_M9, flux_M9 = read_nIRspectra()

        results = {}    # creating empty dictionary for the results
        wav_plot_M0 = []   # creating empty lists for the plots
        flux_plot_M0 = []
        wav_plot_M3 = []
        flux_plot_M3 = []
        wav_plot_M6 = []
        flux_plot_M6 = []
        wav_plot_M9 = []
        flux_plot_M9 = []

        iterations = itertools.product(spectral_types, vsini, resolution, sampling)
        # for star in spectral_types:
        #     for vel in vsini:
        #         for res in resolution:
        #             for smpl in sampling:
        for (star, vel, res, smpl) in iterations:
            file_to_read = ("Spectrum_{0}-PHOENIX-ACES_{1}band_vsini{2}_R{3}"
                            "_res{4}.txt").format(star, band, vel, res, smpl)
            # print("Working on "+file_to_read+".")
            try:
                wav_stellar, flux_stellar = IO.pdread_2col(resampled_dir + file_to_read)
            except file_error_to_catch:
                # Trun list of strings into strings without symbols  ["J", "K"] -> J K
                spectral_str = re.sub(r"[\[\]\"\'\,]", "", str(spectral_types))
                band_str = re.sub(r"[\[\]\"\'\,]", "", str(bands))
                vsini_str = re.sub(r"[\[\]\"\'\,]", "", str(vsini))
                res_str = re.sub(r"[\[\]\"\'\,]", "", str(resolution))
                sampling_str = re.sub(r"[\[\]\"\'\,]", "", str(sampling))

                print(("\nFor just this file I suggest you run\n\tpython nIR_run.py -s {0} -b {1} -v {2} -R {3} "
                       "--sample_rate {4}\nOr for all the combinations you ran here\n\tpython nIR_run.py -s {5}"
                       " -b {6} -v {7} -R {8} --sample_rate {9}"
                       "").format(star, band, vel, res, smpl, spectral_str, band_str, vsini_str, res_str, sampling_str))
                raise
            # Removing boundary effects
            wav_stellar = wav_stellar[2:-2]
            flux_stellar = flux_stellar[2:-2]

            # sample was left aside because only one value existed
            id_string = "{0}-{1}-{2}-{3}".format(star, band, vel, res)

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
            if(id_string in ["M0-J-1.0-100k", "M3-J-1.0-100k",
                             "M6-J-1.0-100k", "M9-J-1.0-100k"]):
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

            # Precision as given by the second_method
            wav_stellar_chunks, flux_stellar_chunks = Qcalculator.bug_fixed_clumping_method(wav_stellar, flux_stellar, mask_atm_selected)

            prec_2_old = Qcalculator.RVprec_calc_masked(wav_stellar_chunks, flux_stellar_chunks)
            prec_2 = Qcalculator.RVprec_calc_masked(wav_stellar, flux_stellar, mask_atm_selected)

            assert np.all(prec_2_old == prec_2)

            """
            # histogram checking
            lengths = [len(chunk) for chunk in flux_stellar_chunks_unformated]
            n, bins, patches = plt.hist(lengths, 500, range=[0.5, 500.5], histtype='stepfilled')
            plt.title(id_string)
            plt.show()
            """

            # Precision as given by the third_method
            prec_3 = Qcalculator.RV_prec_calc_Trans(wav_stellar, flux_stellar, flux_atm_selected)

            # Adding Precision results to the dictionary
            results[id_string] = [prec_1, prec_2, prec_3]

            # Prepare/Do for the ploting.
            if(plot_ste or plot_ste == id_string):
                plt_functions.plot_stellar_spectum(wav_stellar, flux_stellar,
                                                   wav_atm_selected, mask_atm_selected)

            plot_ids = ["M3-Z-1.0-100k", "M3-Y-1.0-100k", "M3-J-1.0-100k",
                        "M3-H-1.0-100k", "M3-K-1.0-100k"]

            if(plot_flux and id_string in plot_ids):
                wav_plot_M0.append(wav_stellar)
                flux_plot_M0.append(flux_stellar)
            if(plot_flux and id_string in plot_ids):
                wav_plot_M3.append(wav_stellar)
                flux_plot_M3.append(flux_stellar)
            if(plot_flux and id_string in plot_ids):
                wav_plot_M6.append(wav_stellar)
                flux_plot_M6.append(flux_stellar)
            if(plot_flux and id_string in plot_ids):
                wav_plot_M9.append(wav_stellar)
                flux_plot_M9.append(flux_stellar)

    if(plot_flux):
        plt_functions.plot_nIR_flux()

    if paper_plots:
        plt_functions.plot_paper_plots()

    else:
        return results


###############################################################################
def compare_output():
    """Function that compares a spectrum prior to convolution, after, and after resampling
    """

    pre_convolution = "PHOENIX_ACES_spectra/lte03900-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes_wave_CUT_nIR.dat"
    pre_wav, pre_flux = IO.pdread_2col(pre_convolution)
    pre_wav = np.array(pre_wav, dtype="float64")*1.0e-4  # conversion to microns
    pre_flux = np.array(pre_flux, dtype="float64")*pre_wav

    convolved = "results_new/Spectrum_M6-PHOENIX-ACES_Jband_vsini1.0_R100k.txt"
    sampled = "resampled_new/Spectrum_M6-PHOENIX-ACES_Jband_vsini1.0_R100k_res3.txt"

    conv_wav, theor_flux, conv_flux = IO.pdread_3col(convolved)
    sampled_wav, sampled_flux = IO.pdread_2col(sampled)

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


def calculate_all_masked(wav_atm, mask_atm):
    """Auxiliary function to calculate masked pixels in banded parts.

    Needs the code to load the atmopsheric data in for each band.

    Need to load a average_tapas file.
    barycenter correct.
    concatenate result.
    """

    # calculating the number of pixels inside the mask
    wav_Z, mask_Z = utils.band_selector(wav_atm, mask_atm, "Z")
    wav_Y, mask_Y = utils.band_selector(wav_atm, mask_atm, "Y")
    wav_J, mask_J = utils.band_selector(wav_atm, mask_atm, "J")
    wav_H, mask_H = utils.band_selector(wav_atm, mask_atm, "H")
    wav_K, mask_K = utils.band_selector(wav_atm, mask_atm, "K")

    bands_masked = np.concatenate((mask_Z, mask_Y, mask_J, mask_H, mask_K))

    print(("Inside the bands, there were {0:.0f} unmasked pixels out of {1:d}"
           ", or {2:.1%}.").format(np.sum(bands_masked), len(bands_masked),
                                   np.sum(bands_masked) / len(bands_masked)))


def RV_cumulative(RV_vector):
    """Function that calculates the cumulative RV vector weighted_error."""

    return [weighted_error(RV_vector[:2]), weighted_error(RV_vector[:3]),
            weighted_error(RV_vector[:4]), weighted_error(RV_vector)]


def weighted_error(RV_vector):
    """Function that calculates the average weighted error from a vector of errors."""

    RV_vector = np.array(RV_vector)
    RV_value = 1.0/(np.sqrt(np.sum((1.0/RV_vector)**2.0)))

    return RV_value


def moving_average(x, window_size):
    """Moving average."""
    window = np.ones(int(window_size))/float(window_size)
    return np.convolve(x, window, 'same')


###############################################################################

if __name__ == "__main__":
    args = vars(_parser())
    opts = {k: args[k] for k in args}
    sys.exit(main(**opts))

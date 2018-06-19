#!/usr/bin/env python
"""Near-Infrared radial velocity precision."""
import argparse
import itertools
import os
import re
import sys

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc

import eniric
import eniric.atmosphere as atm
import eniric.IOmodule as io
import eniric.plotting_functions as plt_functions
import eniric.Qcalculator as Qcalculator
import eniric.snr_normalization as snrnorm
import eniric.utilities as utils

rc("text", usetex=True)  # set stuff for latex usage


def _parser():
    """Take care of all the argparse stuff.

    :returns: the args
    """
    parser = argparse.ArgumentParser(
        description="Calculate radial velocity precision of model spectra."
    )

    parser.add_argument(
        "-b",
        "--bands",
        type=str,
        default="J",
        choices=["ALL", "VIS", "GAP", "Z", "Y", "J", "H", "K", "None"],
        help="Wavelength bands to select. Default=J.",
        nargs="+",
    )
    parser.add_argument(
        "-u",
        "--use_unshifted",
        default=False,
        action="store_true",
        help="Start with the un-doppler-shifted atmmodel.",
    )
    parser.add_argument(
        "-s", "--save", default=False, action="store_true", help="Save results to file."
    )
    parser.add_argument(
        "--snr", help="Mid-band SNR scaling. (Default=100)", default=100, type=float
    )
    parser.add_argument(
        "--ref_band",
        help="SNR reference band. Default=J. (Default=100). "
        "'self' scales each band relative to the SNR itself.",
        choices=["self", "VIS", "GAP", "Z", "Y", "J", "H", "K"],
        default="J",
        type=str,
    )
    return parser.parse_args()


def main(bands="J", use_unshifted=False, save=False, snr=100, ref_band="J"):
    """Main function that calls calc_precision.

    Parameters
    ----------
    bands: str or list of str or None, Default="J"
        Band letters to use. None does the bands Z through K.
    use_unshifted: bool default=False
        Flag to start with the un-Doppler shifted atmmodel.
    save: bool
        Save results to file.

    """
    os.makedirs(eniric.paths["precision"], exist_ok=True)

    spectral_types = ["M0", "M3", "M6", "M9"]
    if ("ALL" in bands) or ("None" in bands):
        bands = ["Z", "Y", "J", "H", "K"]
    elif isinstance(bands, str) and (bands != "None"):
        bands = [bands]

    vsini = ["1.0", "5.0", "10.0"]
    resolution = ["60k", "80k", "100k"]
    sampling = ["3"]

    results = calculate_prec(
        spectral_types,
        bands,
        vsini,
        resolution,
        sampling,
        plot_atm=False,
        plot_ste=False,
        plot_flux=False,
        paper_plots=False,
        rv_offset=0.0,
        use_unshifted=use_unshifted,
        snr=snr,
        ref_band=ref_band,
    )

    print("{Combination\t\tPrec_1\t\tPrec_2\t\tPrec_3")
    print("-" * 20)
    for key in results:
        print(
            "{0:s}\t\t{1:0.4f}\t{2:0.4f}\t{3:0.4f}".format(
                key, results[key][0], results[key][1], results[key][2]
            )
        )
    # Save precision results
    if save:
        output_filename = os.path.join(
            eniric.paths["precision"],
            "precision_results_2017_ref_band-{0}_snr-{1}.txt".format(ref_band, snr),
        )
        ids = []
        prec_1s = []
        prec_2s = []
        prec_3s = []
        for star in spectral_types:
            for band in bands:
                for vel in vsini:
                    for res in resolution:
                        id_string = "{0:s}-{1:s}-{2:.1f}-{3:s}".format(
                            star, band, float(vel), res
                        )
                        ids.append(id_string)
                        prec_1s.append(results[id_string][0].value)
                        prec_2s.append(results[id_string][1].value)
                        prec_3s.append(results[id_string][2].value)

        io.pdwrite_cols(
            output_filename,
            ids,
            prec_1s,
            prec_2s,
            prec_3s,
            header=["# id", r"prec_1 [m/s]", r"prec_2 [m/s]", r"prec_3 [m/s]"],
            float_format="%5.01f",
        )
        print("saved results to {}".format(output_filename))
    # return results


def strip_result_quantities(results):
    """Remove the units from Quantity results."""
    for key in results:
        results[key] = [
            results[key][0].value,
            results[key][1].value,
            results[key][2].value,
        ]
    return results


def calculate_prec(
    spectral_types,
    bands,
    vsini,
    resolution,
    sampling,
    plot_atm=False,
    plot_ste=False,
    plot_flux=True,
    paper_plots=True,
    rv_offset=0.0,
    use_unshifted=False,
    snr=100,
    ref_band="J",
    new=True,
):
    """Calculate precisions for given combinations."""
    # TODO: iterate over band last so that the J band normalization value can be
    # obtained first and applied to each band.

    print("using new config.yaml file here!!!!!!!!!!!!!!")
    results = {}  # creating empty dictionary for the results
    wav_plot_m0 = []  # creating empty lists for the plots
    flux_plot_m0 = []
    wav_plot_m3 = []
    flux_plot_m3 = []
    wav_plot_m6 = []
    flux_plot_m6 = []
    wav_plot_m9 = []
    flux_plot_m9 = []

    for band in bands:

        if use_unshifted:
            atmmodel = os.path.join(
                eniric.paths["atmmodel"], "Average_TAPAS_2014_{}.txt".format(band)
            )
            print("Reading atmospheric model...")
            wav_atm, flux_atm, std_flux_atm, mask_atm = atm.prepare_atmosphere(atmmodel)
            print(
                (
                    "There were {0:d} unmasked pixels out of {1:d}., or {2:.1%}." ""
                ).format(
                    np.sum(mask_atm), len(mask_atm), np.sum(mask_atm) / len(mask_atm)
                )
            )

            print(
                "The model ranges from {0:4.2f} to {1:4.2f} micron.".format(
                    wav_atm[0], wav_atm[-1]
                )
            )
            print("Done.")
            print("Calculating impact of Barycentric movement on mask...")
            # mask_atm = atm.old_barycenter_shift(wav_atm, mask_atm, rv_offset=rv_offset)
            mask_atm = atm.barycenter_shift(wav_atm, mask_atm, rv_offset=rv_offset)
        else:
            shifted_atmmodel = os.path.join(
                eniric.paths["atmmodel"], "Average_TAPAS_2014_{}_bary.txt".format(band)
            )
            print("Reading pre-doppler-shifted atmospheric model...")
            wav_atm, flux_atm, std_flux_atm, mask_atm = atm.prepare_atmosphere(
                shifted_atmmodel
            )
        print("Done.")

        print(
            ("There were {0:d} unmasked pixels out of {1:d}, or {2:.1%}." "").format(
                np.sum(mask_atm), len(mask_atm), np.sum(mask_atm) / len(mask_atm)
            )
        )

        if plot_atm:
            # moved plotting code to separate code, eniric.plotting_functions.py
            plt_functions.plot_atmosphere_model(wav_atm, flux_atm, mask_atm)

        # theoretical ratios calculation
        # wav_m0, flux_m0, wav_m3, flux_m3, wav_m6, flux_m6, wav_m9, flux_m9 = read_nIRspectra()

        iterations = itertools.product(spectral_types, vsini, resolution, sampling)
        # for star in spectral_types:
        #     for vel in vsini:
        #         for res in resolution:
        #             for smpl in sampling:
        for (star, vel, res, smpl) in iterations:
            file_to_read = (
                "Spectrum_{0}-PHOENIX-ACES_{1}band_vsini{2}_R{3}" "_res{4:2.01f}.txt"
            ).format(star, band, vel, res, float(smpl))
            # print("Working on "+file_to_read+".")
            try:
                wav_stellar, flux_stellar = io.pdread_2col(
                    os.path.join(eniric.paths["resampled"], file_to_read)
                )
            except FileNotFoundError:
                # Turn list of strings into strings without symbols  ["J", "K"] -> J K
                spectral_str = re.sub(r"[\[\]\"\',]", "", str(spectral_types))
                band_str = re.sub(r"[\[\]\"\',]", "", str(bands))
                vsini_str = re.sub(r"[\[\]\"\',]", "", str(vsini))
                res_str = re.sub(r"[\[\]\"\',]", "", str(resolution))
                sampling_str = re.sub(r"[\[\]\"\',]", "", str(sampling))

                print(
                    (
                        "\nFor just this file I suggest you run\n\tpython nIR_run.py -s {0} -b {1} -v {2} -R {3} "
                        "--sample_rate {4}\nOr for all the combinations you ran here\n\tpython nIR_run.py -s {5}"
                        " -b {6} -v {7} -R {8} --sample_rate {9}"
                        ""
                    ).format(
                        star,
                        band,
                        vel,
                        res,
                        smpl,
                        spectral_str,
                        band_str,
                        vsini_str,
                        res_str,
                        sampling_str,
                    )
                )
                raise
            # Removing boundary effects
            wav_stellar = wav_stellar[2:-2]
            flux_stellar = flux_stellar[2:-2]

            # sample was left aside because only one value existed
            # TODO: Add metallicity and logg into id string
            id_string = "{0:s}-{1:s}-{2:.1f}-{3:s}".format(star, band, float(vel), res)

            # Getting the wav, flux and mask values from the atm model
            # that are the closest to the stellar wav values, see
            # https://stackoverflow.com/questions/2566412/find-nearest-value-in-numpy-array
            index_atm = np.searchsorted(wav_atm, wav_stellar)
            # replace indexes outside the array, at the very end, by the value at the very end
            # index_atm = [index if(index < len(wav_atm)) else len(wav_atm)-1 for index in index_atm]
            indx_mask = index_atm >= len(wav_atm)  # find broken indexs
            index_atm[indx_mask] = len(wav_atm) - 1  # replace with index of end.

            wav_atm_selected = wav_atm[index_atm]
            flux_atm_selected = flux_atm[index_atm]
            mask_atm_selected = mask_atm[index_atm]

            # Check mask masks out deep atmosphere absorption
            if np.any(flux_atm_selected[mask_atm_selected] < 0.98):
                print(
                    "####WARNGING####\nThis absorption mask does not mask out deep atmosphere transmission!"
                )
                print(
                    "Min flux_atm_selected[mask_atm_selected] = {} < 0.98\n####".format(
                        np.min(flux_atm_selected[mask_atm_selected])
                    )
                )

            # Normalize to SNR 100 in middle of J band 1.25 micron!
            # flux_stellar = normalize_flux(flux_stellar, id_string)
            # flux_stellar = snrnorm.normalize_flux(flux_stellar, id_string, new=True)  # snr=100, ref_band="J"
            flux_stellar = snrnorm.normalize_flux(
                flux_stellar, id_string, new=new, snr=snr, ref_band=ref_band
            )  # snr=100, ref_band="J"

            if id_string in [
                "M0-J-1.0-100k",
                "M3-J-1.0-100k",
                "M6-J-1.0-100k",
                "M9-J-1.0-100k",
            ]:
                index_ref = np.searchsorted(
                    wav_stellar, 1.25
                )  # searching for the index closer to 1.25 micron
                snr_estimate = np.sqrt(
                    np.sum(flux_stellar[index_ref - 1 : index_ref + 2])
                )
                print(
                    "\tSanity Check: The S/N for the {0:s} reference model was of {1:4.2f}.".format(
                        id_string, snr_estimate
                    )
                )
            elif "J" in id_string:
                index_ref = np.searchsorted(
                    wav_stellar, 1.25
                )  # searching for the index closer to 1.25 micron
                snr_estimate = np.sqrt(
                    np.sum(flux_stellar[index_ref - 1 : index_ref + 2])
                )
                print(
                    "\tSanity Check: The S/N for the {0:s} non-reference model was of {1:4.2f}.".format(
                        id_string, snr_estimate
                    )
                )

            # Precision given by the first method:
            print("Performing analysis for: ", id_string)
            prec_1 = Qcalculator.RVprec_calc(wav_stellar, flux_stellar)

            # Precision as given by the second_method
            wav_stellar_chunks, flux_stellar_chunks = Qcalculator.mask_clumping(
                wav_stellar, flux_stellar, mask_atm_selected
            )

            prec_2_old = Qcalculator.RVprec_calc_masked(
                wav_stellar_chunks, flux_stellar_chunks
            )
            prec_2 = Qcalculator.RVprec_calc_masked(
                wav_stellar, flux_stellar, mask_atm_selected
            )

            assert np.all(prec_2_old == prec_2)

            """
            # histogram checking
            lengths = [len(chunk) for chunk in flux_stellar_chunks_unformatted]
            n, bins, patches = plt.hist(lengths, 500, range=[0.5, 500.5], histtype='stepfilled')
            plt.title(id_string)
            plt.show()
            """

            # Precision as given by the third_method
            prec_3 = Qcalculator.RV_prec_calc_Trans(
                wav_stellar, flux_stellar, flux_atm_selected
            )

            # Adding Precision results to the dictionary
            results[id_string] = [prec_1, prec_2, prec_3]

            # Prepare/Do for the plotting.
            if plot_ste or plot_ste == id_string:
                plt_functions.plot_stellar_spectum(
                    wav_stellar, flux_stellar, wav_atm_selected, mask_atm_selected
                )

            plot_ids = [
                "M3-Z-1.0-100k",
                "M3-Y-1.0-100k",
                "M3-J-1.0-100k",
                "M3-H-1.0-100k",
                "M3-K-1.0-100k",
            ]

            if plot_flux and id_string in plot_ids:
                wav_plot_m0.append(wav_stellar)
                flux_plot_m0.append(flux_stellar)
            if plot_flux and id_string in plot_ids:
                wav_plot_m3.append(wav_stellar)
                flux_plot_m3.append(flux_stellar)
            if plot_flux and id_string in plot_ids:
                wav_plot_m6.append(wav_stellar)
                flux_plot_m6.append(flux_stellar)
            if plot_flux and id_string in plot_ids:
                wav_plot_m9.append(wav_stellar)
                flux_plot_m9.append(flux_stellar)

    if plot_flux:
        plt_functions.plot_nIR_flux()

    if paper_plots:
        plt_functions.plot_paper_plots()

    else:
        return results


###############################################################################
def compare_output():
    """Function that compares a spectrum prior to convolution, after, and after resampling."""
    pre_convolution = (
        "PHOENIX_ACES_spectra/lte03900-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes_wave_CUT_nIR.dat"
    )
    pre_wav, pre_flux = io.pdread_2col(pre_convolution)
    pre_wav = np.array(pre_wav, dtype="float64") * 1.0e-4  # conversion to microns
    pre_flux = np.array(pre_flux, dtype="float64") * pre_wav

    convolved = "results_new/Spectrum_M6-PHOENIX-ACES_Jband_vsini1.0_R100k.txt"
    sampled = "resampled_new/Spectrum_M6-PHOENIX-ACES_Jband_vsini1.0_R100k_res3.txt"

    conv_wav, theory_flux, conv_flux = io.pdread_3col(convolved)
    sampled_wav, sampled_flux = io.pdread_2col(sampled)

    theory_flux = np.array(theory_flux)
    conv_flux = np.array(conv_flux)

    ratio_flux = moving_average(conv_flux, 300) / moving_average(theory_flux, 300)
    ratio_flux = ratio_flux / ratio_flux[0]

    plt.figure(1)
    plt.xlabel(r"wavelength [$\mu$m])")
    plt.ylabel(r"Flux[ ] ")
    plt.plot(conv_wav, np.array(theory_flux) / theory_flux[0], color="k")
    plt.plot(conv_wav, np.array(conv_flux) / conv_flux[0], color="b")
    plt.plot(conv_wav, ratio_flux, color="g", linestyle="--")
    plt.show()
    plt.close()

    conv_flux_corrected = conv_flux / ratio_flux

    plt.figure(1)
    plt.xlabel(r"wavelength [$\mu$m])")
    plt.ylabel(r"Flux corrected[ ] ")
    plt.plot(conv_wav, np.array(theory_flux) / theory_flux[0], color="k")
    plt.plot(
        conv_wav, np.array(conv_flux_corrected) / conv_flux_corrected[0], color="b"
    )
    plt.show()
    plt.close()


def calculate_all_masked(wav_atm, mask_atm):
    """Auxiliary function to calculate masked pixels in banded parts.

    Needs the code to load the atmospheric data in for each band.

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

    print(
        (
            "Inside the bands, there were {0:.0f} unmasked pixels out of {1:d}"
            ", or {2:.1%}."
        ).format(
            np.sum(bands_masked),
            len(bands_masked),
            np.sum(bands_masked) / len(bands_masked),
        )
    )


def rv_cumulative(rv_vector):
    """Function that calculates the cumulative RV vector weighted_error."""
    return [
        weighted_error(rv_vector[:2]),
        weighted_error(rv_vector[:3]),
        weighted_error(rv_vector[:4]),
        weighted_error(rv_vector),
    ]


def weighted_error(rv_vector):
    """Function that calculates the average weighted error from a vector of errors."""
    rv_vector = np.array(rv_vector)
    rv_value = 1.0 / (np.sqrt(np.sum((1.0 / rv_vector) ** 2.0)))

    return rv_value


def moving_average(x, window_size):
    """Moving average."""
    window = np.ones(int(window_size)) / float(window_size)
    return np.convolve(x, window, "same")


###############################################################################

if __name__ == "__main__":
    args = vars(_parser())
    opts = {k: args[k] for k in args}
    sys.exit(main(**opts))

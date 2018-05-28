import argparse
import itertools

import numpy as np

from eniric.Qcalculator import quality, RVprec_calc
from eniric.nIRanalysis import convolution
from eniric.resample import log_resample
from eniric.snr_normalization import snr_constant_band
from eniric.utilities import load_aces_spectrum


def _parser():
    """Take care of all the argparse stuff.

    :returns: the args
    """
    parser = argparse.ArgumentParser(description='Calculate quality for any library spectra.')
    parser.add_argument("-t", "--temp", type=int,
                        help="Temperature", nargs="+")
    parser.add_argument("-l", "--logg", type=float, default=[4.5],
                        help="Logg, default = [4.5]", nargs="+")
    parser.add_argument("-m", "--metal", type=float, default=[0.0],
                        help="Metallicity, default=[0.0]", nargs="+")
    parser.add_argument("-a", "--alpha", type=float, default=[0.0],
                        help="Alpha, default=[0.0]", nargs="+")
    parser.add_argument("-s", "--sampling", type=float, default=[3.0],
                        help="Sampling", nargs="+")
    parser.add_argument("-r", "--resolution", type=float, default=[50000],
                        help="Instrumental resolution", nargs="+")
    parser.add_argument("-v", "--vsini", type=float, default=[1.0],
                        help="Rotational Velocity", nargs="+")
    parser.add_argument("-b", "--bands", type=str, default="J",
                        choices=["ALL", "VIS", "GAP", "Z", "Y", "J", "H", "K", 'None'],
                        help="Wavelength bands to select. Default=J.", nargs="+")
    parser.add_argument("--model", type=str, default="phoenix",
                        choices=['phoenix', 'btsettl'],
                        help="Spectral models to use. Default=phoenix.")
    parser.add_argument("--save", default=False, action="store_true",
                        help="Save results to file.")
    # parser.add_argument("--snr", help="Mid-band SNR scaling. (Default=100)", default=100, type=float)
    # parser.add_argument("--ref_band",
    #                    help="SNR reference band. Default=J. (Default=100). "
    #                         "'self' scales each band relative to the SNR itself.",
    #                    choices=["self", "VIS", "GAP", "Z", "Y", "J", "H", "K"], default="J", type=str)
    return parser.parse_args()


def do_analysis(star_params, vsini, R, band, sampling=3.0, conv_kwargs=None):
    """Precision and Quality for specific parameter set."""
    if conv_kwargs is None:
        conv_kwargs = {"epsilon": 0.6, "fwhm_lim": 5.0, "num_procs": 1, "normalize": True}

    # Full photon count spectrum
    wav, flux = load_aces_spectrum(star_params, photons=True)

    # Convolutions
    # This limits the spectrum to the required band
    # wav_band, flux_band, convolved_flux = convolution(wav, flux, vsini, R, band, **conv_kwargs)
    #
    # # Re-sample to sampling per resolution element.
    # wav_grid = log_resample(wav_band, sampling, R)
    #
    # sampled_flux = np.interp(wav_grid, wav_band, convolved_flux)
    wav_grid, sampled_flux = convolve_and_resample(wav, flux, vsini, R, band, sampling, conv_kwargs)

    # Spectral Quality
    q = quality(wav_grid, sampled_flux)

    # Scale normalization for precision
    ref_band = "J"
    snr = 100
    wav_ref, sampled_ref = convolve_and_resample(wav, flux, vsini, R, ref_band, sampling, conv_kwargs)
    snr_normalize = snr_constant_band(wav_ref, sampled_ref, snr=snr, band=ref_band)
    sampled_flux = sampled_flux / snr_normalize

    q2 = quality(wav_grid, sampled_flux)  # After scaling

    if band == "J":
        index_ref = np.searchsorted(wav_grid, 1.25)  # searching for the index closer to 1.25 micron
        snr_estimate = np.sqrt(np.sum(sampled_flux[index_ref - 1:index_ref + 2]))
        print("\tSanity Check: The S/N at 1.25 micron = {0:4.2f}, (should be {1:g}).".format(snr_estimate, snr))

    # Spectral Precision
    prec1 = RVprec_calc(wav_grid, sampled_flux)
    # prec2 =
    # prec3 =

    return [q, q2, prec1]  # ], prec2, prec3]


def convolve_and_resample(wav, flux, vsini, R, band, sampling, conv_kwargs):
    wav_band, flux_band, convolved_flux = convolution(wav, flux, vsini, R, band, **conv_kwargs)
    # Re-sample to sampling per resolution element.
    wav_grid = log_resample(wav_band, sampling, R)
    sampled_flux = np.interp(wav_grid, wav_band, convolved_flux)
    return wav_grid, sampled_flux


if __name__ == "__main__":
    args = _parser()
    print(args.bands)

    try:
        num_procs = args.num_procs
    except AttributeError:
        num_procs = 2
    try:
        normalize = args.normalzie
    except AttributeError:
        normalize = True

    conv_kwargs = {"epsilon": 0.6, "fwhm_lim": 5.0, "num_procs": num_procs, "normalize": normalize}

    # Load the relevant spectra
    models_list = itertools.product(args.temp, args.logg, args.metal, args.alpha)

    with open("Result_text_file.txt", "w") as f:
        f.write("Star Model \t Parameters \t Results\n")
        for model in models_list:
            print("model", model)
            
            # Create generator for params_list
            params_list = itertools.product(args.resolution, args.bands, args.vsini, args.sampling)

            for (R, band, vsini, sample) in params_list:
                print(model, pars, [res.value for res in result])

                f.write("{0}\t{1}\t{2}\n".format(model, pars, [res.value for res in result]))

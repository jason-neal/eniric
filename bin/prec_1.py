"""Quick and dirty precision 1.

Doesn't involve atmopshere model so can perform realively easily to check
precision is working.
"""


import numpy as np
import sys
# import matplotlib.pyplot as plt

# to remove labels in one tick
# from matplotlib.ticker import MaxNLocator

import eniric.IOmodule as IOmodule
import eniric.Qcalculator as Qcalculator

# from eniric.utilities import band_selector

# from eniric.plotting_functions import plot_atmopshere_model, plot_stellar_spectum
from bin.nIR_precision import normalize_flux

# from matplotlib import rc
# set stuff for latex usage
# rc('text', usetex=True)
import argparse


def _parser():
    """Take care of all the argparse stuff.

    :returns: the args
    """
    parser = argparse.ArgumentParser(description='Calculate perfect precision of all convolved spectra.')

    parser.add_argument('-s', '--startype', help='Spectral Type e.g "MO"', type=str, nargs="*", default=None)
    parser.add_argument("-v", "--vsini", help="Rotational velocity of source",
                        type=float, nargs="*", default=None)
    parser.add_argument("-R", "--resolution", help="Observational resolution",
                        type=float, nargs="*", default=None)
    parser.add_argument("-b", "--band", type=str, default="ALL",
                        choices=["ALL", "VIS", "GAP", "Z", "Y", "J", "H", "K"],
                        help="Wavelength band to select", nargs="+")
    parser.add_argument('-d', '--data_dir', help='Data directory', type=str, default=None)
    parser.add_argument('--sample_rate', default=None, type=float, nargs="*",
                        help="Resample rate, pixels per FWHM. Default=3.0")
    parser.add_argument('--results', default=None, type=str,
                        help='Result directory Default=data_dir+"/results/"')
    parser.add_argument('--resamples', default=None, type=str,
                        help='Resample directory. Default=data_dir+"/resampled/"')
    parser.add_argument('--noresample', help='Resample output', default=False,
                        action="store_true")
    parser.add_argument('--normalize', help='Use convolution normalized spectra', default=True,
                        action="store_false")
    parser.add_argument('--org', help='Only use original .dat files, (temporary option)',
                        default=False, action="store_true")
    args = parser.parse_args()
    return args

#atmmodel = "../data/atmmodel/Average_TAPAS_2014.txt"
resampled_dir = "../data/resampled/"
file_error_to_catch = getattr(__builtins__,'FileNotFoundError', IOError)

def calc_prec1(star, band,  vel,  resolution,  smpl, normalize=True, resampled_dir=resampled_dir):
    """ Just caluclate precision for 1st case.

    resolution in short form e.g 100k
    """
    if normalize:
        file_to_read = ("Spectrum_{}-PHOENIX-ACES_{}band_vsini{:.1f}_R{}"
                        "_res{}.txt").format(star, band, vel, resolution,  smpl)
    else:
        file_to_read = ("Spectrum_{}-PHOENIX-ACES_{}band_vsini"
                        "{:.1f}_R{}_unnormalized_res{}.txt").format(star, band, vel,
                                                                resolution,
                                                                smpl)
    # print("Working on " + file_to_read)
    wav_stellar, flux_stellar = IOmodule.pdread_2col(resampled_dir + file_to_read)

    # removing boundary effects
    wav_stellar = wav_stellar[2:-2]
    flux_stellar = flux_stellar[2:-2]

    if normalize:
        id_string = "{}-{}-{}-{}".format(star, band, vel, resolution)   # sample was left aside because only one value existed
    else:
        id_string = "{}-{}-{}-{}-unnorm".format(star, band, vel, resolution)   # sample was left aside because only one value existed

    # Normaize to SNR 100 in middle of J band 1.25 micron!
    flux_stellar = normalize_flux(flux_stellar, id_string)

    if(id_string in ["M0-J-1.0-100k", "M3-J-1.0-100k", "M6-J-1.0-100k", "M9-J-1.0-100k"]):
        index_reference = np.searchsorted(wav_stellar, 1.25)    # searching for the index closer to 1.25 micron
        SN_estimate = np.sqrt(np.sum(flux_stellar[index_reference-1:index_reference+2]))
        print("\tSanity Check: The S/N for the {:s} reference model was of {:4.2f}.".format(id_string, SN_estimate))
    elif("J" in id_string):
        index_reference = np.searchsorted(wav_stellar, 1.25)    # searching for the index closer to 1.25 micron
        SN_estimate = np.sqrt(np.sum(flux_stellar[index_reference-1:index_reference+2]))
        print("\tSanity Check: The S/N for the {:s} non-reference model was of {:4.2f}.".format(id_string, SN_estimate))


    # print("Performing analysis for: ", id_string)
    prec_1 = Qcalculator.RVprec_calc(wav_stellar, flux_stellar)

    # print("{}: \t{}".format(id_string, prec_1))
    return id_string, prec_1


def main(startype=None, vsini=None, resolution=None, band=None, data_dir=None, results=None,
         resamples=None, sample_rate=3.0, noresample=False, normalize=True,
         org=False):
    if data_dir is None:
        data_dir = "../data/"

    if results is None:
        results_dir = data_dir + "results/"
    else:
        results_dir = results

    if resamples is None:
        resampled_dir = data_dir + "resampled/"
    else:
        resampled_dir = resamples

    if startype is None:
        spectral_types = ["M0", "M3", "M6", "M9"]
    else:
        spectral_types = startype
    if band is None:
        bands = ["Z", "Y", "J", "H", "K"]
    else:
        bands = band
    if vsini is None:
        vsini = ["1.0", "5.0", "10.0"]
    if resolution is None:
        resolution = ["60k", "80k", "100k"]
    else:
        resolution = ["{:.0f}k".format(R/1000) for R in resolution]

    if sample_rate is None:
        sampling = ["3"]
    else:
        sampling = sample_rate

    precision = {}  # dict for storing precision value
    for star in spectral_types:
        for band in bands:
            for vel in vsini:
                for R in resolution:
                    for smpl in sampling:
                        for normalize in [True, False]:
                            # print(star, band, vel, R, smpl)
                            try:
                                id_string, prec_1 = calc_prec1(star, band, vel, R, smpl, normalize=normalize)
                                precision[id_string] = prec_1
                            except error_to_catch:
                                pass  # When file not found skip
                            except:
                                print(star, band, vel, R, smpl, "normalized"*normalize, "Failed!")


    print("id_string\t\tprec_1")
    for key in precision:
        print("{:s} \t\t{:02.1f}".format(key, precision[key]))

    return 0


if __name__ == '__main__':
    args = vars(_parser())
    #startype = args.pop("startype")  # positional arguments

    opts = {k: args[k] for k in args}

    sys.exit(main(**opts))

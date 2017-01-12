"""Quick and dirty precision 1.

Doesn't involve atmopshere model so can perform realively easily to check
precision is working.
"""


import numpy as np
# from sys import exit
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

#atmmodel = "../data/atmmodel/Average_TAPAS_2014.txt"
resampled_dir = "../data/resampled/"


def calc_prec1(star, band,  vel,  resolution,  smpl, normalize=True):
    """ Just caluclate precision for 1st cas."""
    if normalize:
        file_to_read = ("Spectrum_{}-PHOENIX-ACES_{}band_vsini{}_R{}"
                        "_res{}.txt").format(star, band, vel, resolution,  smpl)
    else:
        file_to_read = ("Spectrum_{}-PHOENIX-ACES_{}band_vsini"
                        "{}_R{}_unnormalized_res{}.txt").format(star, band, vel,
                                                                resolution,
                                                                smpl)
    # print("Working on "+file_to_read+".")
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

    print("{}: \t{}".format(id_string, prec_1))

spectral_types = ["M0", "M3", "M6", "M9"]
bands = ["Z", "Y", "J", "H", "K"]
vsini = ["1.0", "5.0", "10.0"]
R = ["60k", "80k", "100k"]
sampling = ["3"]

precision = {}  # dict for storing precision value
for star in spectral_types:
    for band in bands:
        for vel in vsini:
            for resolution in R:
                for smpl in sampling:
                    for norm in [False, True]:
                        try:
                            id_string, prec_1 = calc_prec1(star, band, vel, resolution, smpl, norm=norm)
                            if norm:
                                id_string += "-norm"
                            precision[id_string] = prec_1
                        except:
                            pass

print("id_string\t\tprec_1")
for key in precision:
    print("{:s} {:02.2f}".format(key, precision[key]))

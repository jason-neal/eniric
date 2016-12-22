#!/usr/bin/python
# Testing script for nIR analysis
# Run new and old code to test output.S
# Jason Neal
# December 2016
from __future__ import division, print_function
from eniric.nIRanalysis import convolution, resample_allfiles
from eniric.IOmodule import pdread_2col
from eniric.Qcalculator import RVprec_calc
import matplotlib.pyplot as plt
import datetime

spectrum_name = "lte03900-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes_wave.dat"

data_rep = "../data/PHOENIX-ACES_spectra/"
results_dir = "../data/results/"
resampled_dir = "../data/resampled/"

spectrum_path = data_rep + spectrum_name
# Some test parameters
bands = ["GAP", "Y", "J", "H", "K"]  # VIS takes 24 minutes to compute.
#                                    # Compared to 24 seconds for the rest???
band = "Y"
R = 100000
vsini = 1
epsilon = 0.6
FWHM_lim = 5
plot = False
numProcs = 4

for band in bands:
    #print("Starting band", band)
    pass
if False:
# New version
    wav_band, flux_band = convolution(spectrum_path, band, vsini, R, epsilon,
                                          FWHM_lim, plot, numProcs=numProcs,
                                          results_dir="../data/results/unnorm/",
                                          normalize=False)

    resample_allfiles(results_dir="../data/results/unnorm/",
                      resampled_dir="../data/resampled/unnorm/")

    wav_band_norm, flux_band_norm = convolution(spectrum_path, band, vsini, R,
                                                epsilon, FWHM_lim, plot,
                                                numProcs=numProcs,
                                                results_dir="../data/results/norm/",
                                                normalize=True)

    resample_allfiles(results_dir="../data/results/norm/",
                      resampled_dir="../data/resampled/norm/")


    if plot:
        # Plot results together
        plt.plot(wav_band, flux_band, label='Unnormalized')
        plt.plot(wav_band_norm, flux_band_norm, label='Normalized (res only)')
        plt.legend(loc=0)
        plt.show()



# Calculate RVPresision
print("Radial velocity Presision, vsini = {0}, R = {1}".format(vsini, R))
for band in bands:

    # Argument unpacking with *
    norm_prec = RVprec_calc(*pdread_2col("../data/resampled/norm/"
                                         "Spectrum_M0-PHOENIX-ACES_"
                                         "{}band_vsini{}_R{}k_res3.txt".format(band, int(vsini),
                                         int(R / 1000))))

    unnorm_prec = RVprec_calc(*pdread_2col("../data/resampled/unnorm/"
                               "Spectrum_M0-PHOENIX-ACES_"
                               "{}band_vsini{}_R{}k_res3.txt".format(band, int(vsini),
                                                                     int(R / 1000))))

    print("Unnormalized RV_Presision {} band = {}".format(band, unnorm_prec))
    print("Normalized RV_Presision {} band = {}".format(band, norm_prec))

#!/usr/bin/python
# Testing script for nIR analysis
# Run new and old code to test output.S
# Jason Neal
# December 2016
from __future__ import division, print_function
from eniric.nIRanalysis import convolve_spectra, resample_allfiles
from eniric.original_code.nIRanalysis import convolution as old_convolution
from eniric.original_code.nIRanalysis import resample_allfiles as old_resample_allfiles
import matplotlib.pyplot as plt
import datetime

spectrum_name = "lte03900-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes_wave.dat"

data_rep = "../data/PHOENIX-ACES_spectra/"
results_dir = "../data/results/"
resampled_dir = "../data/resampled/"

spectrum_path = data_rep + spectrum_name
# Some test parameters
band = "GAP"
R = 100000
vsini = 1
epsilon = 0.6
FWHM_lim = 5
plot = False
numProcs = 4
do_old = False

# for band in ["GAP", "Y", "J", "K"]:
for band in ["K"]:
    # New version
    start_time = datetime.datetime.now()
    print("Time at start of {} band, {}".format(band, start_time))
    wav_band, flux_band, flux_conv_res = convolve_spectra(spectrum_path, band, vsini, R, epsilon, FWHM_lim, plot, numProcs=numProcs)
    end_time = datetime.datetime.now()
    print("Time at end, ", end_time)
    print("Time to run {} band  convolution = {}".format(band, (end_time-start_time)))

    resample_allfiles()

# for band in ["GAP", "Y", "J", "K"]:
for band in ["K"]:
    # The unchanged version
    if do_old:
        old_start_time = datetime.datetime.now()
        print("Time at start of {} band, {}".format(band, old_start_time))
        old_wav_band, old_flux_conv_res = old_convolution(spectrum_path, band, vsini, R, epsilon, FWHM_lim, plot)  # takes a very long time. good progress indicator though
        old_end_time = datetime.datetime.now()
        print("Time at end,  ", old_end_time)
        print("Time to run old {} band convolution = {}".format(band, (old_end_time-old_start_time)))

        old_resample_allfiles()

if plot:
    # Plot results together
    plt.plot(old_wav_band, old_flux_conv_res, label='Old code')
    plt.plot(wav_band, flux_conv_res, label='New code')
    plt.legend(loc=0)
    plt.show()

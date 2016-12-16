# Testing script for nIR analysis
# Jason Neal
# December 2016
from __future__ import division, print_function
from eniric.nIRanalysis import convolution
from eniric.original_code.nIRanalysis import convolution as old_convolution
import matplotlib.pyplot as plt
import datetime

spectrum_name = "lte03900-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes_wave.dat"

data_rep = "../data/nIRmodels/"
results_dir = "../data/results/"
resampled_dir = "../data/resampled/"

spectrum_path = data_rep+"PHOENIX-ACES/PHOENIX-ACES-AGSS-COND-2011-HiRes/" + spectrum_name
# Some test parameters
band = "Y"
R = 105000
vsini = 1
epsilon = 0.6
FWHM_lim = 5
plot = True
# print("readin =", read_spectrum(spectrum))  # takes a bit of time

# New version
start_time = datetime.datetime.now()
wav_band, flux_conv_res = convolution(spectrum_path, band, vsini, R, epsilon, FWHM_lim, plot)  # takes a very long time. good progress indicator though
end_time = datetime.datetime.now()
print("Time to run new convolution = {}".format((end_time-start_time)))


# The unchanged version
old_start_time = datetime.datetime.now()
old_wav_band, old_flux_conv_res = old_convolution(spectrum_path, band, vsini, R, epsilon, FWHM_lim, plot)  # takes a very long time. good progress indicator though
old_end_time = datetime.datetime.now()
print("Time to run old convolution = {}".format((end_time-start_time)))

# Plot results together
plt.plot(old_wav_band, old_flux_conv_res, label='Old code')
plt.plot(wav_band, flux_conv_res, label='New code')
plt.legend(loc=0)
plt.show()

# Test rotation kernal
# test_rotation_kernal(vsini, epsilon)

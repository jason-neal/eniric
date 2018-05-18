#!/usr/bin/env python
""" Script to time convolution using different number of processors.
# Jason Neal
# December 2016
"""
from __future__ import division, print_function

import datetime
import os

import eniric
from eniric.nIRanalysis import convolve_spectra

spectrum_name = "lte03900-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes_wave.dat"

data_rep = "../data/nIRmodels/"
results_dir = eniric.paths["results"]

spectrum_path = os.path.join(data_rep, "PHOENIX-ACES", "PHOENIX-ACES-AGSS-COND-2011-HiRes", spectrum_name)
# Some test parameters
band = "K"
R = 100000
vsini = 1
epsilon = 0.6
fwhm_lim = 5
plot = False

numprocs = [None, 0, 1, 2, 3, 4]


def time_diff_procs(num_procs):
    """Time the convolution with different number of processors"""
    conv_times = dict()
    for proc in num_procs:
        start_time = datetime.datetime.now()
        convolve_spectra(spectrum_path, band, vsini, R, epsilon, fwhm_lim, plot, num_procs=proc)
        end_time = datetime.datetime.now()
        conv_times[proc] = end_time - start_time
    return conv_times


convolution_times = time_diff_procs(numprocs)

print("Num Processors\t Time")

for key in numprocs:
    print("{0}\t{1}".format(key, convolution_times[key]))

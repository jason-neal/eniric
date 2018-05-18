#!/usr/bin/env python
"""
Script to run all convolutions for the new NIR analysis.
Includes adding metallicity and logg parameters

Jason Neal, January 2017

"""
import itertools

import numpy as np

fehs = np.arange(-2, 1, 0.5)
spec_type = ["M0", "M3", "M6", "M9"]
loggs = np.arange(3.0, 6.5, 0.5)
alphas = [0.0]
spectrum_iter = itertools.product(spec_type, loggs, fehs, alphas)

resolution = [60000, 80000, 100000]
vsini = [1.0, 5.0, 10.0]
sampling = [1, 3, 6]  # effect of sampling on precision

convolution_iter = itertools.product(resolution, vsini)

if __name__ == __main__:
    # do convolutions
    pass

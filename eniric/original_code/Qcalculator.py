"""
Created on Mon Dec 29 00:14:56 2014

@author: pfigueira
"""

from eniric.original_code.IOmodule import read_2col

import numpy as np

c = 299792458  # m/s


def RVprec_calc(spectrum_file="resampled/Spectrum_M0-PHOENIX-ACES_Hband_vsini1.0_R60k_res3.txt"):
    """
    function that claculates the RV precision achievable on a spectrum
    """
    wavelength, flux = read_2col(spectrum_file)

    return [c / SqrtSumWis(wavelength, flux)]


def SqrtSumWis(wavelength, flux):
    """
    Calculation of the SquareRoot of the sum of the Wis for a spectrum
    """

    delta_F = (np.array(flux[1:]) - np.array(flux[:-1]))
    delta_l = np.array(wavelength[:-1])

    derivF_over_lambda = delta_F/delta_l

    return np.sqrt( np.sum( np.array(wavelength[:-1])**2.0 * derivF_over_lambda **2.0 / np.array(flux[:-1]) ) )

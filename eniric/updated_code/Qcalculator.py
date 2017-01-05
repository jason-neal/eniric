"""
Created on Mon Dec 29 00:14:56 2014

@author: pfigueira
"""

from IOmodule import read_2col

import numpy as np

c = 299792458 # m/s

def RVprec_test(spectrum_file= "resampled_new/Spectrum_M0-PHOENIX-ACES_Jband_vsini1.0_R100k_res3.txt"):
    wavelength, flux = read_2col(spectrum_file)
    wavelength = np.array(wavelength, dtype="float64")
    flux = np.array(flux, dtype="float64")/1000.0

    lenght = len(wavelength)
    print("normal_prec:", RVprec_calc(wavelength, flux))
    print("first half prec", RVprec_calc(wavelength[:lenght/2], flux[:lenght/2]))
    print("second half prec", RVprec_calc(wavelength[lenght/2:], flux[lenght/2:]))

    wavs=[]
    flxs=[]
    for i in range(lenght/10, lenght, lenght/10):
        wavs.append(wavelength[(i-lenght/10):i])
        flxs.append(flux[(i-lenght/10): i])

    print("result of division into ten pieces", RVprec_calc_chunks(wavs, flxs))

def RVprec_calc(wavelength, flux):
    """
    function that calculates the RV precision achievable on a spectrum
    """
    return c / SqrtSumWis(wavelength, flux)

def SqrtSumWis(wavelength, flux):
    """
    Calculation of the SquareRoot of the sum of the Wis for a spectrum
    """

    delta_F = (np.array(flux[1:]) - np.array(flux[:-1]))
    delta_l = (np.array(wavelength[1:]) - np.array(wavelength[:-1]))

    derivF_over_lambda = delta_F/delta_l

    return np.sqrt( np.sum( np.array(wavelength[:-1])**2.0 * derivF_over_lambda **2.0 / np.array(flux[:-1])))

def RVprec_calc_chunks(wavelength, flux):
    """
    the same as RVprec_calc, but now wavelength and flux are organized into chunks and the weighted average formula is used
    """

    RV_vector = np.array([RVprec_calc(wav_chunk, flux_chunk) for wav_chunk, flux_chunk in zip(wavelength, flux)])
    RV_value = 1.0/(np.sqrt(np.sum( (1.0/RV_vector)**2.0)))

    return RV_value

###############################################################################

def RV_prec_calc_Trans(wavelength, flux, transmission):
    """
    The same as RV_prec_calc, but considering a transmission different than zero
    """
    return c / SqrtSumWisTrans(wavelength, flux, transmission)

def SqrtSumWisTrans(wavelength, flux, transmission):
    """
    Calculation of the SquareRoot of the sum of the Wis for a spectrum, considering transmission
    """

    delta_F = (np.array(flux[1:]) - np.array(flux[:-1]))
    delta_l = (np.array(wavelength[1:]) - np.array(wavelength[:-1]))

    derivF_over_lambda = delta_F/delta_l

    transmission = np.array(transmission)

    return np.sqrt( np.sum( np.array(wavelength[:-1])**2.0 * derivF_over_lambda **2.0 / (np.array(flux[:-1])/(np.array(transmission[:-1])**2.0))))

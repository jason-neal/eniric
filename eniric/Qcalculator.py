"""
Created on Mon Dec 29 00:14:56 2014

@author: pfigueira

Editied Thur Dec 15 13:00 2016 by Jason Neal for eniric.
"""

# from eniric.IOmodule import read_2col

import numpy as np
import pandas as pd
c = 299792458  # m/s

def RVprec_test(spectrum_file= "resampled/Spectrum_M0-PHOENIX-ACES_Hband_vsini1.0_R60k_res3.txt"):
    data = pd.read_table(spectrum_file, comment='#', header=None, names=["wavelength", "flux"], dtype=np.float64, delim_whitespace=True)
    wavelength, flux = data["wavelength"].values, data["flux"].values

    return RVprec_calc(wavelength, flux)

def RVprec_calc(wavelength, flux):
    """ Calculate the RV precision achievable on a spectrum.
    """

    return c / SqrtSumWis(wavelength, flux)


def SqrtSumWis(wavelength, flux):
    """
    Calculation of the SquareRoot of the sum of the Wis for a spectrum
    """

    delta_F = np.diff(flux)
    delta_l = np.diff(wavelength)

    derivF_over_lambda = delta_F / delta_l

    return np.sqrt(np.sum(wavelength[:-1]**2.0 * derivF_over_lambda**2.0 / flux[:-1]))


def RVprec_calc_chunks(wavelength, flux):
    """
    the same as RVprec_calc, but now wavelength and flux are organized into chunks and the weighted average formula is used
    """

    RV_vector = np.array([RVprec_calc(wav_chunk, flux_chunk) for wav_chunk, flux_chunk in zip(wavelength, flux)])

    RV_value = 1.0/(np.sqrt(np.sum( (1.0/RV_vector)**2.0 )))

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

    Parameters
    ----------
    wavelength: array-like
        Wavelength array
    flux: array-like
        Flux array
    transmission: array-like
        Transmission array

    Returns
    -------
    SqrtSumWisTrans: array-like
        Squarerooted sum of pixel weigths including transmission.
    """

    delta_F = np.diff(flux)
    delta_l = np.diff(wavelength)

    derivF_over_lambda = delta_F/delta_l

    return np.sqrt(np.sum(wavelength[:-1]**2.0 * derivF_over_lambda**2.0 / (flux[:-1])/(transmission[:-1]**2.0))))

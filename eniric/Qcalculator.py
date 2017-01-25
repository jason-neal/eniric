"""
Created on Mon Dec 29 00:14:56 2014

@author: pfigueira

Editied Thur Dec 15 13:00 2016 by Jason Neal for eniric.
"""

# from eniric.IOmodule import read_2col

import numpy as np
import pandas as pd
# c = 299792458  # m/s
from astropy.constants import c

def RVprec_test(spectrum_file= "resampled/Spectrum_M0-PHOENIX-ACES_Hband_vsini1.0_R60k_res3.txt"):
    """Test a RVprec_calc for a singal specturm.
    """
    data = pd.read_table(spectrum_file, comment='#', header=None, names=["wavelength", "flux"], dtype=np.float64, delim_whitespace=True)
    wavelength, flux = data["wavelength"].values, data["flux"].values

    return RVprec_calc(wavelength, flux)


def RVprec_calc(wavelength, flux):
    """ Calculate the RV precision achievable on a spectrum.

    Parameters
    ----------
    wavelength: array-like
        Wavelength of spectrum.
    flux: array-like
        Flux of spectrum.

    Returns
    -------
    RVrms : float
       Radial velocity precision of spectra in m/s.

    Notes
    -----
    Extract from https://arxiv.org/pdf/1511.07468v1.pdf
    From Eq. (11) and (12) of
    https://arxiv.org/pdf/1511.07468v1.pdf#cite.2001A%26A...374..733B
    it follows that the RV uncertainty associated with the information content
    of a given spectra is given by

        RV_rms = c /(Q *Sum{Ne}) = c / (sqrt{sum{W(i)}})

    in which c is the speed of light in vacuum, Q the quality factor of a
    spectrum, and N_e the total number of photoelectrons collected inside the
    wavelength range of interest. However, the precision can only be calculated
    using the concept of optimal pixel weight W(i) for each of the pixels i
    that compose the spectra.

        W(i) = lambda(i)**2 (d'A_0(i) / d'lambda(i))**2 / (A_0(i) + sigma_D**2)

    in which lambda(i) and A_0(i) are the values of each pixel wave-length and
    flux, respectively. The weight will be proportional to the information
    content of the spectrum, given by the derivative of the amplitude, and
    calculated following Connes (1985).
    The denominator of the previous equation is the variance of the flux of
    the pixel A_0(i), depending on the flux value itself and on the
    detector noise sigma_D. In this paper we exclusively consider the high
    signal-to-noise ratio regime, so we can approximate A_0(i) + sigma_D**2
    to A_0(i).
    """

    return c / SqrtSumWis(wavelength, flux)


def SqrtSumWis(wavelength, flux):
    """ Calculation of the SquareRoot of the sum of the weigths(Wis) for a spectrum.

    Parameters
    ----------
    wavelength: array-like
        Wavelength of spectrum.
    flux: array-like
        Flux of spectrum.

    Returns
    -------
    sqrt{sum{W(i)}}: float
       Squareroot of the sum of the pixel weigths(Wis)

    Notes
    -----
    Extract from https://arxiv.org/pdf/1511.07468v1.pdf

    The precision can only be calculated using the concept of optimal pixel weight W(i) for each of the pixels i
    that compose the spectra.

        W(i) = lambda(i)**2 (d'A_0(i) / d'lambda(i))**2 / (A_0(i) + sigma_D**2)

    in which lambda(i) and A_0(i) are the values of each pixel wave-length and
    flux, respectively. The weight will be proportional to the information
    content of the spectrum, given by the derivative of the amplitude, and
    calculated following Connes (1985).

"""
    delta_F = np.diff(flux)
    delta_l = np.diff(wavelength)

    derivF_over_lambda = delta_F / delta_l

    return np.sqrt(np.sum(wavelength[:-1]**2.0 * derivF_over_lambda**2.0 /
                          flux[:-1]))


def RVprec_calc_chunks(wavelength, flux):
    """ The same as RVprec_calc, but now wavelength and flux are organized into
    chunks and the weighted average formula is used

    When considering the average RV as delivered by several slices of a
    spectrum, the error on the average is given by the error on a weighted
    average.
        mean(RV_rms) = 1 / sqrt(sum_i((1 / RV_rms(i))**2))

    Parameters
    ----------
    wavelength: list(array-like)
        Wavelength values of chunks.
    flux: list(array-like)
        Flux values of the chunks.
    Returns
    -------
    RV_value: float
        Weighted average RV value of spectral chunks.

    """


    slice_rvs = np.empty(len(wavelength), dtype=float) # Radial velocity of each slice
    for i, (wav_slice, flux_slice) in enumerate(zip(wavelength, flux)):
        wav_slice = np.array(wav_slice)
        flux_slice = np.array(flux_slice)
        rv_val = RVprec_calc(wav_slice, flux_slice)
        slice_rvs[i] = rv_val # .value   # strip constants unit
    rv_value = 1.0 / (np.sqrt(np.sum((1.0 / slice_rvs)**2.0)))

    return rv_value

###############################################################################

def RV_prec_calc_Trans(wavelength, flux, transmission):
    """
    The same as RV_prec_calc, but considering a transmission different than zero

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
    RVrms: float
        Radial velocity precision for a spectrum affected by atmospheric transmission
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
        Squarerooted sum of pixel weigths including effects of transmission.
    """

    delta_F = np.diff(flux)
    delta_l = np.diff(wavelength)

    derivF_over_lambda = delta_F/delta_l

    return np.sqrt(np.sum(wavelength[:-1]**2.0 * derivF_over_lambda**2.0 /
                          (flux[:-1] / transmission[:-1]**2.0)))

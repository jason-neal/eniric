"""
Created on Mon Dec 29 00:14:56 2014

@author: pfigueira

Editied Thur Dec 15 13:00 2016 by Jason Neal for eniric.
"""

import numpy as np
import pandas as pd
# c = 299792458  # m/s
from astropy.constants import c
import astropy.units as u
from astropy.units import Quantity


def RVprec_test(spectrum_file="resampled/Spectrum_M0-PHOENIX-ACES_Hband_vsini1.0_R60k_res3.txt"):
    """Test a RVprec_calc for a singal specturm.
    """
    data = pd.read_table(spectrum_file, comment='#', header=None,
                         names=["wavelength", "flux"], dtype=np.float64,
                         delim_whitespace=True)
    wavelength, flux = data["wavelength"].values, data["flux"].values

    return RVprec_calc(wavelength, flux)


def RVprec_calc(wavelength, flux):
    """Calculate the RV precision achievable on a spectrum.

    Parameters
    ----------
    wavelength: array-like or Quantity array
        Wavelength of spectrum.
    flux: array-like or Quantity array
        Flux of spectrum.

    Returns
    -------
    RVrms : Quantity scalar
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
    """Calculation of the SquareRoot of the sum of the weigths(Wis) for a spectrum.

    Parameters
    ----------
    wavelength: array-like or Quantity array
        Wavelength of spectrum.
    flux: array-like or Quantity array
        Flux of spectrum.

    Returns
    -------
    sqrt{sum{W(i)}}: float or Quantity scalar
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
    if not isinstance(wavelength, np.ndarray):
        print("Your wavelength and flux should really be numpy arrays! Converting them here.")
        wavelength = np.asarray(wavelength)
    if not isinstance(flux, np.ndarray):
        flux = np.asarray(flux)

    delta_F = np.diff(flux)
    delta_l = np.diff(wavelength)

    derivF_over_lambda = delta_F / delta_l

    if isinstance(flux, u.Quantity):
        """Units of variance are squared """
        flux_variance = flux.value * (flux.unit)**2
    else:
        flux_variance = flux

    return np.sqrt(np.sum(wavelength[:-1]**2.0 * derivF_over_lambda**2.0 /
                          flux_variance[:-1]))


def RVprec_calc_masked(wavelength, flux, mask=None):
    """The same as RVprec_calc, but now wavelength and flux are organized into
    chunks according to the mask and the weighted average formula is used to
    calculate the combined precision.

    When considering the average RV as delivered by several slices of a
    spectrum, the error on the average is given by the error on a weighted
    average.
        mean(RV_rms) = 1 / sqrt(sum_i((1 / RV_rms(i))**2))

    Parameters
    ----------
    wavelength: array-like or list(array-like)
        Wavelength values of chunks.
    flux: array-like or list(array-like)
        Flux values of the chunks.
    mask: array-like of bool or None
        Mask of transmission cuts. Zero values are excluded and used to cut up
        the spectrum.

    Returns
    -------
    RV_value: Quantity scalar
        Weighted average RV value of spectral chunks.


    Notes
    -----
    A "clump" is defined as a contiguous region of the array.
    Solution for clumping comes from
    https://stackoverflow.com/questions/14605734/numpy-split-1d-array-of-chunks-separated-by-nans-into-a-list-of-the-chunks

    There was a bug in the original clumping code which ment that chose the
    clump depending on the first element of mask.
    A test for the failing condition is added so that if ever encountered we
    can investigate the effect on he previously published results.
    """
    if mask is not None:
        if mask[0] is False:  # First value of mask is False was a bug in original code
            print(("{0:s}\nWarning\nA condition that would have given bad "
                   "precision the by broken clumping function was found.\nNeed "
                   "to find the model parameters for this!\n{0:s}\n").format("#"*40))
        # Turn wavelength and flux into masked arrays
        wavelength_clumps, flux_clumps = mask_clumping(wavelength, flux, mask)

    else:
        # When given a already clumped solution and no mask given.
        assert isinstance(wavelength, list)
        assert isinstance(flux, list)
        wavelength_clumps = wavelength
        flux_clumps = flux

    # Turn ndarray into quantity array.
    # Need to use np.zeros instead of np.empty. Unassigned zeros are removed after with nonzero.
    # The "empty" values (1e-300) do not get removed and effect precision
    slice_rvs = Quantity(np.zeros(len(wavelength), dtype=float),
                         unit=u.meter / u.second)  # Radial velocity of each slice

    for i, (wav_slice, flux_slice) in enumerate(zip(wavelength_clumps, flux_clumps)):
        if len(wav_slice) == 1:
            """Results in infinate rv, can not determine the slope of single point."""
            continue

        else:
            wav_slice = np.asarray(wav_slice)
            flux_slice = np.asarray(flux_slice)
            slice_rvs[i] = RVprec_calc(wav_slice, flux_slice)

    # Zeros created from the inital empty array, when skipping single element chunks)
    slice_rvs = slice_rvs[np.nonzero(slice_rvs)]  # Only use nonzero values.

    rv_value = 1.0 / (np.sqrt(np.sum((1.0 / slice_rvs)**2.0)))

    return rv_value


def mask_clumping(wave, flux, mask):
    """Clump contiguous wavelength and flux sections into list.

    Note: Our value of mask (0 = bad points) is opposite to usage in
    np.ma.masked_array (1 = bad)
    Separate function to enable through testing.

    Parameters
    ----------
    wave: array-like of floats
        The wavelength array to clump.
    flux: array-like of floats
        The glux array to clump.
    mask: array-like of bool
        Boolean array with True indicating the values to use/keep.

    Returns
    -------
    wave_clumps: list(array)
        List of the valid wavelength sections.
    flux_clumps: list(array)
       List of the valid flux sections.

    """
    # Turn into masked array to use clump_unmasked method.
    mask = np.asarray(mask, dtype=bool)  # Make it bool so ~ works correctly

    masked_wave = np.ma.masked_array(wave, mask=~mask)   # ma mask is inverted
    masked_flux = np.ma.masked_array(flux, mask=~mask)   # ma mask is inverted

    wave_clumps = [wave[s] for s in np.ma.clump_unmasked(masked_wave)]
    flux_clumps = [flux[s] for s in np.ma.clump_unmasked(masked_flux)]

    return wave_clumps, flux_clumps


def bug_fixed_clumping_method(wav, flux, mask):
    """Old clumping method that is difficult to understand ...[0]+1)[::2].

    There was a signifcant bug which was fixed.
    The returned values were dependant on the first value in the mask. """
    if mask[0] is False:  # First value of mask is False was a bug in original code
        print(("{0:s}\nWarning\nA condition that would have given bad "
               "precision the by broken clumping function was found.\nNeed "
               "to find the model parameters for this!\n{0:s}\n").format("#"*40))

    if mask[0] == 1:
        wav_chunks_unformated = np.array_split(wav, np.where(np.diff(mask))[0]+1)[::2]
        flux_chunks_unformated = np.array_split(flux, np.where(np.diff(mask))[0]+1)[::2]
    else:
        wav_chunks_unformated = np.array_split(wav, np.where(np.diff(mask))[0]+1)[1::2]
        flux_chunks_unformated = np.array_split(flux, np.where(np.diff(mask))[0]+1)[1::2]

    wav_chunks = [list(chunk) for chunk in wav_chunks_unformated]
    flux_chunks = [list(chunk) for chunk in flux_chunks_unformated]

    return wav_chunks, flux_chunks


def bugged_clumping_method(wav, flux, mask):
    """Old clumping method that is difficult to understand ...[0]+1)[::2].
    There was a signifcant bug in which the returned values depend on the first value in mask."""
    wav_chunks_unformated = np.array_split(wav, np.where(np.diff(mask))[0]+1)[::2]
    wav_chunks = [list(chunk) for chunk in wav_chunks_unformated]

    flux_chunks_unformated = np.array_split(flux, np.where(np.diff(mask))[0]+1)[::2]
    flux_chunks = [list(chunk) for chunk in flux_chunks_unformated]

    return wav_chunks, flux_chunks


###############################################################################
def RV_prec_calc_Trans(wavelength, flux, transmission):
    """The same as RV_prec_calc, but considering a transmission different than zero

    Parameters
    ----------
    wavelength: array-like or Quantity array
        Wavelength array
    flux: array-like or Quantity array
        Flux array
    transmission: array-like
        Transmission array

    Returns
    -------
    RVrms: Quantity scalar
        Radial velocity precision for a spectrum affected by atmospheric transmission
    """

    return c / SqrtSumWisTrans(wavelength, flux, transmission)


def SqrtSumWisTrans(wavelength, flux, transmission):
    """Calculation of the SquareRoot of the sum of the Wis for a spectrum, considering transmission.

    The transmission reduces the flux so has an affect on the variance.

    Parameters
    ----------
    wavelength: array-like or Quantity array
        Wavelength array
    flux: array-like or Quantity array
        Flux array
    transmission: array-like
        Transmission array

    Returns
    -------
    SqrtSumWisTrans: array-like or Quantity
        Squarerooted sum of pixel weigths including effects of transmission.
    """

    if not isinstance(wavelength, np.ndarray):
        print("Your wavelength and flux should really be numpy arrays! Converting them here.")
        wavelength = np.asarray(wavelength)
    if not isinstance(flux, np.ndarray):
        flux = np.asarray(flux)
    if not isinstance(transmission, np.ndarray):
        transmission = np.asarray(transmission)

    # Check for units of transmission
    if isinstance(transmission, u.Quantity):
        if not transmission.unit == u.dimensionless_unscaled:
            raise TypeError("transmission has a unit that is not dimensionless and unscaled!")

        # Check for values of quantity transmission
        if np.any(transmission.value > 1) or np.any(transmission.value < 0):
            raise ValueError("Transmission should range from 0 to 1 only.")
    else:
        # Check for values of transmission
        if np.any(transmission > 1) or np.any(transmission < 0):
            raise ValueError("Transmission should range from 0 to 1 only.")

    delta_F = np.diff(flux)
    delta_l = np.diff(wavelength)

    derivF_over_lambda = delta_F / delta_l

    if isinstance(flux, u.Quantity):
        """Units of variance are squared"""
        flux_variance = flux.value * (flux.unit)**2
    else:
        flux_variance = flux

    return np.sqrt(np.sum(wavelength[:-1]**2.0 * derivF_over_lambda**2.0 /
                          (flux_variance[:-1] / transmission[:-1]**2.0)))

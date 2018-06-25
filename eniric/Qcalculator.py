"""Calculate the Radial Velocity Precision of NIR Spectra.

Using the Quality factor of the spectra.
"""
import warnings
from typing import Any, List, Optional, Tuple, Union

import astropy.units as u
import numpy as np
import pandas as pd
from astropy.constants import c
from astropy.units.quantity import Quantity
from numpy import float64, int32, ndarray


def RVprec_test(spectrum_file: str = "resampled/Spectrum_M0-PHOENIX-ACES_Hband_vsini1.0_R60k_res3.txt") -> Quantity:
    """Test a RVprec_calc for a single spectrum."""
    data = pd.read_table(spectrum_file, comment='#', header=None,
                         names=["wavelength", "flux"], dtype=np.float64,
                         delim_whitespace=True)
    wavelength, flux = data["wavelength"].values, data["flux"].values

    return RVprec_calc(wavelength, flux)


def RVprec_calc(wavelength: Union[Quantity, ndarray], flux: Union[Quantity, ndarray]) -> Quantity:
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
    return c / sqrt_sum_wis(wavelength, flux)


def quality(wavelength: Union[Quantity, ndarray], flux: Union[Quantity, ndarray]) -> Union[
    float64, Quantity]:
    """Calculation of the spectral Quality, Q, for a spectrum.

    Parameters
    ----------
    wavelength: array-like or Quantity array
        Wavelength of spectrum.
    flux: array-like or Quantity array
        Flux of spectrum.

    Returns
    -------
    sqrt{sum{W(i)}}: float or Quantity scalar
       Spectral quality

    Notes
    -----
    Extract from https://arxiv.org/pdf/1511.07468v1.pdf

        Q = sqrt{sum{W(i)}} / sqrt{sum{A_0{i}}

    where, W(i), is the optimal pixel weights

        W(i) = lambda(i)**2 (d'A_0(i) / d'lambda(i))**2 / (A_0(i) + sigma_D**2)

    in which lambda(i) and A_0(i) are the values of each pixel wave-length and
    flux, respectively. The weight will be proportional to the information
    content of the spectrum, given by the derivative of the amplitude, and
    calculated following Connes (1985).

    The spectral quality, Q, is independent of the flux level and is only
    a function of the spectral profile.

    """
    if not isinstance(wavelength, np.ndarray):
        print("Your wavelength and flux should really be numpy arrays! Converting them here.")
        wavelength = np.asarray(wavelength)
    if not isinstance(flux, np.ndarray):
        flux = np.asarray(flux)

    flux = flux * u.dimensionless_unscaled # Turn into Quantity if not already
    flux = flux / flux.unit  # Remove units from flux (sqrt(N_e) is unitless)

    wis = sqrt_sum_wis(wavelength, flux)

    return wis / np.sqrt(np.nansum(flux))


def sqrt_sum_wis(wavelength: Union[Quantity, ndarray], flux: Union[Quantity, ndarray]) -> Union[
    float64, Quantity]:
    """Calculation of the Square root of the sum of the weights(Wis) for a spectrum.

    Parameters
    ----------
    wavelength: array-like or Quantity array
        Wavelength of spectrum.
    flux: array-like or Quantity array
        Flux of spectrum.

    Returns
    -------
    sqrt{sum{W(i)}}: float or Quantity scalar
       Square root of the sum of the pixel weights(Wis)

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

    delta_flux = np.diff(flux)
    delta_lambda = np.diff(wavelength)

    derivf_over_lambda = delta_flux / delta_lambda

    if isinstance(flux, u.Quantity):
        """Units of variance are squared"""
        flux_variance = flux.value * flux.unit * flux.unit
    else:
        flux_variance = flux

    wis = np.sqrt(np.nansum(wavelength[:-1] ** 2.0 * derivf_over_lambda ** 2.0 /
                         flux_variance[:-1]))
    if not np.isfinite(wis):
        warnings.warn("Weight sum is not finite = {}".format(wis))
    return wis


def RVprec_calc_masked(wavelength: Union[List[List[Any]], ndarray],
                       flux: Union[ndarray, List[List[Any]]], mask: Optional[ndarray] = None) -> Quantity:
    """RV precision for split apart spectra.

    The same as RVprec_calc, but now wavelength and flux are organized into
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

    There was a bug in the original clumping code which meant that chose the
    clump depending on the first element of mask.
    A test for the failing condition is added so that if ever encountered we
    can investigate the effect on the previously published results.
    """
    if mask is not None:
        # Turn wavelength and flux into masked arrays
        wavelength_clumps, flux_clumps = mask_clumping(wavelength, flux, mask)

    else:
        # When given a already clumped solution and no mask given.
        assert isinstance(wavelength, list)
        assert isinstance(flux, list)
        wavelength_clumps = wavelength
        flux_clumps = flux

    # Turn an ndarray into quantity array.
    # Need to use np.zeros instead of np.empty. Unassigned zeros are removed after with nonzero.
    # The "empty" values (1e-300) do not get removed and effect precision
    slice_rvs = Quantity(np.zeros(len(wavelength), dtype=float),
                         unit=u.meter / u.second)  # Radial velocity of each slice

    for i, (wav_slice, flux_slice) in enumerate(zip(wavelength_clumps, flux_clumps)):
        if len(wav_slice) == 1:
            """Results in infinite rv, can not determine the slope of single point."""
            continue

        else:
            wav_slice = np.asarray(wav_slice)
            flux_slice = np.asarray(flux_slice)
            slice_rvs[i] = RVprec_calc(wav_slice, flux_slice)

    # Zeros created from the initial empty array, when skipping single element chunks)
    slice_rvs = slice_rvs[np.nonzero(slice_rvs)]  # Only use nonzero values.
    rv_value = 1.0 / (np.sqrt(np.nansum((1.0 / slice_rvs) ** 2.0)))

    return rv_value


def mask_clumping(wave: ndarray, flux: ndarray, mask: ndarray) -> Tuple[List[ndarray], List[ndarray]]:
    """Clump contiguous wavelength and flux sections into list.

    Note: Our value of mask (0 = bad points) is opposite to usage in
    np.ma.masked_array (1 = bad)
    Separate function to enable thorough testing.

    Parameters
    ----------
    wave: array-like of floats
        The wavelength array to clump.
    flux: array-like of floats
        The flux array to clump.
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

    masked_wave = np.ma.masked_array(wave, mask=~mask)  # ma mask is inverted
    masked_flux = np.ma.masked_array(flux, mask=~mask)  # ma mask is inverted

    wave_clumps = [wave[s] for s in np.ma.clump_unmasked(masked_wave)]
    flux_clumps = [flux[s] for s in np.ma.clump_unmasked(masked_flux)]

    return wave_clumps, flux_clumps


def bug_fixed_clumping_method(wav: ndarray, flux: ndarray, mask: ndarray) -> Tuple[List[Any], List[Any]]:
    """Old clumping method that is difficult to understand ...[0] + 1)[::2].

    There was a significant bug which was fixed.
    The returned values were dependant on the first value in the mask.
    """
    if mask[0] == 1:
        wav_chunks_unformatted = np.array_split(wav, np.where(np.diff(mask))[0] + 1)[::2]
        flux_chunks_unformatted = np.array_split(flux, np.where(np.diff(mask))[0] + 1)[::2]
    else:
        wav_chunks_unformatted = np.array_split(wav, np.where(np.diff(mask))[0] + 1)[1::2]
        flux_chunks_unformatted = np.array_split(flux, np.where(np.diff(mask))[0] + 1)[1::2]

    wav_chunks = [list(chunk) for chunk in wav_chunks_unformatted]
    flux_chunks = [list(chunk) for chunk in flux_chunks_unformatted]

    return wav_chunks, flux_chunks


###############################################################################
def RV_prec_calc_Trans(wavelength: ndarray, flux: ndarray, transmission: ndarray) -> Quantity:
    """The same as RV_prec_calc, but considering a transmission different than zero.

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
    return c / sqrt_sum_wis_trans(wavelength, flux, transmission)


def sqrt_sum_wis_trans(wavelength: Union[Quantity, ndarray], flux: Union[Quantity, ndarray],
                       transmission: Union[Quantity, ndarray]) -> Union[float64, Quantity]:
    """Calculation of the Square root of the sum of the Weights for a spectrum, considering transmission.

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
    sqrt_sum_wis_trans: array-like or Quantity
        Square root sum of pixel weights including effects of transmission.

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

    delta_flux = np.diff(flux)
    delta_lambda = np.diff(wavelength)

    derivf_over_lambda = delta_flux / delta_lambda

    if isinstance(flux, u.Quantity):
        """Units of variance are squared"""
        flux_variance = flux.value * flux.unit * flux.unit
    else:
        flux_variance = flux

    return np.sqrt(np.nansum(wavelength[:-1] ** 2.0 * derivf_over_lambda ** 2.0 /
                          (flux_variance[:-1] / transmission[:-1] ** 2.0)))

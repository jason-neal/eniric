"""Legacy functions kept for accessing changes to precision."""
from typing import Any, List, Optional, Tuple, Union

import numpy as np
from astropy import units as u
from astropy.units import Quantity
from numpy.core.multiarray import ndarray

from eniric.Qcalculator import rv_precision


# Legacy function in which the spectra is split into chunks first and then the pixel weights are calculated.
def RVprec_calc_masked(
    wavelength: Union[List[List[Any]], ndarray],
    flux: Union[ndarray, List[List[Any]]],
    mask: Optional[ndarray] = None,
    **kwargs,
) -> Quantity:
    """RV precision for split apart spectra.

    The same as rv_precision, but now wavelength and flux are organized into
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
    kwargs:
        Kwargs for sqrt_sum_wis

    Returns
    -------
    RV_value: Quantity scalar
        Weighted average RV value of spectral chunks.

    Notes
    -----
    A "clump" is defined as a contiguous region of the array.
    Solution for clumping comes from
    https://stackoverflow.com/questions/14605734/numpy-split-1d-array-of-chunks-separated-by-nans-into-a-list-of-the-chunks
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
    slice_rvs = Quantity(
        np.zeros(len(wavelength_clumps), dtype=float), unit=u.meter / u.second
    )  # Radial velocity of each slice

    for i, (wav_slice, flux_slice) in enumerate(zip(wavelength_clumps, flux_clumps)):
        if len(wav_slice) == 1:
            # Results in infinite rv, can not determine the slope of single point.
            continue

        else:
            wav_slice = np.asarray(wav_slice)
            flux_slice = np.asarray(flux_slice)
            slice_rvs[i] = rv_precision(wav_slice, flux_slice, **kwargs)

    # Zeros created from the initial empty array, when skipping single element chunks)
    slice_rvs = slice_rvs[np.nonzero(slice_rvs)]  # Only use nonzero values.
    rv_value = 1.0 / (np.sqrt(np.nansum((1.0 / slice_rvs) ** 2.0)))

    return rv_value


def mask_clumping(
    wave: ndarray, flux: ndarray, mask: ndarray
) -> Tuple[List[ndarray], List[ndarray]]:
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


def RVprec_calc_weights_masked(
    wavelength: ndarray, flux: ndarray, mask: Optional[ndarray] = None, **kwargs
) -> Quantity:
    """RV precision setting weights of telluric lines to zero.

    Instead of splitting the spectra after every telluric line and
    individually calculating precision and taking the weighted average
    this just sets the pixel weights to zero.

    Parameters
    ----------
    wavelength: ndarray
        Wavelength array
    flux: ndarray
        Flux array
    mask: array or None
        Mask of transmission cuts. Zero values are excluded and used to cut up
        the spectrum.

    Returns
    -------
    RV_value: Quantity scalar
        RV precision.

    RVprec_calc_weights_masked is just a wrapper around rv_precision.

    """
    return rv_precision(wavelength, flux, mask=mask, **kwargs)

"""Calculate the Radial Velocity Precision of NIR Spectra.

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

    Spectral Quality:
        Q = sqrt{sum{W(i)}} / sqrt{sum{A_0{i}}

    The spectral quality, Q, is independent of the flux level and is only
    a function of the spectral profile.

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
import warnings
from typing import Any, List, Optional, Tuple, Union

import astropy.units as u
import numpy as np
import astropy.constants as const
from astropy.units.quantity import Quantity
from numpy import float64, ndarray

c = const.c


def RVprec_calc(
    wavelength: Union[Quantity, ndarray], flux: Union[Quantity, ndarray], **kwargs
) -> Quantity:
    """Calculate the theoretical RV precision achievable on a spectrum.

    Parameters
    ----------
    wavelength: array-like or Quantity array
        Wavelength of spectrum.
    flux: array-like or Quantity array
        Flux of spectrum.
    kwargs:
        Kwargs for sqrt_sum_wis
    Returns
    -------
    RVrms : Quantity scalar
       Radial velocity precision of spectra in m/s.

    """
    return c / sqrt_sum_wis(wavelength, flux, **kwargs)


def quality(
    wavelength: Union[Quantity, ndarray], flux: Union[Quantity, ndarray], **kwargs
) -> Union[float64, Quantity]:
    """Calculation of the spectral Quality, Q, for a spectrum.

    Parameters
    ----------
    wavelength: array-like or Quantity array
        Wavelength of spectrum.
    flux: array-like or Quantity array
        Flux of spectrum.
    kwargs:
        Kwargs for sqrt_sum_wis

    Returns
    -------
    sqrt{sum{W(i)}}: float or Quantity scalar
       Spectral quality

    Notes
    -----
        Q = sqrt{sum{W(i)}} / sqrt{sum{A_0{i}}

        The spectral quality, Q, is independent of the flux level and is only
    a function of the spectral profile.

    """
    if not isinstance(wavelength, np.ndarray):
        print(
            "Your wavelength and flux should really be numpy arrays! Converting them here."
        )
        wavelength = np.asarray(wavelength)
    if not isinstance(flux, np.ndarray):
        flux = np.asarray(flux)

    flux = flux * u.dimensionless_unscaled  # Turn into Quantity if not already
    flux = flux / flux.unit  # Remove units from flux (sqrt(N_e) is unitless)

    wis = sqrt_sum_wis(wavelength, flux, **kwargs)

    return wis / np.sqrt(np.nansum(flux))


def RVprec_calc_Trans(
    wavelength: ndarray, flux: ndarray, transmission: ndarray, **kwargs
) -> Quantity:
    """The same as RV_prec_calc, but considering a transmission different than zero.

    Parameters
    ----------
    wavelength: array-like or Quantity array
        Wavelength array
    flux: array-like or Quantity array
        Flux array
    transmission: array-like
        Transmission array
    kwargs:
        Kwargs for sqrt_sum_wis

    Returns
    -------
    RVrms: Quantity scalar
        Radial velocity precision for a spectrum affected by atmospheric transmission

    """
    return c / sqrt_sum_wis(wavelength, flux, mask=transmission, **kwargs)


def sqrt_sum_wis(
    wavelength: Union[Quantity, ndarray],
    flux: Union[Quantity, ndarray],
    mask: Optional[Union[Quantity, ndarray]] = None,
    grad: bool = True,
) -> Union[float64, Quantity]:
    """Calculation of the Square root of the sum of the weights(Wis) for a spectrum.

    Mask is used to apply a masking function to the weights (to mask out telluric lines for example)

        W(i) = W(i) * m(i)

    Parameters
    ----------
    wavelength: array-like or Quantity array
        Wavelength of spectrum.
    flux: array-like or Quantity array
        Flux of spectrum.
    mask: Optional ndarray
        Weighting mask function. Default is all ones.
    grad: bool
        Flag to use np.gradient. Default=True. Original publication used less precise method.

    Returns
    -------
    sqrt{sum{W(i)}}: float or Quantity scalar
       Square root of the sum of the pixel weights(Wis)

    """
    if mask is None:
        # Don't use np.ones_like() as it will take units of a Quantity.
        mask = np.ones(len(wavelength))

    mask_check(mask)

    # Square mask for telluric masking.
    # Boolean mask is not affected by square 0->0, 1->1.
    mask = mask ** 2
    pixel_wis = pixel_weights(wavelength, flux, grad=grad)

    # Apply masking function
    if grad:
        masked_wis = pixel_wis * mask
    else:
        masked_wis = pixel_wis * mask[:-1]

    sqrt_sum = np.sqrt(np.nansum(masked_wis))
    if not np.isfinite(sqrt_sum):
        warnings.warn("Weight sum is not finite = {}".format(sqrt_sum))
    return sqrt_sum


def mask_check(mask):
    """Checks for mask array."""
    if isinstance(mask, u.Quantity):
        if not (mask.unit == u.dimensionless_unscaled):
            raise TypeError(
                "Mask should not be a non-dimensionless and unscaled Quantity!"
            )
        mask = mask.value

    # Check for values of mask
    if np.any(mask > 1) or np.any(mask < 0):
        raise ValueError("Mask should be within range from 0 to 1 only.")


def RVprec_calc_masked(
    wavelength: Union[List[List[Any]], ndarray],
    flux: Union[ndarray, List[List[Any]]],
    mask: Optional[ndarray] = None,
    **kwargs
) -> Quantity:
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
            """Results in infinite rv, can not determine the slope of single point."""
            continue

        else:
            wav_slice = np.asarray(wav_slice)
            flux_slice = np.asarray(flux_slice)
            slice_rvs[i] = RVprec_calc(wav_slice, flux_slice, **kwargs)

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

    """
    if mask is None:
        # Turn wavelength and flux into masked arrays
        mask = np.ones_like(wavelength)

    if (len(mask) != len(wavelength)) or (len(wavelength) != len(flux)):
        raise ValueError("Input values are not correct length")

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
        np.zeros(len(wavelength), dtype=float), unit=u.meter / u.second
    )  # Radial velocity of each slice

    for i, (wav_slice, flux_slice) in enumerate(zip(wavelength_clumps, flux_clumps)):
        if len(wav_slice) == 1:
            """Results in infinite rv, can not determine the slope of single point."""
            continue

        else:
            wav_slice = np.asarray(wav_slice)
            flux_slice = np.asarray(flux_slice)
            slice_rvs[i] = RVprec_calc(wav_slice, flux_slice, **kwargs)

    # Zeros created from the initial empty array, when skipping single element chunks)
    slice_rvs = slice_rvs[np.nonzero(slice_rvs)]  # Only use nonzero values.
    rv_value = 1.0 / (np.sqrt(np.nansum((1.0 / slice_rvs) ** 2.0)))

    return rv_value


def slope(wavelength, flux):
    """Finite difference derivative which looses one value of array.

        f' = (f(x+h)-f(x)) / h
    """
    delta_flux = np.diff(flux)
    delta_lambda = np.diff(wavelength)

    return delta_flux / delta_lambda


def pixel_weights(
    wavelength: Union[ndarray, Quantity],
    flux: Union[ndarray, Quantity],
    grad: bool = True,
):
    """Calculate individual pixel weights.
    w(i) = \lambda(i)^2 (\partial A(i)/\partial\lambda)^2 / A(i)

    Parameters
    ----------
    grad: bool
        Toggle function for spectral slope. Default False + forward finite difference.
    """
    if isinstance(flux, u.Quantity):
        """Units of variance are squared"""
        flux_variance = flux.value * flux.unit * flux.unit
    else:
        flux_variance = flux

    dydx_unit = 1
    if grad:
        # Hack for quantities with numpy gradient
        if isinstance(flux, Quantity):
            dydx_unit *= flux.unit
            flux = flux.value
        if isinstance(wavelength, Quantity):
            dydx_unit /= wavelength.unit
            wave = wavelength.value
        else:
            wave = wavelength
        derivf_over_lambda = np.gradient(flux, wave) * dydx_unit

        return (wavelength * derivf_over_lambda) ** 2.0 / flux_variance
    else:
        derivf_over_lambda = slope(wavelength, flux)
        return (wavelength[:-1] * derivf_over_lambda) ** 2.0 / flux_variance[:-1]

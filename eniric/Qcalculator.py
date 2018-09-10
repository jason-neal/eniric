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

    The precision can only be calculated using the concept of
    optimal pixel weight W(i) for each of the pixels i that compose the spectra.

        W(i) = lambda(i)**2 (d'A_0(i) / d'lambda(i))**2 / (A_0(i) + sigma_D**2)

    in which lambda(i) and A_0(i) are the values of each pixel wave-length and
    flux, respectively. The weight will be proportional to the information
    content of the spectrum, given by the derivative of the amplitude, and
    calculated following Connes (1985).

"""
import warnings
from typing import Optional, Tuple, Union

import astropy.constants as const
import astropy.units as u
import numpy as np
from astropy.units.quantity import Quantity
from numpy import ndarray

from eniric.resample import log_chunks

c = const.c


def rv_precision(
    wavelength: Union[Quantity, ndarray],
    flux: Union[Quantity, ndarray],
    mask: Optional[ndarray] = None,
    **kwargs,
) -> Quantity:
    """Calculate the theoretical RV precision achievable on a spectrum.

    Parameters
    ----------
    wavelength: array-like or Quantity array
        Wavelength of spectrum.
    flux: array-like or Quantity array
        Flux of spectrum.
    mask:Optional[ndarray]
        Masking function array to apply to the pixel weights.
    kwargs:
        Kwargs for sqrt_sum_wis

    Returns
    -------
    RVrms : Quantity scalar
       Radial velocity precision of spectra in m/s.

    """
    return c / sqrt_sum_wis(wavelength, flux, mask=mask, **kwargs)


def quality(
    wavelength: Union[Quantity, ndarray], flux: Union[Quantity, ndarray], **kwargs
) -> Union[float, Quantity]:
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
    sqrt{sum{W(i)}}: float
       Spectral quality

    Notes
    -----
        Q = sqrt{sum{W(i)}} / sqrt{sum{A_0{i}}

    The spectral quality, Q, is independent of the flux level and is only
    a function of the spectral profile.

    """
    flux = flux * u.dimensionless_unscaled  # Turn into Quantity if not already
    flux = flux / flux.unit  # Remove units from flux (sqrt(N_e) is unitless)

    wis = sqrt_sum_wis(wavelength, flux, **kwargs)
    q = wis / np.sqrt(np.nansum(flux))
    return q.value


def sqrt_sum_wis(
    wavelength: Union[Quantity, ndarray],
    flux: Union[Quantity, ndarray],
    mask: Optional[Union[Quantity, ndarray]] = None,
    grad: bool = True,
) -> Union[float, Quantity]:
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

    pixel_wis = pixel_weights(wavelength, flux, grad=grad)

    # Apply masking function
    if grad:
        masked_wis = pixel_wis * mask
    else:
        masked_wis = pixel_wis * mask[:-1]

    sqrt_sum = np.sqrt(np.nansum(masked_wis))
    if not np.isfinite(sqrt_sum):
        warnings.warn("Weight sum is not finite = {}".format(sqrt_sum))
    if sqrt_sum == 0:
        warnings.warn(
            "Sum of weights sum is = {}. This will cause infinite errors.".format(
                sqrt_sum
            )
        )
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


def slope(wavelength, flux):
    """Forward Finite difference derivative which looses one value of array.

        f' = (f(x+h)-f(x)) / h.
        f'[i] = (flux[i+1] - flux[i])/ (wavelength[i+1] - wavelength[i])

    Returns
    -------
    Array with n-1 points.
    """
    return np.diff(flux) / np.diff(wavelength)


def pixel_weights(
    wavelength: Union[ndarray, Quantity],
    flux: Union[ndarray, Quantity],
    grad: bool = True,
):
    """Calculate individual pixel weights.
    w(i) = \lambda(i)^2 (\partial A(i)/\partial\lambda)^2 / A(i)

    Parameters
    ----------
    wavelength: Union[ndarray, Quantity]
        Wavelength array.
    flux: Union[ndarray, Quantity]
     Flux array.
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


def incremental_quality(
    wavelength: ndarray, flux: ndarray, percent: Union[int, float] = 10, **kwargs
) -> Tuple[ndarray, ndarray]:
    """Determine spectral quality in incremental sections.

    Parameters
    ----------
    wavelength: array-like or Quantity array
        Wavelength of spectrum.
    flux: array-like or Quantity array
        Flux of spectrum.
    percent: Union[int,float]  (default=10)
        The percent size of chunk around each wavelength position.
    kwargs:
        Kwargs passed onto quality().

    Returns
    -------
    x: ndarray
     Central wavelength value for quality section.
    q: ndarray
       Spectral quality for each section.
    """
    positions = log_chunks(wavelength, percent)
    qualities = []
    for pos1, pos2 in zip(positions[:-1], positions[1:]):
        mask = (wavelength >= pos1) & (wavelength < pos2)
        x = wavelength[mask]
        y = flux[mask]
        q = quality(x, y, **kwargs)
        qualities.append([np.mean(x), q])
    x, q = np.asarray(qualities).T
    return x, q


def incremental_rv(
    wavelength: ndarray, flux: ndarray, mask: ndarray, percent: float = 10, **kwargs
) -> Tuple[ndarray, ndarray]:
    """Determine spectral RV precision in incremental sections.
    Parameters
    ----------
    wavelength: array-like or Quantity array
        Wavelength of spectrum.
    flux: array-like or Quantity array
        Flux of spectrum.
    flux: array-like or Quantity array
        Pixel weight mask.
    percent: Union[int,float]  (default=10)
        The percent size of chunk around each wavelength position.
    kwargs:
        Kwargs passed onto quality().

    Returns
    -------
    x: ndarray
     Central wavelength value for quality section.
    rv: ndarray
       Spectral RV precision for each section.
    """
    positions = log_chunks(wavelength, percent)
    velocities = []
    for pos1, pos2 in zip(positions[:-1], positions[1:]):
        pos_mask = (wavelength >= pos1) & (wavelength < pos2)
        x = wavelength[pos_mask]
        y = flux[pos_mask]
        if mask is not None:
            z = mask[pos_mask]
        else:
            z = mask  # None
        rv_calc = rv_precision(x, y, z, **kwargs)
        velocities.append([np.mean(x), rv_calc.value])

    x, rv = np.asarray(velocities).T
    return x, rv

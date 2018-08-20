from typing import Optional, Union

import astropy.units as u
import numpy as np
from astropy import constants as const
from astropy.units.quantity import Quantity
from numpy import ndarray

c = const.c


class RvPrecision(object):
    """

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

    def __init__(
        self,
        wavelength: Union[Quantity, ndarray],
        flux: Union[Quantity, ndarray],
        mask: Optional[ndarray] = None,
        telluric: Optional[ndarray] = None,
        grad=True,
    ):
        if mask is None:
            mask = np.ones_like(flux, dtype=np.bool)
        if telluric is None:
            telluric = np.ones_like(flux, dtype=np.bool)

        if (
            len(wavelength) != len(flux)
            or (len(flux) != len(mask))
            or (len(flux) != len(telluric))
        ):
            raise ValueError("Input values are not the correct")

        self.wavelength = wavelength
        self.flux = flux
        self.mask = mask
        self.telluric = telluric
        self.grad = grad

    def pixel_weights(self, mask=None):
        """Optimal pixel weights

        W(i) = lambda(i)**2 (d'A_0(i) / d'lambda(i))**2 / (A_0(i) + sigma_D**2)

        in which lambda(i) and A_0(i) are the values of each pixel wavelength and
        flux, respectively. The weight will be proportional to the information
        content of the spectrum, given by the derivative of the amplitude, and
        calculated following Connes (1985).

        This handle all three cases from Figueira et al. 2016.
        For the case of transmission scaling the variance division is equivalent to multiplication.

        For the telluric masking the mask must be boolian 0, 1 only (the square does not affect it.

        For full inclusion the mask  must be all ones.

        If no mask then
        """
        if mask is None:
            mask = np.ones_like(self.wavelength)

        delta_flux = np.diff(self.flux)
        delta_lambda = np.diff(self.wavelength)

        derivf_over_lambda = delta_flux / delta_lambda

        if isinstance(self.flux, u.Quantity):
            """Units of variance are squared"""
            flux_variance = self.flux.value * (self.flux.unit ** 2)
        else:
            flux_variance = self.flux

        weights = (
            self.wavelength[:-1] ** 2.0 * derivf_over_lambda ** 2.0 / flux_variance[:-1]
        )
        return weights * mask[:-1] ** 2

    def sqrt_sum_weights(self, mask=None):
        """Square root of sum of pixel weights."""
        return np.sqrt(np.nansum(self.pixel_weights(mask=mask)))

    def total_flux(self):
        """Sum of flux.

        Count of photoelectrons N_{e^{-}}}."""
        return np.nansum(self.flux)

    def Q(self):
        """Calculate spectral quality, which is dependant on spectral profile only.

        Q = sqrt{sum{W(i)}} / sqrt{sum{A_0{i}}.
        """
        # TODO
        print("Check if the spectra quality need mask function applied.")
        return self.sqrt_sum_weights(mask=None) / np.sqrt(self.total_flux())

    def condition1(self):
        """Precision using pixel weights without masking."""
        return c / self.sqrt_sum_weights(mask=None)

    def condition2(self):
        return c / self.sqrt_sum_weights(mask=self.mask)

    def condition3(self):
        return c / self.sqrt_sum_weights(mask=self.telluric)

    @staticmethod
    def weighted_error(rv_errors):
        """Function that calculates the average weighted error from a vector of errors."""
        # rv_errors = np.asarray(rv_errors)
        rv_value = 1.0 / (np.sqrt(np.sum((1.0 / rv_errors) ** 2.0)))

        return rv_value

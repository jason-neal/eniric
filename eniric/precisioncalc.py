import warnings

from typing import Any, List, Optional, Tuple, Union

import astropy.units as u
import numpy as np
from astropy.constants import c
from astropy.units.quantity import Quantity
from numpy import float64, ndarray


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

    def __init__(self,
                 wavelength: Union[Quantity, ndarray],
                 flux: Union[Quantity, ndarray],
                 mask: Optional[ndarray] = None,
                 telluric: Optional[ndarray] = None,
                 ):
        if mask is None:
            mask = np.ones_like(flux, dtype=np.bool)
        if telluric is None:
            telluric = np.ones_like(flux, dtype=np.bool)

        if len(wavelength) != len(flux) or (len(flux) != len(mask)) or (len(flux) != len(telluric)):
            raise ValueError("Input values are not the correct")

        self.wavelength = wavelength
        self.flux = flux
        self.mask = mask
        self.telluric = telluric

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

        weigths = self.wavelength[:-1] ** 2.0 * derivf_over_lambda ** 2.0 / flux_variance[:-1]
        return weigths * mask[:-1] ** 2

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

    def rv_prec(self):
        """Precision using pixel weights without masking."""
        return c / self.sqrt_sum_weights(mask=None)

    def rv_prec_mask(self):
        """Using mask directly on pixel weights."""
        # masks should be 1/or zero only
        assert np.all((self.mask == 0) | (self.mask == 1))  # all values are 1 or zero
        #tweaked mask   account for increased pixel loss due to multiple diffs when clumping
        mask = self.mask.copy()
        print("mask ratio", "{}/{}".format(sum(mask),len(mask)))
        mask[:-1] = mask[:-1]*mask[1:]
        assert not np.allclose(self.mask, mask)
        print("new mask ratio", "{}/{}".format(sum(mask),len(mask)))
        return c / self.sqrt_sum_weights(mask=self.mask)

    def rv_prec_clumped(self):
        # return c / self.sqrt_sum_weights(mask=False)
        from eniric.Qcalculator import mask_clumping
        wavelength_clumps, flux_clumps = mask_clumping(self.wavelength, self.flux, self.mask)

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
                slice_rvprec = RvPrecision(wav_slice, flux_slice)
                slice_rvs[i] = slice_rvprec.rv_prec()

        # Zeros created from the initial empty array, when skipping single element chunks)
        slice_rvs = slice_rvs[np.nonzero(slice_rvs)]  # Only use nonzero values.
        slice_rvs = slice_rvs[np.isfinite(slice_rvs)]
        return self.weighted_error(slice_rvs)

    def rv_prec_trans(self):
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
        return c / self.sqrt_sum_weights(mask=self.telluric)


    @staticmethod
    def weighted_error(rv_errors):
        """Function that calculates the average weighted error from a vector of errors."""
        # rv_errors = np.asarray(rv_errors)
        rv_value = 1.0 / (np.sqrt(np.sum((1.0 / rv_errors) ** 2.0)))

        return rv_value

    def rv_shift_mask(self, rv=30):
        """Extend transmission mask.
        RV shift the telluric spectra +- 30km/s"""
        #See how it is done in atmosphere
        # TODO
        # assert False

    def mask_deep_telluric(self, percent=2):
        """Mask out telluric lines deeper than given percent.

        Parameters
        ----------
        percent: float  (default 2)
            Depth percentage limit of telluric lines to mask.

        Transforms self.mask from a transmission spectrum to a boolean mask.
        Idempotent if already a boolean mask.
        """
        if np.all(self.telluric == 1):
            warnings.warn("Telluric us all unity, have you added a telluric spectra.")

        self.mask = self.telluric >= (1 - percent / 100.0)

"""Spectral broadening functions:

Used to convolve the model spectra for
    - Stellar rotation broadening.
    - Instrumental profile broadening.

Uses joblib.Memory to cache convolution results to skip repeated computation.
"""
from typing import Optional, Union

import joblib
import numpy as np
from astropy import constants as const
from joblib import Memory, Parallel, delayed
from numpy.core.multiarray import ndarray
from tqdm import tqdm

from eniric import config
from eniric.utilities import band_selector, cpu_minus_one, mask_between, wav_selector

# Cache convolution results.
memory = Memory(location=config.cache["location"], verbose=0)

c_kmps = const.c.to("km/s").value

num_cpu_minus_1 = cpu_minus_one()


@memory.cache(ignore=["num_procs", "verbose"])
def rotational_convolution(
    wavelength: ndarray,
    extended_wav: ndarray,
    extended_flux: ndarray,
    vsini: float,
    *,
    epsilon: float = 0.6,
    normalize: bool = True,
    num_procs: Optional[Union[int, joblib.parallel.Parallel]] = None,
    verbose: bool = True,
) -> ndarray:
    r"""Perform Rotational convolution.

    Parameters
    ----------
    wavelength: ndarray
        Wavelength array.
    extended_wav: ndarray
       Extended wavelength array to avoid boundary issues.
    extended_flux: ndarray
        Photon flux for the extended wavelength array.
    vsini: float
        Rotational velocity in km/s.
    epsilon: float
        Limb darkening coefficient. Default is 0.6.
    normalize: bool
        Area normalize the broadening kernel. This corrects for kernel area with unequal wavelength spacing.
        Default is True.
    num_procs: int, None or joblib.parallel.Parallel.
        Number of processes to use, n_job parameter in joblib.
        If num_procs =  1, then a single core is used.
        Can also be a joblib.parallel.Parallel instance. Default is None.
    verbose: bool
        Show the tqdm progress bar. Default is True.

    Returns
    -------
    convolved_flux: ndarray
        The convolved flux evaluated at the wavelength array.

    """

    def element_rot_convolution(single_wav: float) -> float:
        r"""Embarrassingly parallel part of rotational convolution.

        Calculates the convolution value for a single pixel.
        The parameters extended_wav, extended_flux, vsini, epsilon and normalize
        are obtained from the outer scope.

        Parameters
        ----------
        single_wav: float
            Wavelength to calculate the convolution at.

        Returns
        -------
        sum_val: float
            Sum of flux convolved for this wavelength.

        """
        # Select all values such that they are within the fwhm limits
        delta_lambda_l = single_wav * vsini / c_kmps

        index_mask = mask_between(
            extended_wav, single_wav - delta_lambda_l, single_wav + delta_lambda_l
        )

        flux_2convolve = extended_flux[index_mask]
        rotation_profile = rotation_kernel(
            extended_wav[index_mask] - single_wav,
            delta_lambda_l,
            vsini=vsini,
            epsilon=epsilon,
        )

        sum_val = np.sum(rotation_profile * flux_2convolve)

        if normalize:
            # Correct for the effect of non-equidistant sampling
            unitary_rot_val = np.sum(rotation_profile)
            sum_val /= unitary_rot_val

        return sum_val

    if num_procs is None:
        num_procs = num_cpu_minus_1

    if vsini != 0:
        tqdm_wav = tqdm(wavelength, disable=not verbose)

        if isinstance(num_procs, int):
            with Parallel(n_jobs=num_procs) as parallel:
                convolved_flux = np.asarray(
                    parallel(delayed(element_rot_convolution)(wav) for wav in tqdm_wav)
                )
        else:
            try:
                # Assume num_procs is a joblib.parallel.Parallel.
                convolved_flux = np.asarray(
                    num_procs(delayed(element_rot_convolution)(wav) for wav in tqdm_wav)
                )
            except TypeError:
                raise TypeError(
                    "num_proc must be an int or joblib.parallel.Parallel. Not '{}'".format(
                        type(num_procs)
                    )
                )
    else:
        # Skip convolution for vsini=0
        if wavelength is extended_wav:
            convolved_flux = extended_flux  # No change
        else:
            # Interpolate to the new wavelength vector.
            convolved_flux = np.interp(wavelength, extended_wav, extended_flux)
    return convolved_flux


@memory.cache(ignore=["num_procs", "verbose"])
def resolution_convolution(
    wavelength: ndarray,
    extended_wav: ndarray,
    extended_flux: ndarray,
    R: float,
    *,
    fwhm_lim: float = 5.0,
    normalize: bool = True,
    num_procs: Optional[Union[int, joblib.parallel.Parallel]] = None,
    verbose: bool = True,
) -> ndarray:
    r"""Perform Resolution convolution.

    Parameters
    ----------
    wavelength: ndarray
        Wavelength array.
    extended_wav: ndarray
       Extended wavelength array to avoid boundary issues.
    extended_flux: ndarray
        Photon flux for the extended wavelength array.
    R: float
        Resolution of Guassian instrumental profile.
    fwhm_lim: float
        FWHM limit for instrument broadening. Default is 5.0.
    normalize: bool
        Area normalize the broadening kernel. This corrects for kernel area with unequal wavelength spacing.
        Default is True.
    num_procs: int, None or joblib.parallel.Parallel.
        Number of processes to use, n_job parameter in joblib.
        If num_procs =  1, then a single core is used.
        Can also be a joblib.parallel.Parallel instance.
    verbose: bool
        Show the tqdm progress bar. Default is True.

    Returns
    -------
    convolved_flux: ndarray
        The convolved flux evaluated at the wavelength array.
    """

    def element_res_convolution(single_wav: float) -> float:
        r"""Embarrassingly parallel component of resolution convolution.

        Calculates the convolution value for a single pixel.

        The parameters extended_wav, fwhm_lim, R and normalize
        are obtained from the outer scope.

        Parameters
        ----------
        single_wav: float
            Wavelength value to calculate convolution at.

        Returns
        -------
        sum_val: float
            Sum of flux convolved for this pixel/wavelength.
        """
        fwhm = single_wav / R
        # Mask of wavelength range within fwhm_lim* fwhm of wav
        fwhm_space = fwhm_lim * fwhm
        index_mask = mask_between(
            extended_wav, single_wav - fwhm_space, single_wav + fwhm_space
        )

        flux_2convolve = extended_flux[index_mask]
        # Gaussian Instrument Profile for given resolution and wavelength
        instrument_profile = unitary_gaussian(
            extended_wav[index_mask], single_wav, fwhm=fwhm
        )

        sum_val = np.sum(instrument_profile * flux_2convolve)
        if normalize:
            # Correct for the effect of convolution with non-equidistant positions
            unitary_val = np.sum(instrument_profile)
            sum_val /= unitary_val

        return sum_val

    tqdm_wav = tqdm(wavelength, disable=not verbose)

    if num_procs is None:
        num_procs = num_cpu_minus_1

    if isinstance(num_procs, int):
        with Parallel(n_jobs=num_procs) as parallel:
            convolved_flux = np.asarray(
                parallel(delayed(element_res_convolution)(wav) for wav in tqdm_wav)
            )
    else:
        # Assume num_procs is joblib.parallel.Parallel.
        try:
            convolved_flux = np.asarray(
                num_procs(delayed(element_res_convolution)(wav) for wav in tqdm_wav)
            )
        except TypeError:
            raise TypeError(
                "num_proc must be an int or joblib.parallel.Parallel. Not '{}'".format(
                    type(num_procs)
                )
            )
    return convolved_flux


@memory.cache(ignore=["num_procs", "verbose"])
def convolution(
    wav: ndarray,
    flux: ndarray,
    vsini: float,
    R: float,
    band: str = "All",
    *,
    epsilon: float = 0.6,
    fwhm_lim: float = 5.0,
    num_procs: Optional[Union[int, joblib.parallel.Parallel]] = None,
    normalize: bool = True,
    verbose: bool = True,
):
    r"""Perform rotational then Instrumental broadening with convolutions.

    Parameters
    ----------
    wav: ndarray
        Wavelength array.
    flux: ndarray
        Flux array.
    vsini: float
        Rotational velocity in km/s.
    R: int
        Resolution of instrumental profile.
    band: str
        Wavelength band to choose, default is "All".
    epsilon: float
        Limb darkening coefficient. Default is 0.6.
    fwhm_lim: float
        FWHM limit for instrument broadening. Default is 5.0.
    normalize: bool
        Area normalize the broadening kernel. This corrects for kernel area with unequal wavelength spacing.
        Default is True.
    num_procs: int, None or joblib.parallel.Parallel.
        Number of processes to use, n_job parameter in joblib.
        If num_procs =  1, then a single core is used.
        Can also be a joblib.parallel.Parallel instance.
    verbose: bool
        Show the tqdm progress bar. Default is True.

    Returns
    -------
    wav_band: ndarray
        Wavelength for the selected band.
    flux_band: ndarray
        Original flux for the selected band.
    flux_conv: ndarray
        Convolved flux for the selected band.
    """

    wav_band, flux_band = band_selector(wav, flux, band)

    # Calculate FWHM at each end for the convolution
    fwhm_min = wav_band[0] / R  # fwhm at the extremes of vector
    fwhm_max = wav_band[-1] / R

    # performing convolution with rotation kernel
    if verbose:
        print("Starting the Rotation convolution for vsini={0:.2f}...".format(vsini))

    delta_lambda_min = wav_band[0] * vsini / c_kmps
    delta_lambda_max = wav_band[-1] * vsini / c_kmps

    # widest wavelength bin for the rotation convolution
    lower_lim = wav_band[0] - delta_lambda_min - fwhm_lim * fwhm_min
    upper_lim = wav_band[-1] + delta_lambda_max + fwhm_lim * fwhm_max
    wav_ext_rotation, flux_ext_rotation = wav_selector(wav, flux, lower_lim, upper_lim)

    # wide wavelength bin for the resolution_convolution
    lower_lim = wav_band[0] - fwhm_lim * fwhm_min
    upper_lim = wav_band[-1] + fwhm_lim * fwhm_max
    extended_wav, __ = wav_selector(wav, flux, lower_lim, upper_lim)

    # rotational convolution
    flux_conv_rot = rotational_convolution(
        extended_wav,
        wav_ext_rotation,
        flux_ext_rotation,
        vsini,
        epsilon=epsilon,
        num_procs=num_procs,
        normalize=normalize,
        verbose=verbose,
    )
    if verbose:
        print("Starting the Resolution convolution...")

    flux_conv_res = resolution_convolution(
        wav_band,
        extended_wav,
        flux_conv_rot,
        R,
        fwhm_lim=fwhm_lim,
        num_procs=num_procs,
        normalize=normalize,
        verbose=verbose,
    )

    return wav_band, flux_band, flux_conv_res


def unitary_gaussian(
    x: Union[range, int, ndarray],
    center: Union[float, int, str],
    fwhm: Union[float, int, str],
) -> ndarray:
    """Gaussian kernal of area 1.

    Parameters
    ----------
    x: array-like
        Position array.
    center: float
        Central position of Gaussian.
    fwhm: float
        Full Width at Half Maximum.

    Returns
    -------
    kernel: array-like
        Gaussian kernel.
    """
    if not isinstance(fwhm, (np.float, np.int)):
        raise TypeError("The fwhm value is not a number, {0}".format(type(fwhm)))
    if not isinstance(center, (np.float, np.int)):
        raise TypeError("The center value is not a number, {0}".format(type(center)))
    if not isinstance(x, np.ndarray):
        raise TypeError

    sigma = np.abs(fwhm) / (2 * np.sqrt(2 * np.log(2)))
    amp = 1.0 / (sigma * np.sqrt(2 * np.pi))
    tau = -((x - center) ** 2) / (2 * (sigma ** 2))
    kernel = amp * np.exp(tau)

    return kernel


def rotation_kernel(
    delta_lambdas: ndarray, delta_lambda_l: float, vsini: float, epsilon: float
) -> ndarray:
    r"""Rotation kernel for a given line center.

    Parameters
    ----------
    delta_lambdas: array
        Wavelength difference from the line center lambda.
    delta_lambda_l: float
        Maximum line shift of line center by vsini.
    vsini: float
        Projected rotational velocity in km/s.
    epsilon: float
        Linear limb-darkening coefficient, between 0 and 1.

    Returns
    -------
    kernel: array
        Rotational kernel.

    Notes
    -----
    Gray, D. F. (2005). The Observation and Analysis of Stellar Photospheres. 3rd ed. Cambridge University Press.

    """
    denominator = np.pi * vsini * (1.0 - epsilon / 3.0)
    lambda_ratio_sqr = (delta_lambdas / delta_lambda_l) ** 2.0

    c1 = 2.0 * (1.0 - epsilon) / denominator
    c2 = 0.5 * np.pi * epsilon / denominator
    kernel = c1 * np.sqrt(1.0 - lambda_ratio_sqr) + c2 * (1.0 - lambda_ratio_sqr)

    return kernel


def oned_circle_kernel(x: ndarray, center: float, fwhm: float):
    """Calculate the convolution kernel for a circular fiber.

    Parameters
    ----------
    x: array
        Value to evaluate kernel at.
    center: float
        Center of kernel.
    fwhm: float
        FWHM of desired kernel.

    Returns
    -------
    Collapsed circle kernel.

    Notes
    -----
    Tries to represent the broadening by the fiber of a fiber feed spectrograph.

    Artigau 2018 - stated this is mathematically equivalent to a cosine between -pi/2 and pi/2.
    """
    fwhm_scale = 2.094_395_1  # Numerically derived

    A = 1  # Amplitude
    B = fwhm_scale / fwhm  # Scale to give specific fwhm

    kernel = A * np.cos(B * (x - center))

    # Limit to main cosine lobe only
    upper_xi = center + np.pi / 2 / B
    lower_xi = center - np.pi / 2 / B
    mask = mask_between(x, lower_xi, upper_xi)
    kernel[~mask] = 0

    return kernel

"""Broadening functions

Used to convolve the spectra for
    - stellar rotation
    - instrumental resolution

"""
from typing import Optional, Union

import multiprocess as mprocess
import numpy as np
from joblib import Memory
from numpy.core.multiarray import ndarray
from tqdm import tqdm

import eniric
from eniric.utilities import band_selector, mask_between, wav_selector

# Cache convolution results.
memory = Memory(cachedir=eniric.cache["cachedir"], verbose=0)


@memory.cache(ignore=["num_procs"])
def rotational_convolution(
    wav_extended,
    wav_ext_rotation,
    flux_ext_rotation,
    vsini,
    epsilon,
    num_procs: Optional[int] = None,
    normalize: bool = True,
):
    """Perform Rotational convolution part of convolution.
    """

    def wrapper_rot_parallel_convolution(args):
        """Wrapper for rot_parallel_convolution needed to unpack the arguments for
        fast_convolve as multiprocess.Pool.map does not accept multiple
        arguments
        """
        return element_rot_convolution(*args)

    def element_rot_convolution(
        wav,
        wav_extended,
        wav_ext_rotation,
        flux_ext_rotation,
        vsini: float,
        epsilon: float,
        normalize: bool,
    ):
        """Embarrassingly parallel part of rotational convolution"""
        # select all values such that they are within the fwhm limits
        delta_lambda_l = wav * vsini / 3.0e5

        index_mask = mask_between(
            wav_ext_rotation, wav - delta_lambda_l, wav + delta_lambda_l
        )

        flux_2convolve = flux_ext_rotation[index_mask]
        rotation_profile = rotation_kernel(
            wav_ext_rotation[index_mask] - wav, delta_lambda_l, vsini, epsilon
        )

        sum_val = np.sum(rotation_profile * flux_2convolve)

        if normalize:
            # Correct for the effect of non-equidistant sampling
            unitary_rot_val = np.sum(rotation_profile)  # Affects precision
            return sum_val / unitary_rot_val
        else:
            return sum_val

    if num_procs != 0:
        if num_procs is None:
            num_procs = mprocess.cpu_count() - 1

        mproc_pool = mprocess.Pool(processes=num_procs)

        args_generator = tqdm(
            [
                [
                    wav,
                    wav_extended,
                    wav_ext_rotation,
                    flux_ext_rotation,
                    vsini,
                    epsilon,
                    normalize,
                ]
                for wav in wav_extended
            ]
        )

        flux_conv_rot = np.array(
            mproc_pool.map(wrapper_rot_parallel_convolution, args_generator)
        )

        mproc_pool.close()

    else:  # num_procs == 0
        flux_conv_rot = np.empty_like(wav_extended)  # Memory assignment
        for ii, wav in enumerate(tqdm(wav_extended)):
            flux_conv_rot[ii] = element_rot_convolution(
                wav,
                wav_extended,
                wav_ext_rotation,
                flux_ext_rotation,
                vsini,
                epsilon,
                normalize=normalize,
            )
        print("Done.\n")
    return flux_conv_rot


@memory.cache(ignore=["num_procs"])
def resolution_convolution(
    wav_band,
    wav_extended,
    flux_conv_rot,
    R,
    fwhm_lim,
    num_procs: Optional[int] = 1,
    normalize: bool = True,
):
    """Perform Resolution convolution part of convolution."""

    # Define inner convolution functions
    def element_res_convolution(
        wav, R, wav_extended, flux_conv_rot, fwhm_lim, normalize: bool = True
    ):
        """Embarrassingly parallel component of resolution convolution"""
        fwhm = wav / R
        # Mask of wavelength range within fwhm_lim* fwhm of wav
        fwhm_space = fwhm_lim * fwhm
        index_mask = mask_between(wav_extended, wav - fwhm_space, wav + fwhm_space)

        flux_2convolve = flux_conv_rot[index_mask]
        # Gaussian Instrument Profile for given resolution and wavelength
        IP = unitary_gaussian(wav_extended[index_mask], wav, fwhm)

        sum_val = np.sum(IP * flux_2convolve)
        if normalize:
            # Correct for the effect of convolution with non-equidistant positions
            unitary_val = np.sum(IP)  # Affects precision
            return sum_val / unitary_val
        else:
            return sum_val

    def wrapper_res_parallel_convolution(args):
        """Wrapper for res_parallel_convolution needed to unpack the arguments
        for fast_convolve as multiprocess.Pool.map does not accept multiple
        arguments
        """
        return element_res_convolution(*args)

    if num_procs != 0:
        if num_procs is None:
            num_procs = mprocess.cpu_count() - 1

        mproc_pool = mprocess.Pool(processes=num_procs)
        # Need to update the values here
        args_generator = tqdm(
            [
                [wav, R, wav_extended, flux_conv_rot, fwhm_lim, normalize]
                for wav in wav_band
            ]
        )
        flux_conv_res = np.array(
            mproc_pool.map(wrapper_res_parallel_convolution, args_generator)
        )
        mproc_pool.close()

    else:  # num_procs == 0
        flux_conv_res = np.empty_like(wav_band)  # Memory assignment
        for jj, wav in enumerate(tqdm(wav_band)):
            flux_conv_res[jj] = element_res_convolution(
                wav, R, wav_extended, flux_conv_rot, fwhm_lim, normalize=normalize
            )
        print("Done.\n")
    return flux_conv_res


@memory.cache(ignore=["num_procs"])
def convolution(
    wav,
    flux,
    vsini,
    R,
    band: str = "All",
    epsilon: float = 0.6,
    fwhm_lim: float = 5.0,
    num_procs: Optional[int] = None,
    normalize: bool = True,
):
    """Perform convolution of spectrum.

    Rotational convolution followed by a Gaussian of a specified resolution R.

    Parameters
    ----------
    wav: ndarray
        Wavelength in microns
    flux: ndarray
        Photon flux
    vsini: float
        Rotational velocity in km/s.
    R: int
        Resolution of instrumental profile.
    band: str
        Wavelength band to choose, default="All"
    num_procs: int, None
        Number of processes to use with multiprocess. If None it is assigned to 1 less then total number of cores.
        If num_procs = 0, then multiprocess is not used.

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

    # We need to calculate the fwhm at this value in order to set the starting point for the convolution
    fwhm_min = wav_band[0] / R  # fwhm at the extremes of vector
    fwhm_max = wav_band[-1] / R

    # performing convolution with rotation kernel
    print("Starting the Rotation convolution for vsini={0:.2f}...".format(vsini))

    delta_lambda_min = wav_band[0] * vsini / 3.0e5
    delta_lambda_max = wav_band[-1] * vsini / 3.0e5

    # widest wavelength bin for the rotation convolution
    lower_lim = wav_band[0] - delta_lambda_min - fwhm_lim * fwhm_min
    upper_lim = wav_band[-1] + delta_lambda_max + fwhm_lim * fwhm_max
    wav_ext_rotation, flux_ext_rotation = wav_selector(wav, flux, lower_lim, upper_lim)

    # wide wavelength bin for the resolution_convolution
    lower_lim = wav_band[0] - fwhm_lim * fwhm_min
    upper_lim = wav_band[-1] + fwhm_lim * fwhm_max
    wav_extended, flux_extended = wav_selector(wav, flux, lower_lim, upper_lim)

    # rotational convolution
    flux_conv_rot = rotational_convolution(
        wav_extended,
        wav_ext_rotation,
        flux_ext_rotation,
        vsini,
        epsilon,
        num_procs=num_procs,
        normalize=normalize,
    )

    print("Starting the Resolution convolution...")

    flux_conv_res = resolution_convolution(
        wav_band,
        wav_extended,
        flux_conv_rot,
        R,
        fwhm_lim,
        num_procs=num_procs,
        normalize=normalize,
    )

    return wav_band, flux_band, flux_conv_res


def unitary_gaussian(
    x: Union[range, int, ndarray],
    center: Union[float, int, str],
    fwhm: Union[float, int, str],
) -> ndarray:
    """Gaussian function of area = 1.

    Parameters
    ----------
    x: array-like
        Position array
    center: float
        Central position of Gaussian
    fwhm: float
        Full Width at Half Maximum

    Returns
    -------
    result: array-like
        Result of gaussian function sampled at x values.
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
    result = amp * np.exp(tau)

    return result


def rotation_kernel(
    delta_lambdas: ndarray, delta_lambda_l: float, vsini: float, epsilon: float
) -> ndarray:
    """Calculate the rotation kernel for a given wavelength

    Parameters
    ----------
    delta_lambdas: array
        Wavelength values selected within delta_lambda_l around central value. (check)
    delta_lambda_l: float
        FWHM of rotational broadening. (check)
    vsini: float
        Projected rotational velocity [km/s]
    epsilon: float
        Linear limb-darkening coefficient (0-1).

    Returns
    -------
        Rotational kernel

    Notes:
    Equations * from .... book.

    """
    denominator = np.pi * vsini * (1.0 - epsilon / 3.0)
    lambda_ratio_sqr = (delta_lambdas / delta_lambda_l) ** 2.0

    c1 = 2.0 * (1.0 - epsilon) / denominator
    c2 = 0.5 * np.pi * epsilon / denominator

    return c1 * np.sqrt(1.0 - lambda_ratio_sqr) + c2 * (1.0 - lambda_ratio_sqr)


def oned_circle_kernel(x, center, fwhm):
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
        Collapsed circle kernel

    Notes:
    Tries to represent the broadening by the fiber of a fiber feed spectrograph.

    Artigau 2018 - stated mathematically equivalent to a cosine between -pi/2 and pi/2. This is what has tried to be created.
    """
    fwhm_scale = 2.0943951  # Numerically derived

    A = 1  # Amplitude
    B = fwhm_scale / fwhm  # Scale to give specific fwhm

    result = A * np.cos(B * (x - center))

    # Limit to main cos lobe only
    upper_xi = center + np.pi / 2 / B
    lower_xi = center - np.pi / 2 / B
    mask = mask_between(x, lower_xi, upper_xi)
    result[~mask] = 0

    return result

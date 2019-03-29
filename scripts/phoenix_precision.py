#!/usr/bin/env python
"""
phoenix_precision.py
--------------------
Script to generate RV precision of synthetic spectra, see :ref:`Calculating-Precisions`.

"""
import argparse
import itertools
import os
import warnings
from typing import Any, List, Optional, Tuple

import numpy as np
from astropy import units as u
from astropy.units import Quantity
from numpy import ndarray

from eniric import config
from eniric.atmosphere import Atmosphere
from eniric.broaden import convolution
from eniric.corrections import correct_artigau_2018
from eniric.precision import quality, rv_precision
from eniric.resample import log_resample
from eniric.snr_normalization import snr_constant_band
from eniric.utilities import (
    band_middle,
    cpu_minus_one,
    doppler_shift_flux,
    load_aces_spectrum,
    load_btsettl_spectrum,
)

num_cpu_minus_1 = cpu_minus_one()

ref_choices = ["SELF"]
ref_choices.extend(config.bands["all"])


def _parser():
    """Take care of all the argparse stuff."""
    parser = argparse.ArgumentParser(
        description="Calculate precision and quality for synthetic spectra."
    )
    parser.add_argument(
        "-t",
        "--temp",
        type=int,
        default=[3900],
        help="Temperature, default=[3900].",
        nargs="+",
    )
    parser.add_argument(
        "-l",
        "--logg",
        type=float,
        default=[4.5],
        help="Logg, default = [4.5].",
        nargs="+",
    )
    parser.add_argument(
        "-m",
        "--metal",
        type=float,
        default=[0.0],
        help="Metallicity, default=[0.0].",
        nargs="+",
    )
    parser.add_argument(
        "-a",
        "--alpha",
        type=float,
        default=[0.0],
        help="Alpha, default=[0.0].",
        nargs="+",
    )
    parser.add_argument(
        "-s",
        "--sampling",
        type=float,
        default=[3.0],
        help="Sampling, default=[3.0].",
        nargs="+",
    )
    parser.add_argument(
        "-r",
        "--resolution",
        type=float,
        default=[50000],
        help="Instrumental resolution, default=[50000]",
        nargs="+",
    )
    parser.add_argument(
        "-v",
        "--vsini",
        type=float,
        default=[1.0],
        help="Rotational Velocity, default = [1.0]",
        nargs="+",
    )
    parser.add_argument(
        "-b",
        "--bands",
        type=str,
        default=["J"],
        choices=config.bands["all"],
        help="Wavelength bands to select, default=['J'].",
        nargs="+",
    )
    parser.add_argument(
        "--model",
        type=str,
        default="aces",
        choices=["aces", "btsettl", "phoenix"],
        help="Spectral models to use, default='aces'.",
    )
    parser.add_argument(
        "--snr", help="Mid-band SNR scaling, default=100.", default=100, type=float
    )
    parser.add_argument(
        "--ref_band",
        help="SNR reference band, default='J'. "
        "'self' scales each band relative to the SNR itself.",
        choices=ref_choices,
        default="J",
        type=str,
    )
    parser.add_argument(
        "--num_procs",
        help="Number of processors to use, default = (Total cores - 1)",
        default=num_cpu_minus_1,
        type=int,
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Filename for result file, default='precisions.csv'.",
        default="precisions.csv",
        type=str,
    )
    parser.add_argument(
        "--rv",
        help="Radial velocity value, default=[0.0]",
        default=[0.0],
        type=float,
        nargs="+",
    )
    parser.add_argument(
        "--add_rv", help="Include a radial velocity shift.", action="store_true"
    )
    parser.add_argument(
        "--air", help="Convert to air wavelengths.", action="store_true"
    )
    parser.add_argument(
        "--correct",
        help="Apply Artigau et al. (2018) RV band corrections.",
        action="store_true",
    )

    parser.add_argument("-V", "--verbose", help="Enable verbose.", action="store_true")
    parser.add_argument(
        "--disable_normalization",
        help="Disable the convolution normalization.",
        action="store_true",
    )
    return parser.parse_args()


def do_analysis(
    star_params,
    vsini: float,
    R: float,
    band: str,
    sampling: float = 3.0,
    conv_kwargs=None,
    snr: float = 100.0,
    ref_band: str = "J",
    rv: float = 0.0,
    air: bool = False,
    model: str = "aces",
    verbose: bool = False,
) -> Tuple[Quantity, ...]:
    """Calculate RV precision and Quality for specific parameter set.

    Parameters
    ----------
    star_param:
        Stellar parameters [temp, logg, feh, alpha] for phoenix model libraries.
    vsini: float
       Stellar equatorial rotation.
    R: float
        Instrumental resolution.
    band: str
        Spectral band.
    sampling: float (default=False)
        Per pixel sampling (after convolutions)
    conv_kwargs: Dict (default=None)
        Arguments specific for the convolutions,
        'epsilon', 'fwhm_lim', 'num_procs', 'normalize', 'verbose'.
    snr: float (default=100)
        SNR normalization level. SNR per pixel and the center of the ref_band.
    ref_band: str (default="J")
        Reference band for SNR normalization.
    rv: float
        Radial velocity in km/s (default = 0.0).
    air: bool
        Get model in air wavelengths (default=False).
    model: str
        Name of synthetic library (aces, btsettl) to use. Default = 'aces'.
    verbose:
        Enable verbose (default=False).

    Returns
    -------
    q: astropy.Quality
     Spectral quality.
    result_1: astropy.Quality
        RV precision under condition 1.
    result_2 : astropy.Quality
        RV precision under condition 2.
    result_3: astropy.Quality
        RV precision under condition 3.

    Notes
    -----
        We apply the radial velocity doppler shift after
            - convolution (rotation and resolution)
            - resampling
            - SNR normalization.

        in this way the RV only effects the precision due to the telluric mask interaction.
        Physically the RV should be applied between the rotational and instrumental convolution
        but we assume this effect is negligible.

    """
    if conv_kwargs is None:
        conv_kwargs = {
            "epsilon": 0.6,
            "fwhm_lim": 5.0,
            "num_procs": num_cpu_minus_1,
            "normalize": True,
            "verbose": verbose,
        }

    if ref_band.upper() == "SELF":
        ref_band = band

    model = check_model(model)

    if model == "aces":
        wav, flux = load_aces_spectrum(star_params, photons=True, air=air)
    elif model == "btsettl":
        wav, flux = load_btsettl_spectrum(star_params, photons=True, air=air)
    else:
        raise Exception("Invalid model name reached.")

    wav_grid, sampled_flux = convolve_and_resample(
        wav, flux, vsini, R, band, sampling, **conv_kwargs
    )

    # Doppler shift
    try:
        if rv != 0:
            sampled_flux = doppler_shift_flux(wav_grid, sampled_flux, vel=rv)
    except Exception as e:
        print("Doppler shift was unsuccessful")
        raise e

    # Scale normalization for precision
    wav_ref, sampled_ref = convolve_and_resample(
        wav, flux, vsini, R, ref_band, sampling, **conv_kwargs
    )
    snr_normalize = snr_constant_band(
        wav_ref, sampled_ref, snr=snr, band=ref_band, sampling=sampling, verbose=verbose
    )
    sampled_flux = sampled_flux / snr_normalize

    if (ref_band == band) and verbose:
        mid_point = band_middle(ref_band)
        index_ref = np.searchsorted(
            wav_grid, mid_point
        )  # searching for the index closer to 1.25 micron
        snr_estimate = np.sqrt(np.sum(sampled_flux[index_ref - 1 : index_ref + 2]))
        print(
            "\tSanity Check: The S/N at {0:4.02} micron = {1:4.2f}, (should be {2:g}).".format(
                mid_point, snr_estimate, snr
            )
        )

    # Load Atmosphere for this band.
    atm = Atmosphere.from_band(band=band, bary=True).at(wav_grid)
    assert np.allclose(atm.wl, wav_grid), "The atmosphere does not cover the wav_grid"

    # Spectral Quality/Precision
    q = quality(wav_grid, sampled_flux)

    # Precision given by the first condition:
    result_1 = rv_precision(wav_grid, sampled_flux, mask=None)

    # Precision as given by the second condition
    result_2 = rv_precision(wav_grid, sampled_flux, mask=atm.mask)

    # Precision as given by the third condition: M = T**2
    result_3 = rv_precision(wav_grid, sampled_flux, mask=atm.transmission ** 2)

    # Turn quality back into a Quantity (to give it a .value method)
    q = q * u.dimensionless_unscaled
    return q, result_1, result_2, result_3


def convolve_and_resample(
    wav: ndarray,
    flux: ndarray,
    vsini: float,
    R: float,
    band: str,
    sampling: float,
    **conv_kwargs,
) -> Tuple[ndarray, ndarray]:
    """Convolve and resample functions together.

    Returns
    -------
    wav_grid: ndarray
        Resampled wavelength array
    sampled_flux: ndarray
        Convolved and resampled flux array

    """
    wav_band, __, convolved_flux = convolution(wav, flux, vsini, R, band, **conv_kwargs)
    # Re-sample to sampling per resolution element.
    wav_grid = log_resample(wav_band, sampling, R)
    sampled_flux = np.interp(wav_grid, wav_band, convolved_flux)
    return wav_grid, sampled_flux


def model_format_args(model, pars):
    """Format the model and parameter args to save in output file.
    Change to int, str and float.

    model in [temp, logg, fe/h, alpha]
    pars in order (R, band, vsini, sample).

    Can now also optionally handle a 5th parameter RV.

    """
    temp = int(model[0])
    logg = float(model[1])
    fe_h = float(model[2])
    alpha = float(model[3])
    band = pars[1]
    res = int(pars[0] / 1000)
    vsini = float(pars[2])
    sample = float(pars[3])
    try:
        rv = float(pars[4])
        return temp, logg, fe_h, alpha, band, res, vsini, sample, rv
    except IndexError:
        return temp, logg, fe_h, alpha, band, res, vsini, sample, 0.0


def header_row(add_rv=False) -> str:
    """Header row for output file."""
    if add_rv:
        header = (
            "temp,logg,feh,alpha,band,resolution,vsini,sampling,"
            "rv,correctflag,quality,cond1,cond2,cond3\n"
        )
    else:
        header = (
            "temp,logg,feh,alpha,band,resolution,vsini,sampling,"
            "correctflag,quality,cond1,cond2,cond3\n"
        )
    return header


def get_already_computed(filename: str, add_rv: bool = False) -> List[str]:
    """Get the string of already computed model/parameters from the result file."""
    param_lines = []
    with open(filename, "r") as f:
        for line in f:
            if add_rv:
                # First 10 columns
                param_lines.append(select_csv_columns(line, ncols=10))
            else:
                # First 9 columns
                param_lines.append(select_csv_columns(line, ncols=9))
        # Strip any spaces
        param_lines = [strip_whitespace(value) for value in param_lines]
    return param_lines


def is_already_computed(
    computed_values: List[str],
    model,
    pars,
    add_rv: bool = False,
    correct: bool = False,
    verbose=False,
) -> bool:
    """Check if any combinations have already been preformed.
    Correct is boolean for applied Artigau correction."""
    model_par_str_args = model_format_args(model, pars)
    rv_template = "{0:5d},{1:3.01f},{2:4.01f},{3:3.01f},{4:s},{5:3d}k,{6:4.01f},{7:3.01f},{8:3.01f},{9:1d}"
    no_rv_template = (
        "{0:5d},{1:3.01f},{2:4.01f},{3:3.01f},{4:s},{5:3d}k,{6:4.01f},{7:3.01f},{8:1d}"
    )
    if add_rv:
        if len(model_par_str_args) != 9:
            raise ValueError("model_par_str_args is incorrect length")

        idenifying_line = strip_whitespace(
            rv_template.format(*model_par_str_args, int(correct))
        )
    else:
        model_par_str_args = model_par_str_args[:8]
        if len(model_par_str_args) != 8:
            raise ValueError("model_par_str_args is incorrect length")

        idenifying_line = strip_whitespace(
            no_rv_template.format(*model_par_str_args, int(correct))
        )

    result = idenifying_line in computed_values

    if result and verbose:
        print(model_par_str_args, "already computed")
    return result


def strip_whitespace(line: str) -> str:
    return "".join(line.split())


def select_csv_columns(line: str, ncols: int = 8) -> str:
    """Select first ncols in a line from a csv.

    Parameters
    ----------
    ncols: int
        Number of column to select.

    Returns
    -------
    selected_cols: str
        Selected ncols of csv.

    """
    return ",".join(line.split(",")[:ncols])


def check_model(model: str) -> str:
    """Check model is 'aces' or 'btsettl'.

    Parameters
    ----------
    model: str
        Model name. Should be either 'aces' or 'btsettl'.

    Returns
    -------
    model: str
       Valid model output; either 'aces' or 'btsettl'.

"""
    if model == "phoenix":
        warnings.warn(
            "The model name 'phoenix' is depreciated, use 'aces' instead.",
            DeprecationWarning,
        )
        model = "aces"
    if model not in ["aces", "btsettl"]:
        raise ValueError(
            "Model name error of '{}'. Valid choices are 'aces' and 'btsettl'".format(
                model
            )
        )
    return model


if __name__ == "__main__":
    args = _parser()

    # check bt-settl parameters
    if args.model == "btsettl":
        if (args.metal != [0]) or (args.alpha != [0]):
            raise ValueError(
                "You cannot vary metallicity and alpha for BT-Settl, remove these flags."
            )
    try:
        normalize = not args.disable_normalization
    except AttributeError:
        normalize = True

    try:
        num_procs = args.num_procs
    except AttributeError:
        num_procs = num_cpu_minus_1

    conv_kwargs = {
        "epsilon": 0.6,
        "fwhm_lim": 5.0,
        "num_procs": num_procs,
        "normalize": normalize,
        "verbose": args.verbose,
    }

    snr = args.snr
    air = args.air
    if "ALL" in args.bands:
        args.bands.extend(config.bands["all"])
        args.bands = set(args.bands)  # Unique
    ref_band = args.ref_band

    # Load the relevant spectra
    if args.model == "aces":
        models_list = itertools.product(args.temp, args.logg, args.metal, args.alpha)
    else:
        models_list = itertools.product(args.temp, args.logg, [0], [0])

    if (args.rv != [0.0]) and (not args.add_rv):
        raise ValueError("Need to include --add_rv flag if a non-zero RV is used.")

    if not os.path.exists(args.output):
        with open(args.output, "a") as f:
            header = header_row(add_rv=args.add_rv)
            f.write(header)

    # Find all model/parameter combinations already computed.
    # To later skip recalculation.
    computed_values = get_already_computed(args.output, args.add_rv)

    # Iteration over stellar model (temp. logg, feh, alpha) then an inner loop over parameters (R, vsini, sampling, rv)
    for model in models_list:
        # Create generator for params_list
        params_list = itertools.product(
            args.resolution, args.bands, args.vsini, args.sampling, args.rv
        )
        for (R, band, vsini, sample, rv) in params_list:
            pars = (R, band, vsini, sample, rv)

            if is_already_computed(
                computed_values,
                model,
                pars,
                add_rv=args.add_rv,
                correct=args.correct,
                verbose=args.verbose,
            ):
                # skipping the recalculation
                continue
            else:
                precision_result: Tuple[Any, ...] = do_analysis(
                    model,
                    vsini=vsini,
                    R=R,
                    band=band,
                    conv_kwargs=conv_kwargs,
                    snr=snr,
                    ref_band=ref_band,
                    sampling=sample,
                    air=air,
                    model=args.model,
                )
                result: List[Optional[float]] = [
                    round(res.value, 1) if res is not None else None
                    for res in precision_result
                ]

                if args.correct:
                    # Apply Artigau 2018 Corrections
                    corr_value = correct_artigau_2018(band)
                    for ii, res in enumerate(result):
                        if (ii > 0) and (result[ii] is not None):  # Not the quality
                            result[ii] = res * corr_value

                result[0] = int(result[0]) if result[0] is not None else None

                if args.add_rv:
                    output_template = (
                        "{0:5d},{1:3.01f},{2:4.01f},{3:3.01f},{4:s},{5:3d}k,"
                        "{6:4.01f},{7:3.01f},{8:3.01f},{9:1d},{10:6d},{11:5.01f},"
                        "{12:5.01f},{13:5.01f}\n"
                    )
                    output_model_args = model_format_args(model, pars)
                else:
                    output_template = (
                        "{0:5d},{1:3.01f},{2:4.01f},{3:3.01f},{4:s},{5:3d}k,"
                        "{6:4.01f},{7:3.01f},{8:1d},{9:6d},{10:5.01f},{11:5.01f},"
                        "{12:5.01f}\n"
                    )
                    output_model_args = model_format_args(model, pars)[:8]

                linetowite = output_template.format(
                    *output_model_args,
                    int(args.correct),
                    result[0],
                    result[1],
                    result[2],
                    result[3],
                )
                with open(args.output, "a") as f:
                    f.write(strip_whitespace(linetowite) + "\n")  # Make csv only

    if args.verbose:
        print(
            "`{1}` completed successfully:\n"
            "\tYou shall find you results in '{0}'.".format(
                args.output, os.path.basename(__file__)
            )
        )

import argparse
import itertools
import os
from os.path import join
from typing import List, Tuple

import multiprocess as mprocess
import numpy as np
from astropy import units as u
from numpy import ndarray

import eniric
from eniric.atmosphere import Atmosphere
from eniric.broaden import convolution
from eniric.corrections import correct_artigau_2018
from eniric.Qcalculator import quality, rv_precision
from eniric.resample import log_resample
from eniric.snr_normalization import snr_constant_band
from eniric.utilities import (
    band_middle,
    doppler_shift_flux,
    load_aces_spectrum,
    load_btsettl_spectrum,
)

num_procs_minus_1 = mprocess.cpu_count() - 1

ref_choices = ["SELF"]
ref_choices.extend(eniric.bands["all"])


def _parser():
    """Take care of all the argparse stuff.

    :returns: the args
    """
    parser = argparse.ArgumentParser(
        description="Calculate quality for any library spectra."
    )
    parser.add_argument("-t", "--temp", type=int, help="Temperature", nargs="+")
    parser.add_argument(
        "-l",
        "--logg",
        type=float,
        default=[4.5],
        help="Logg, default = [4.5]",
        nargs="+",
    )
    parser.add_argument(
        "-m",
        "--metal",
        type=float,
        default=[0.0],
        help="Metallicity, default=[0.0]",
        nargs="+",
    )
    parser.add_argument(
        "-a",
        "--alpha",
        type=float,
        default=[0.0],
        help="Alpha, default=[0.0]",
        nargs="+",
    )
    parser.add_argument(
        "-s", "--sampling", type=float, default=[3.0], help="Sampling", nargs="+"
    )
    parser.add_argument(
        "-r",
        "--resolution",
        type=float,
        default=[50000],
        help="Instrumental resolution",
        nargs="+",
    )
    parser.add_argument(
        "-v",
        "--vsini",
        type=float,
        default=[1.0],
        help="Rotational Velocity",
        nargs="+",
    )
    parser.add_argument(
        "-b",
        "--bands",
        type=str,
        default="J",
        choices=eniric.bands["all"],
        help="Wavelength bands to select. Default=J.",
        nargs="+",
    )
    parser.add_argument(
        "--model",
        type=str,
        default="phoenix",
        choices=["phoenix", "btsettl"],
        help="Spectral models to use. Default=phoenix.",
    )
    parser.add_argument(
        "--save", default=False, action="store_true", help="Save results to file."
    )
    parser.add_argument(
        "--snr", help="Mid-band SNR scaling. (Default=100)", default=100, type=float
    )
    parser.add_argument(
        "--ref_band",
        help="SNR reference band. Default=J. (Default=100). "
        "'self' scales each band relative to the SNR itself.",
        choices=ref_choices,
        default="J",
        type=str,
    )
    parser.add_argument(
        "--num_procs",
        help="Number of processors to use. Default = number of cpus -1",
        default=num_procs_minus_1,
        type=int,
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Filename for results",
        default="quality_results.csv",
        type=str,
    )
    parser.add_argument(
        "--rv", help="Radial velocity value.", default=[0.0], type=float, nargs="+"
    )
    parser.add_argument(
        "--add_rv", help="Include radial velocity shift.", action="store_true"
    )
    parser.add_argument(
        "--air", help="Convert wavelengths from vacuum to air", action="store_true"
    )
    parser.add_argument("--correct", help="Apply RV corrections", action="store_true")
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
    model="phoenix",
):
    """Precision and Quality for specific parameter set.

    Parameters
    ----------
    air: bool
        Get model in air wavelengths.
    model: str
        Name of synthetic library to use. (phoenix, btsettl).
    rv: float
        Radial velocity.

    Notes:
        We apply the radial velocity doppler shift after
           - convolution (rotation and resolution)
           - resampling
           - SNR normalization.
        in this way the RV only effects the precision due to the telluric mask interaction.
        The RV should maybe come between the rotational and instrumental convolution
        but we assume this effect is negligible.
    """
    if conv_kwargs is None:
        conv_kwargs = {
            "epsilon": 0.6,
            "fwhm_lim": 5.0,
            "num_procs": 1,
            "normalize": True,
        }

    if ref_band.upper() == "SELF":
        ref_band = band

    if model == "phoenix":
        # Full photon count spectrum
        wav, flux = load_aces_spectrum(star_params, photons=True, air=air)
    elif model == "btsettl":
        wav, flux = load_btsettl_spectrum(star_params, photons=True, air=air)
    else:
        raise ValueError(
            "Model name error in '{}'. Valid choices are 'phoenix and 'btsettl'".format(
                model
            )
        )

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

    # Spectral Quality
    q = quality(wav_grid, sampled_flux)

    # Scale normalization for precision
    wav_ref, sampled_ref = convolve_and_resample(
        wav, flux, vsini, R, ref_band, sampling, **conv_kwargs
    )
    snr_normalize = snr_constant_band(
        wav_ref, sampled_ref, snr=snr, band=ref_band, sampling=sampling
    )
    sampled_flux = sampled_flux / snr_normalize

    if ref_band == band:
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
    prec1 = rv_precision(wav_grid, sampled_flux, mask=None)

    # Precision as given by the second condition
    prec2 = rv_precision(wav_grid, sampled_flux, mask=atm.mask)

    # Precision as given by the third condition: M = T**2
    prec3 = rv_precision(wav_grid, sampled_flux, mask=atm.transmission ** 2)

    # Turn quality back into a Quantity (to give it a .value method)
    q = q * u.dimensionless_unscaled
    return [q, prec1, prec2, prec3]


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
    wav_band, flux_band, convolved_flux = convolution(
        wav, flux, vsini, R, band, **conv_kwargs
    )
    # Re-sample to sampling per resolution element.
    wav_grid = log_resample(wav_band, sampling, R)
    sampled_flux = np.interp(wav_grid, wav_band, convolved_flux)
    return wav_grid, sampled_flux


def model_format_args(model, pars):
    """Format the model and parameter args to save in output file.
    Change to int, str and float.

    model in [temp, logg, fe/h, alpha]
    pars in order (R, band, vsini, sample).

    Can now also optionally handle a 5th parameter rv.
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
        return temp, logg, fe_h, alpha, band, res, vsini, sample


def header_row(add_rv=False):
    """Header row for output file."""
    if add_rv:
        header = (
            "Temp, logg, [Fe/H], Alpha, Band, Resolution, vsini, Sampling, "
            "RV, Quality, Cond. 1, Cond. 2, Cond. 3, correct flag\n"
        )
    else:
        header = (
            "Temp, logg, [Fe/H], Alpha, Band, Resolution, vsini, Sampling, "
            "Quality, Cond. 1, Cond. 2, Cond. 3, correct flag\n"
        )
    return header


def get_already_computed(filename: str, add_rv: bool = False):
    """Get the string of already computer model/parameters from the result file."""
    computed_values = []
    with open(filename, "r") as f:
        for line in f:
            if add_rv:
                # should be longer
                computed_values.append(line[:47])
                computed_values.append(line[:40])
            else:
                computed_values.append(line[:41])
                computed_values.append(line[:34])
    return computed_values


def is_already_computed(computed_values: List[str], model, pars, add_rv: bool = False):
    """Check if any combinations have already been preformed"""
    model_par_str_args = model_format_args(model, pars)

    if add_rv:
        if len(model_par_str_args) != 9:
            raise ValueError("model_par_str_args is incorrect length")

        params_one = (
            "{0:5d}, {1:3.01f}, {2:4.01f}, {3:3.01f}, {4:s}, {5:3d}k,"
            " {6:4.01f}, {7:3.01f}"
        ).format(*model_par_str_args)
        # may change output to have less spaces in future
        params_two = (
            "{0:5d},{1:3.01f},{2:4.01f},{3:3.01f},{4:s},{5:3d}k," "{6:4.01f},{7:3.01f}"
        ).format(*model_par_str_args)

    else:

        model_par_str_args = model_par_str_args[:8]

        if len(model_par_str_args) != 8:
            raise ValueError("model_par_str_args is incorrect length")

        params_one = (
            "{0:5d}, {1:3.01f}, {2:4.01f}, {3:3.01f}, {4:s}, {5:3d}k,"
            " {6:4.01f}, {7:3.01f}"
        ).format(*model_par_str_args)
        # may change output to have less spaces in future
        params_two = (
            "{0:5d},{1:3.01f},{2:4.01f},{3:3.01f},{4:s},{5:3d}k," "{6:4.01f},{7:3.01f}"
        ).format(*model_par_str_args)

    result = (params_one in computed_values) or (params_two in computed_values)
    if result:
        print(model_par_str_args, "model already computed")
        print(params_one)

    return result


if __name__ == "__main__":
    args = _parser()

    # check bt-settl parameters
    if args.model == "btsettl":
        if (args.metal != [0]) or (args.alpha != [0]):
            raise ValueError(
                "You cannot vary metallicity and alpha for BT-Settl, remove these flags."
            )

    try:
        num_procs = args.num_procs
    except AttributeError:
        num_procs = num_procs_minus_1

    try:
        normalize = args.normalzie
    except AttributeError:
        normalize = True

    snr = args.snr
    air = args.air
    if "ALL" in args.bands:
        args.bands.extend(eniric.bands["all"])
        args.bands = set(args.bands)  # Unique
    ref_band = args.ref_band

    conv_kwargs = {
        "epsilon": 0.6,
        "fwhm_lim": 5.0,
        "num_procs": num_procs,
        "normalize": normalize,
    }

    # Load the relevant spectra
    if args.model == "phoenix":
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

    with open(args.output, "a") as f:

        for model in models_list:
            # Create generator for params_list
            params_list = itertools.product(
                args.resolution, args.bands, args.vsini, args.sampling, args.rv
            )

            for (R, band, vsini, sample, rv) in params_list:
                pars = (R, band, vsini, sample, rv)

                if is_already_computed(
                    computed_values, model, pars, add_rv=args.add_rv
                ):
                    # skipping the recalculation
                    continue
                else:
                    result = do_analysis(
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
                    result = [
                        round(res.value, 1) if res is not None else None
                        for res in result
                    ]

                    if args.correct:
                        # Apply Artigau 2018 Corrections
                        corr_value = correct_artigau_2018(band)
                        for ii, res in enumerate(result):
                            if ii > 0:  # Not the quality
                                result[ii] = res * corr_value

                    result[0] = int(result[0]) if result[0] is not None else None

                    if args.add_rv:
                        output_template = (
                            "{0:5d}, {1:3.01f}, {2:4.01f}, {3:3.01f}, {4:s}, {5:3d}k,"
                            " {6:4.01f}, {7:3.01f}, {8:3.01f}, {9:6d}, {10:5.01f}, "
                            "{11:5.01f}, {12:5.01f}, {13:1d}\n"
                        )
                        output_model_args = model_format_args(model, pars)
                    else:
                        output_template = (
                            "{0:5d}, {1:3.01f}, {2:4.01f}, {3:3.01f}, {4:s}, {5:3d}k,"
                            " {6:4.01f}, {7:3.01f}, {8:6d}, {9:5.01f}, {10:5.01f}, "
                            "{11:5.01f}, {12:1d}\n"
                        )
                        output_model_args = model_format_args(model, pars)[:8]
                    f.write(
                        output_template.format(
                            *output_model_args,
                            result[0],
                            result[1],
                            result[2],
                            result[3],
                            int(args.correct),
                        )
                    )
    print("Done")

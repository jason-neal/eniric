#!/usr/bin/env python
# Script to perform a convolution on a spectrum.
# Can take a number of parameters if needed
import argparse
import os
import sys
from datetime import datetime as dt
from typing import List, Optional, Union, Sequence

import eniric
from eniric.nIRanalysis import convolve_spectra
from eniric.resample import resampler
from eniric.utilities import get_spectrum_name, resolutions2ints


def _parser():
    """Take care of all the argparse stuff.

    :returns: the args
    """
    parser = argparse.ArgumentParser(
        description="Perform Convolution in preparation for nIR_precision."
    )
    parser.add_argument(
        "-s", "--startype", help='Spectral Type e.g "MO"', type=str, nargs="+"
    )
    parser.add_argument(
        "-v", "--vsini", help="Rotational velocity of source", type=float, nargs="+"
    )
    parser.add_argument(
        "-R",
        "--resolution",
        help="Observational resolution e.g. 100000 or 100k",
        type=str,
        nargs="+",
    )
    parser.add_argument(
        "-b",
        "--band",
        type=str,
        default="ALL",
        choices=["ALL", "VIS", "GAP", "Z", "Y", "J", "H", "K"],
        help="Wavelength band to select",
        nargs="+",
    )
    parser.add_argument(
        "--sample_rate",
        default=[3.0],
        type=float,
        nargs="+",
        help="Resample rate, pixels per FWHM. Default=3.0",
    )
    parser.add_argument(
        "--noresample", help="Don't Resample output", default=False, action="store_true"
    )
    parser.add_argument(
        "--unnormalized", help="Normalize for wavelength step", action="store_true"
    )
    parser.add_argument(
        "--org",
        help="Only use original .dat files, (temporary option)",
        default=False,
        action="store_true",
    )
    parser.add_argument(
        "-r",
        "--replace",
        action="store_true",
        help="Replace data files if already created.",
    )
    return parser.parse_args()


def main(
        startype: Sequence[str],
        vsini: Sequence[float],
        resolution: Sequence[str],
        band: Sequence[str],
        sample_rate: Optional[List[Union[int, float]]] = None,
        noresample: bool = False,
        unnormalized: bool = False,
        org: bool = False,
        replace: bool = False,
):
    """Run convolutions of NIR spectra for the range of given parameters.

    Multiple values of startype, vsini, resolution, band, and sample_rate can
    be provided.

    Read files from eniric.paths["phoenix_dat"]"

    Parameters
    ----------
    startype: list of strings
    vsini: list of floats
    resolution: list of floats
    band: list of strings
    sample_rate: list of floats default=[3.0]
    noresample: bool default=False
    unnormalized: bool default=False

    """
    normalize = not unnormalized
    if sample_rate is None:
        # Default sample rate
        sample_rate = [3.0]

    # Check the inputs are correct format. (lists)
    for f_input, f_name in zip(
            [startype, band, vsini, resolution, sample_rate],
            ["startype", "band", "vsini", "resolution", "sample_rate"],
    ):
        if not isinstance(f_input, list):
            print(f_name, type(f_input), type(f_name))
            raise TypeError("Input {0} is not list".format(f_name))

    # vsini, resolution, band and sample_rate can all be a series of values
    vsini = [float(v) for v in vsini]  # turn to floats

    # Handle K in Resolution
    res = resolutions2ints(resolution)

    start_time = dt.now()

    phoenix_path = eniric.paths["phoenix_dat"]

    results_dir = eniric.paths["results"]
    os.makedirs(results_dir, exist_ok=True)

    resampled_dir = eniric.paths["resampled"]
    os.makedirs(resampled_dir, exist_ok=True)

    if not normalize:
        norm_ = "_unnormalized"
    else:
        norm_ = ""

    counter = 0
    for star in startype:
        spectrum_name = os.path.join(phoenix_path, get_spectrum_name(star, org=org))

        for b in band:
            for vel in vsini:
                for R in res:
                    for sample in sample_rate:
                        r = int(float(R) / 1000)
                        result_name = "Spectrum_{0}-PHOENIX-ACES_{1}band_vsini{2:.01f}_R{3}k{4}.txt".format(
                            star, b, vel, r, norm_
                        )

                        if (
                                os.path.exists(os.path.join(results_dir, result_name))
                                and not replace
                        ):
                            print(
                                "Skipping convolution as {} already exists".format(
                                    result_name
                                )
                            )
                            continue

                        print("Name to be the result file", result_name)
                        convolve_spectra(
                            spectrum_name,
                            b,
                            vel,
                            R,
                            epsilon=0.6,
                            plot=False,
                            fwhm_lim=5.0,
                            num_procs=None,
                            results_dir=results_dir,
                            normalize=normalize,
                            output_name=result_name,
                        )

                        # Resample only the file just made
                        if not noresample:
                            resampler(
                                result_name,
                                results_dir=results_dir,
                                resampled_dir=resampled_dir,
                                sampling=sample,
                            )
                        counter += 1

    print(
        "Time to convolve {0:d} combinations = {1}".format(
            counter, dt.now() - start_time
        )
    )
    return 0


if __name__ == "__main__":
    args = vars(_parser())
    opts = {k: args[k] for k in args}

    sys.exit(main(**opts))

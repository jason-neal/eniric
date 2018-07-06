#!/usr/bin/env python
"""Script to check if all the results files are ready for precision calculation

 e.g. that they exist before trying to calculate the precision.

 """
import sys
import argparse
import itertools

import eniric


def _parser():
    """Take care of all the argparse stuff.

    :returns: the args
    """
    parser = argparse.ArgumentParser(
        description="Calculate perfect precision of all convolved spectra."
    )

    parser.add_argument(
        "-s",
        "--startype",
        help='Spectral Type e.g "MO"',
        type=str,
        nargs="*",
        default=None,
    )
    parser.add_argument(
        "-v",
        "--vsini",
        help="Rotational velocity of source",
        type=float,
        nargs="*",
        default=None,
    )
    parser.add_argument(
        "-R",
        "--resolution",
        help="Observational resolution",
        type=float,
        nargs="*",
        default=None,
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
        default=None,
        type=float,
        nargs="*",
        help="Resample rate, pixels per FWHM. Default=3.0",
    )
    parser.add_argument(
        "--normalize",
        help="Use convolution normalized spectra",
        default=True,
        action="store_false",
    )
    return parser.parse_args()


def main(
    startype=None,
    vsini=None,
    resolution=None,
    band=None,
    sample_rate=None,
    normalize=True,
):
    """Check if all the results files are ready for precision calculation."""
    if startype is None:
        spectral_types = ["M0", "M3", "M6", "M9"]
    else:
        spectral_types = startype
    if band is None:
        bands = ["Z", "Y", "J", "H", "K"]
    else:
        bands = band
    if vsini is None:
        vsini = ["1.0", "5.0", "10.0"]

    if resolution is None:
        resolution = ["60k", "80k", "100k"]
    else:
        resolution = ["{0:.0f}k".format(R / 1000) for R in resolution]

    if sample_rate is None:
        sampling = ["3"]
    else:
        sampling = sample_rate

    if not normalize:
        norm_ = "_unnormalized"
    else:
        norm_ = ""

    iterations = itertools.product(spectral_types, bands, vsini, resolution, sampling)
    for (star, band, vel, res, smpl) in iterations:

        file_to_read = (
            "Spectrum_{0:s}-PHOENIX-ACES_{1:s}band_vsini{2:.01f}_R{3:s}{5:s}_res{4:3.01f}.txt"
            ""
        ).format(star, band, float(vel), res, float(smpl), norm_)

        full_name = os.path.join(eniric.paths["resampled"], file_to_read)

        # Find if the file exists.
        if not os.path.exists(full_name):
            print("{} does not exist!".format(filename))

    return 0


if __name__ == "__main__":
    args = vars(_parser())
    opts = {k: args[k] for k in args}
    sys.exit(main(**opts))

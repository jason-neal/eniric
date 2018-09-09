#!/usr/bin/env python
"""Pre-doppler-shift the Tapas atmosphere model for RV Precision.

To make RV precision calculations faster.
"""

import argparse
import sys
from os.path import join
from typing import List, Optional

import numpy as np

import eniric
from eniric.atmosphere import Atmosphere

choices = [None, "ALL"]
choices.extend(eniric.bands["all"])


def _parser():
    """Take care of all the argparse stuff.

    :returns: the args
    """
    parser = argparse.ArgumentParser(
        description=("Pre-doppler shift atmosphere masks.")
    )

    parser.add_argument(
        "-b",
        "--bands",
        type=str,
        default="ALL",
        choices=choices,
        help="Wavelength bands to select. Default=None.",
        nargs="+",
    )
    parser.add_argument(
        "-p",
        "--plot",
        action="store_true",
        default=False,
        help="Plot the atmosphere model with masks.",
    )
    _args = parser.parse_args()
    return _args


def main(bands: Optional[List[str]] = None, plot: bool = False):
    """Preform the barycentric shifting of atmosphere masks and saves result.

    This saves time in the precision determination code.

    Parameters
    ----------
    bands: list of str
        Wavelength bands to perform barycenter shifts on. Default is all bands.
    plot: bool
        Flag to plot test plots of masks.
    """
    if (bands is None) or ("ALL" in bands):
        bands_ = eniric.bands["all"]
    else:
        bands_ = bands

    for band in bands_:
        unshifted_atmmodel = join(
            eniric.paths["atmmodel"],
            "{0}_{1}.txt".format(eniric.atmmodel["base"], band),
        )

        print("Reading atmospheric model...", unshifted_atmmodel)

        atm = Atmosphere.from_file(unshifted_atmmodel)

        print("Calculating impact of Barycentric movement on mask...")
        org_mask = atm.mask
        masked_before = np.sum(org_mask)
        atm.bary_shift_mask(consecutive_test=True)

        masked_after = np.sum(atm.mask)
        print(
            "Masked fraction before = {0:0.03f}".format(
                (len(org_mask) - masked_before) / len(org_mask)
            )
        )
        print(
            "Masked fraction after = {0:0.03f}".format(
                (len(atm.mask) - masked_after) / len(atm.mask)
            )
        )

        shifted_atmmodel = unshifted_atmmodel.replace(".txt", "_bary.txt")
        print("Saving doppler-shifted atmosphere model to {}".format(shifted_atmmodel))

        header = ["# atm_wav(nm)", "atm_flux", "atm_std_flux", "atm_mask"]
        atm.to_file(fname=shifted_atmmodel, header=header, fmt="%11.8f")
    print("Done")


if __name__ == "__main__":
    args = vars(_parser())
    opts = {k: args[k] for k in args}
    sys.exit(main(**opts))

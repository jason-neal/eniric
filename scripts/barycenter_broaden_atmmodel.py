#!/usr/bin/env python
"""
barycenter_broaden_atmmodel
---------------------------
Doppler shift the Tapas atmosphere model and save to files.
This makes the RV precision calculations faster.

"""

import argparse
import sys
from os.path import join
from typing import List, Optional

import numpy as np

from eniric import config
from eniric.atmosphere import Atmosphere

choices = [None, "ALL"]
choices.extend(config.bands["all"])


def _parser():
    """Take care of all the argparse stuff."""
    parser = argparse.ArgumentParser(
        description=("Pre-doppler shift atmosphere masks.")
    )

    parser.add_argument(
        "-b",
        "--bands",
        type=str,
        default="ALL",
        choices=choices,
        help="Wavelength bands to select. Default is None.",
        nargs="+",
    )
    parser.add_argument("-v", "--verbose", help="Turn on verbose.", action="store_true")
    return parser.parse_args()


def main(bands: Optional[List[str]] = None, verbose: bool = False) -> None:
    """Preform the barycentric shifting of atmosphere masks and saves result.

    This saves time in the precision determination code.

    Parameters
    ----------
    bands: list of str
        Wavelength bands to perform barycenter shifts on. Default is all bands.

    """
    if (bands is None) or ("ALL" in bands):
        bands_ = config.bands["all"]
    else:
        bands_ = bands

    for band in bands_:
        unshifted_atmmodel = join(
            config.pathdir,
            config.paths["atmmodel"],
            "{0}_{1}.dat".format(config.atmmodel["base"], band),
        )
        if verbose:
            print("Reading atmospheric model...", unshifted_atmmodel)

        atm = Atmosphere.from_file(unshifted_atmmodel)

        if verbose:
            print("Calculating impact of Barycentric movement on mask...")
        org_mask = atm.mask
        masked_before = np.sum(org_mask)
        atm.barycenter_broaden(consecutive_test=True)

        masked_after = np.sum(atm.mask)
        if verbose:
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

        shifted_atmmodel = unshifted_atmmodel.replace(".dat", "_bary.dat")
        if verbose:
            print(
                "Saving doppler-shifted atmosphere model to {}".format(shifted_atmmodel)
            )

        header = ["# atm_wav(nm)", "atm_flux", "atm_std_flux", "atm_mask"]
        atm.to_file(fname=shifted_atmmodel, header=header, fmt="%11.8f")

    if verbose:
        print("Finished barycentric shifting of atmosphere masks")


if __name__ == "__main__":
    args = vars(_parser())
    opts = {k: args[k] for k in args}
    try:
        main(**opts)
        sys.exit(0)
    except:
        raise

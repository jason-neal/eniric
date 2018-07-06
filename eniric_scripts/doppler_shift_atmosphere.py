#!/usr/bin/env python
"""Pre-doppler-shift the Tapas atmosphere model for RV Precision.

To make RV precision calculations faster.
"""

import argparse
import sys

import numpy as np

import eniric.atmosphere as atm
import eniric.IOmodule as io


def _parser():
    """Take care of all the argparse stuff.

    :returns: the args
    """
    parser = argparse.ArgumentParser(
        description=("Calculate radial velocity" "precision of model spectra.")
    )

    parser.add_argument(
        "-b",
        "--bands",
        type=str,
        default="J",
        choices=["ALL", "VIS", "GAP", "Z", "Y", "J", "H", "K", None, "visible"],
        help="Wavelength bands to select. Default=J.",
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


def main(bands=None, plot=False):
    """Preform the barycentric shifting of atmosphere masks and saves result.

    This saves time in the precision determination code.

    Parameters
    ----------
    bands: list of str
        Wavelength bands to perform shifts on
    plot: bool
        Flag to plot test plots of masks.
    """
    if bands is None:
        bands = ["Z", "Y", "J", "H", "K"]

    for band in bands:
        unshifted_atmmodel = "../data/atmmodel/Average_TAPAS_2014_{}.txt".format(band)
        shifted_atmmodel = unshifted_atmmodel.replace(".txt", "_bary.txt")

        print("Reading atmospheric model...", unshifted_atmmodel)

        wav_atm, flux_atm, std_flux_atm, mask_atm = atm.prepare_atmosphere(
            unshifted_atmmodel
        )
        print(
            ("There were {0:d} unmasked pixels out of {1:d}., or {2:.1%}." "").format(
                np.sum(mask_atm), len(mask_atm), np.sum(mask_atm) / len(mask_atm)
            )
        )
        print(
            "The model ranges from {0:4.2f} to {1:4.2f} micron.".format(
                wav_atm[0], wav_atm[-1]
            )
        )
        print("Done.")

        print("Calculating impact of Barycentric movement on mask...")
        org_mask = mask_atm
        mask_atm = atm.barycenter_shift(wav_atm, mask_atm)

        print("Saving doppler-shifted atmosphere model...", shifted_atmmodel)

        header = ["# atm_wav(nm)", "atm_flux", "atm_std_flux", "atm_mask"]

        # Turn wav_atm back to nanometers for saving.
        io.pdwrite_cols(
            shifted_atmmodel,
            wav_atm * 1000,
            flux_atm,
            std_flux_atm,
            mask_atm,
            header=header,
            float_format="%11.8f",
        )

        if plot:
            atm.plot_atm_masks(
                wav_atm, flux_atm, org_mask, new_mask=mask_atm, block=True
            )


if __name__ == "__main__":
    args = vars(_parser())
    opts = {k: args[k] for k in args}
    sys.exit(main(**opts))

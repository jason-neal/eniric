#!/usr/bin/env python
"""Print the point at which the scale occurs."""
import os

import matplotlib.pyplot as plt

import eniric
import eniric.IOmodule as io
import eniric.utilities as utils

path = "plots/band_centers/"
os.makedirs(path, exist_ok=True)

for band in ["Z", "H", "J", "Y", "K"]:
    band_limits = utils.band_limits(band)
    plt.figure()
    plt.vlines(sum(band_limits) / 2, 0, 1, linestyle="--")
    for sptype in ["M0", "M3", "M6", "M9"]:
        file_to_read = "Spectrum_{0}-PHOENIX-ACES_{1}band_vsini1.0_R100k_res3.txt".format(
            sptype, band
        )
        wave_stellar, flux_stellar = io.pdread_2col(
            os.path.join(eniric.paths["resampled"], file_to_read)
        )

        plt.plot(wave_stellar, flux_stellar / 1e10, label=sptype)

    plt.title("{0} Band SNR scaling point".format(band))
    plt.xlabel("Wavelength (micron)")
    plt.ylabel("flux")
    plt.legend(loc=1)
    plt.xlim([sum(band_limits) / 2 - 0.002, sum(band_limits) / 2 + 0.002])
    plt.savefig(
        os.path.join(path, "{}_band_center_marking.png".format(band)),
        bbox_inches="tight",
    )

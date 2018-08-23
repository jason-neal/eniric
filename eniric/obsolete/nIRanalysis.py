#!/usr/bin/env python
import os
import sys
from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc
from numpy import ndarray

import eniric
import eniric.io_module as io
from eniric.broaden import convolution
from eniric.obsolete.utilities import read_spectrum

# set stuff for latex usage
rc("text", usetex=True)

data_rep = eniric.paths["phoenix_dat"]
results_dir = eniric.paths["results"]
resampled_dir = eniric.paths["resampled"]

# models form PHOENIX-ACES
M0_ACES = data_rep + "lte03900-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes_wave.dat"
M3_ACES = data_rep + "lte03500-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes_wave.dat"
M6_ACES = data_rep + "lte02800-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes_wave.dat"
M9_ACES = data_rep + "lte02600-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes_wave.dat"


def run_convolutions(spectrum_string: str, band: str) -> None:
    """
    Runs the convolutions for a set of spectra in batch
    """
    vsini = [1.0, 5.0, 10.0]
    R = [60000, 80000, 100000]

    spectrum = spectrum_string  # removed exec usage
    print(spectrum)
    print(
        "Running the convolutions for spectra of {0:s} in band {1:s}\n.".format(
            spectrum, band
        )
    )
    for vel in vsini:
        for res in R:
            convolve_spectra(spectrum, band, vel, res)


def save_convolution_results(
    filename: str, wavelength: ndarray, flux: ndarray, convolved_flux: ndarray
) -> int:
    """Saves convolution results to a file.

    Parameters
    ----------
    filename: str
    wavelength: array-like
     flux, convolved_flux
    """
    print("Saving results...")

    # Note: difference in sampling at 1.0 and 1.5 microns makes jumps
    # in the beginning of Y and H bands
    eniric.obsolete.IOmodule.write_e_3col(filename, wavelength, flux, convolved_flux)
    print("Done.")
    return 0


def convolve_spectra(
    spectrum,
    band,
    vsini,
    R,
    epsilon: float = 0.6,
    fwhm_lim: float = 5.0,
    num_procs: Optional[int] = None,
    results_dir: str = results_dir,
    normalize: bool = True,
    output_name: Optional[str] = None,
) -> int:
    """Load Spectrum, apply convolution and then save results."""
    print("Reading the data...")
    wav, flux = read_spectrum(spectrum)  # In microns and  photon flux.
    print("Done.")

    wav_band, flux_band, convolved_flux = convolution(
        wav,
        flux,
        vsini,
        R,
        band,
        epsilon=epsilon,
        fwhm_lim=fwhm_lim,
        num_procs=num_procs,
        normalize=normalize,
    )
    if not normalize:
        norm_ = "_unnormalized"
    else:
        norm_ = ""

    if output_name is None:
        name_model = name_assignment(spectrum)

        filename = "{0}Spectrum_{1}_{2}band_vsini{3:3.1f}_R{4:d}k{5}.txt".format(
            results_dir, name_model, band, vsini, R / 1000, norm_
        )
    else:
        filename = os.path.join(results_dir, output_name)

    save_convolution_results(filename, wav_band, flux_band, convolved_flux)

    return 0


###############################################################################
def name_assignment(spectrum: str):
    """
    assigns a name to the filename in which the spectrum is going to be saved
    """
    # Simplified to temperature and base in spectrum name.
    m0_aces = "lte03900"
    m3_aces = "lte03500"
    m6_aces = "lte02800"
    m9_aces = "lte02600"
    base = "PHOENIX-ACES-AGSS-COND-2011-HiRes_wave.dat"
    if (m0_aces in spectrum) and (base in spectrum):
        name = "M0-PHOENIX-ACES"
    elif (m3_aces in spectrum) and (base in spectrum):
        name = "M3-PHOENIX-ACES"
    elif (m6_aces in spectrum) and (base in spectrum):
        name = "M6-PHOENIX-ACES"
    elif (m9_aces in spectrum) and (base in spectrum):
        name = "M9-PHOENIX-ACES"
    else:
        raise ValueError("Name {0} not found!".format(spectrum))
    return name


if __name__ == "__main__":
    if len(sys.argv) == 3:
        run_convolutions(sys.argv[1], sys.argv[2])
    else:
        print("Arguments not compatible with called function.")

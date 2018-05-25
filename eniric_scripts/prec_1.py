#!/usr/bin/env python
"""Quick and dirty precision 1.

Doesn't involve atmosphere model so can perform relatively easily to check
precision is working.
"""

import argparse
import os
import sys

import numpy as np

import eniric
import eniric.IOmodule as io
import eniric.Qcalculator as Qcalculator
from eniric.snr_normalization import normalize_flux


def _parser():
    """Take care of all the argparse stuff.

    :returns: the args
    """
    parser = argparse.ArgumentParser(description='Calculate perfect precision of all convolved spectra.')

    parser.add_argument('-s', '--startype', help='Spectral Type e.g "MO"', type=str, nargs="*", default=None)
    parser.add_argument("-v", "--vsini", help="Rotational velocity of source",
                        type=float, nargs="*", default=None)
    parser.add_argument("-R", "--resolution", help="Observational resolution",
                        type=float, nargs="*", default=None)
    parser.add_argument("-b", "--band", type=str, default=["ALL"],
                        choices=["ALL", "VIS", "GAP", "Z", "Y", "J", "H", "K"],
                        help="Wavelength band to select", nargs="+")
    parser.add_argument('--sample_rate', default=None, type=float, nargs="*",
                        help="Resample rate, pixels per FWHM. Default=3.0")
    parser.add_argument('--normalize', help='Turn off convolution normalization.', action="store_false")
    return parser.parse_args()


# atmmodel = "../data/atmmodel/Average_TAPAS_2014.txt"
resampled_dir = eniric.paths["resampled"]


def calc_prec1(star, band, vel, resolution, smpl, normalize=True):
    """Just calculate precision for 1st case.

    Resolution in short form e.g 100k

    Loads in the file, and calculates RV precision on full band.
    """
    vel = float(vel)

    if not normalize:
        norm_ = "_unnormalized"
        norm_id = "-unnorm"
    else:
        norm_ = ""
        norm_id = ""
    print(star, band, vel, resolution, smpl, norm_)
    print(type(star), type(band), type(vel), type(resolution), type(smpl), type(norm_))
    file_to_read = ("Spectrum_{0:s}-PHOENIX-ACES_{1:s}band_vsini{2:.01f}_R{3:s}{5:s}_res{4:d}.txt"
                    "").format(star, band, vel, resolution, smpl, norm_)

    # sample was left aside because only one value existed
    id_string = "{0}-{1}-{2:.01f}-{3}{4}".format(star, band, vel, resolution, norm_id)

    wav_stellar, flux_stellar = io.pdread_2col(os.path.join(eniric.paths["resampled"], file_to_read))

    # removing boundary effects
    wav_stellar = wav_stellar[2:-2]
    flux_stellar = flux_stellar[2:-2]

    # Normalize to SNR 100 in middle of J band 1.25 micron!
    flux_stellar = normalize_flux(flux_stellar, id_string)

    if id_string in ["M0-J-1.0-100k", "M3-J-1.0-100k", "M6-J-1.0-100k", "M9-J-1.0-100k"]:
        index_reference = np.searchsorted(wav_stellar, [1.25])[0]  # searching for the index closer to 1.25 micron
        sn_estimate = np.sqrt(np.sum(flux_stellar[index_reference - 1:index_reference + 2]))
        print("\tSanity Check: The S/N for the {0:s} reference model was of {1:4.2f}.".format(id_string, sn_estimate))

    elif "J" in id_string:
        index_reference = np.searchsorted(wav_stellar, [1.25])[0]  # searching for the index closer to 1.25 micron
        sn_estimate = np.sqrt(np.sum(flux_stellar[index_reference - 1:index_reference + 2]))
        print(
            "\tSanity Check: The S/N for the {0:s} non-reference model was of {1:4.2f}.".format(id_string, sn_estimate))

    # print("Performing analysis for: ", id_string)
    prec_1 = Qcalculator.RVprec_calc(wav_stellar, flux_stellar)

    # print("{0}: \t{1}".format(id_string, prec_1))
    return id_string, prec_1


def main(startype=None, vsini=None, resolution=None, bands=None, sample_rate=None, normalize: bool = True):
    """Script that calculates the RV precision without atmosphere.

    Parameters
    ----------
    startype: List[str]
        Spectral type of star.
    vsini: List [str, float]
        Rotation of star.
    resolution: List[str, int]
            Spectral resolutions.
    bands: list[str]
        Spectral bands to use.
    sample_rate: Optional[int]
        Sample rate of spectrum. Default = 3.0.
    normalize: bool
        Normalize the convolution. Default=True.

    """
    if startype is None:
        startype = ["M0", "M3", "M6", "M9"]

    if bands is None:
        bands = ["Z", "Y", "J", "H", "K"]

    if vsini is None:
        vsini = ["1.0", "5.0", "10.0"]

    if resolution is None:
        resolution = ["60k", "80k", "100k"]
    else:
        resolution = ["{0:.0f}k".format(R / 1000) for R in resolution]

    if sample_rate is None:
        sample_rate = ["3"]

    if normalize is None:
        normalize = [True, False]
    else:
        if isinstance(normalize, bool):
            normalize = [normalize]

    # Check the inputs are correct format. (lists)
    for f_input, f_name in zip([startype, bands, vsini, resolution, sample_rate],
                               ["startype", "band", "vsini", "resolution", "sample_rate"]):
        if not isinstance(f_input, list):
            print(f_name, type(f_input), type(f_name))
            raise TypeError("Input {0} is not list".format(f_name))

    precision = {}  # dict for storing precision value

    # TODO: iterate over band last so that the J band normalization value can be
    # Obtained first and applied to each band.

    for star in startype:
        for band in bands:
            for vel in vsini:
                for R in resolution:
                    for smpl in sample_rate:
                        for norm in normalize:
                            try:
                                id_string, prec_1 = calc_prec1(star, band, vel, R, smpl, normalize=norm)
                                precision[id_string] = prec_1
                            except FileNotFoundError:
                                print("File Not found ", star, band, vel, R, smpl, norm)
                                continue

    print("id_string\t\tprec_1")
    for key in precision:
        print("{0:s} \t\t{1:02.1f}".format(key, precision[key]))

    return 0


if __name__ == '__main__':
    print("hello")
    args = vars(_parser())
    opts = {k: args[k] for k in args}
    sys.exit(main())

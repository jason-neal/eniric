"""Script to check if all the results files are ready for precision calculation

 e.g. that they exist before trying to calculate the precision.

 Jason Neal - January 2017
 """

import argparse
import itertools


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
    parser.add_argument("-b", "--band", type=str, default="ALL",
                        choices=["ALL", "VIS", "GAP", "Z", "Y", "J", "H", "K"],
                        help="Wavelength band to select", nargs="+")
    parser.add_argument('-d', '--data_dir', help='Data directory', type=str, default=None)
    parser.add_argument('--sample_rate', default=None, type=float, nargs="*",
                        help="Resample rate, pixels per FWHM. Default=3.0")
    parser.add_argument('--results', default=None, type=str,
                        help='Result directory Default=data_dir+"/results/"')
    parser.add_argument('--resamples', default=None, type=str,
                        help='Resample directory. Default=data_dir+"/resampled/"')
    parser.add_argument('--noresample', help='Resample output', default=False,
                        action="store_true")
    parser.add_argument('--normalize', help='Use convolution normalized spectra', default=True,
                        action="store_false")
    parser.add_argument('--org', help='Only use original .dat files, (temporary option)',
                        default=False, action="store_true")
    return parser.parse_args()


def main(startype=None, vsini=None, resolution=None, band=None, data_dir=None,
         results=None, resamples=None, sample_rate=3.0, noresample=False,
         normalize=True, org=False):
    """Check if all the results files are ready for precision calculation."""
    if data_dir is None:
        data_dir = "../data/"

    if results is None:
        results_dir = data_dir + "results/"
    else:
        results_dir = results

    if resamples is None:
        resampled_dir = data_dir + "resampled/"
    else:
        resampled_dir = resamples

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

    iterations = itertools.product(spectral_types, bands, vsini, resolution, sampling)
    for (star, band, vel, res, smpl) in iterations:

        # Find if the file exists.
        if exists:
            pass
        else:
            print("{} does not exist!".format(filename))

    return 0


if __name__ == "__main__":
    args = vars(_parser())
    opts = {k: args[k] for k in args}
    sys.exit(main(**opts))

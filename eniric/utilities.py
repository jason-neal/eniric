
"""
Auxiliary functions for nIRanalysis

"""
import re
import os
import errno
import numpy as np
import matplotlib.pyplot as plt
from eniric.IOmodule import pdread_2col


def read_spectrum(spec_name):
    """ Function that reads a flux spectra from the database!.

    If a energy flux spectra is read then it converts it to photons.

    Parameters
    ----------
    spec_name: str
        Location and name of model spectrum file.

    Returns
    -------
    wav: array-like, float64
        Wavelength in microns.
    flux_photons: array-like, float64
        Photon flux.

    """
    if "photon" in spec_name:
        wav_micron, flux_photons = pdread_2col(spec_name)
    else:
        wav, flux = pdread_2col(spec_name)

        wav_micron = wav * 1.0e-4  # conversion to microns
        flux_photons = flux * wav_micron   # Convert to photons

    return wav_micron, flux_photons


def get_spectrum_name(startype, logg=4.50, feh=0.0, alpha=None, org=False, flux_type="photon"):
    """ Return correct phoenix spectrum filename for a given spectral type.

    Based off phoenix_utils module.

    Parameters
    ----------
    flux_type: str
        Indicate which filetype to try find. e.g. "photon", "wave", ("fits" Not Implemented yet)

    Returns
    -------
    spectrum_name: str
        The name of spectrum with choosen Parameters


    Ability to select logg and metalicity (feh) later on.
    org = original locations (without Z folder option)
    """
    if feh == 0:
        feh = -0.0    # make zero negative to signed integer.

    temps = {"M0": 3900, "M3": 3500, "M6": 2800, "M9": 2600}
    if (flux_type == "photon") and (not org):
        base = "PHOENIX-ACES-AGSS-COND-2011-HiRes_wave_photon.dat"
    else:
        base = "PHOENIX-ACES-AGSS-COND-2011-HiRes_wave.dat"

    if startype in temps.keys():
        if org:
            phoenix_name = "PHOENIX-ACES_spectra/lte{0:05d}-{1}-{2}.{3}".format(temps[startype], "4.50", "0.0", base)
        elif (alpha is not None) and (alpha != 0.0):
            if abs(alpha) > 0.2:
                print("Warning! Alpha is outside acceptable range -0.2->0.2")

            phoenix_name = ("PHOENIX-ACES_spectra/Z{2:+4.1f}.Alpha={3:+5.2f}/"
                            "lte{0:05d}-{1:4.2f}{2:+4.1f}.Alpha={3:+5.2f}.{4:s}"
                            "").format(temps[startype], logg, feh, alpha, base)
        else:
            phoenix_name = ("PHOENIX-ACES_spectra/Z{2:+4.1f}/lte{0:05d}-{1:4.2f}"
                            "{2:+4.1f}.{3:s}").format(temps[startype], logg, feh, base)

        spectrum_name = phoenix_name
    elif re.match(r"^[OBAFGKML][0-9]$", startype):   # Valid spectral types
        raise NotImplementedError("The spectral type '{:s}' is not implemented yet.".format(startype))
    else:
        raise ValueError("'{:s}' is not a valid spectral type.".format(startype))

    return spectrum_name


def band_selector(wav, flux, band):
    """ Select a specific wavelength band.

    Parameters
    ----------
    wav: array-like
        Wavelength values.
    flux: array-like
        Flux values.
    band: str
        Band letter to select, upper or lower case is supported. Options
        are ("ALL" or ""), "VIS", "GAP", "Z", "Y", "J", "H", "K".
    """
    band = band.upper()

    bands = {"VIS": (0.38, 0.78), "GAP": (0.78, 0.83), "Z": (0.83, 0.93),
             "Y": (1.0, 1.1), "J": (1.17, 1.33), "H": (1.5, 1.75),
             "K": (2.07, 2.35), "CONT": (0.45, 1.05), "NIR": (0.83, 2.35)}
    if(band in ["ALL", ""]):
        return wav, flux
    elif band in bands:
        # select values form the band
        bandmin = bands[band][0]
        bandmax = bands[band][1]

        return wav_selector(wav, flux, bandmin, bandmax)
    else:
        print("Unrecognized band tag.")
        exit(1)


def wav_selector(wav, flux, wav_min, wav_max):
    """
    function that returns wavelength and flux withn a giving range

    Parameters
    ----------
    wav: array-like
        Wavelength array.
    flux: array-like
        Flux array.
    wav_min: float
        Lower bound wavelength value.
    wav_max: float
        Upper bound wavelength value.

    Returns
    -------
    wav_sel: array
        New wavelength array within bounds wav_min, wav_max
    flux_sel: array
        New wavelength array within bounds wav_min, wav_max
        """
    wav = np.asarray(wav, dtype="float64")
    flux = np.asarray(flux, dtype="float64")

    mask = (wav > wav_min) & (wav < wav_max)
    flux_sel = flux[mask]
    wav_sel = wav[mask]

    return wav_sel, flux_sel


def unitary_Gauss(x, center, FWHM):
    """ Gaussian function of area = 1.

    Parameters
    ----------
    x: array-like
        Position array
    center: float
        Central postion of Gaussian
    FHWM: float
        Full Width at Half Maximum

    Returns
    -------
    result: array-like
        Result of gaussian function sampled at x values.
    """

    sigma = np.abs(FWHM) / (2 * np.sqrt(2 * np.log(2)))
    Amp = 1.0 / (sigma * np.sqrt(2 * np.pi))
    tau = -((x - center)**2) / (2 * (sigma**2))
    result = Amp * np.exp(tau)

    return result


def rotation_kernel(delta_lambdas, delta_lambda_L, vsini, epsilon):
    """ Calculate the rotation kernel for a given wavelength

    Parameters
    ----------
    delta_lambdas: array
        Wavelength values selected within delta_lambda_L around central value. (check)
    delta_lambda_L: float
        FWHM of rotational broading. (check)
    vsini: float
        Projected rotational velocity [km/s]
    epsilon: float
        Linear limb-darkening coefficient (0-1).

    Returns
    -------
        Rotational kernal

    Notes:
    Equations * from .... book.

    """
    denominator = (np.pi * vsini * (1.0 - epsilon / 3.0))
    lambda_ratio_sqr = (delta_lambdas / delta_lambda_L)**2.0

    c1 = 2.0 * (1.0 - epsilon) / denominator
    c2 = 0.5 * np.pi * epsilon / denominator

    return (c1 * np.sqrt(1.0-lambda_ratio_sqr) + c2 * (1.0-lambda_ratio_sqr))


def plotter(spectrum, band, vsini=0, R=0):
    """
    Reads and plots the selected spectrum in a given band
    """
    wav, flux = read_spectrum(spectrum)
    wav_band, flux_band = band_selector(wav, flux, band)
    flux_counts = flux_band

    plt.figure(1)
    plt.xlabel(r"wavelength [ $\mu$m ])")
    plt.ylabel(r"flux [counts] ")
    plt.plot(wav_band, flux_counts, color='k', marker="o", linestyle="-")
    plt.show()
    plt.close()

    delta_wav = np.array(wav_band[1:]) - np.array(wav_band[:-1])
    plt.figure(1)
    plt.xlabel(r"wavelength [ $\mu$m ])")
    plt.ylabel(r"Step [$\Delta\,\mu$m] ")
    plt.plot(wav_band[1:], delta_wav, color='k', marker="o", linestyle="-")
    plt.show()
    plt.close()

def calculate_ratios(spectrum):
    """ Calculate ratios between the different bands.

    This was labeled WRONG from original code, but including here for reference.
    """
    wav, flux = read_spectrum(spectrum)
    wav_band, flux_band = band_selector(wav, flux, "ALL")
    # This seams wrong as read_spectrum already turns to photons
    flux_band = flux_band * wav_band  # passing it to photon counts
    index = np.searchsorted(wav_band, [0.5, 1.0, 1.5, 2.07])

    flux_0p5 = np.max(flux_band[index[0]-5: index[0]+5])
    flux_right1 = np.max(flux_band[index[1]: index[1]+10])
    flux_right1p5 = np.max(flux_band[index[2]: index[2]+10])
    flux_right2p07 = np.max(flux_band[index[3]: index[3]+10])

    flux_values = [flux_0p5, flux_right1, flux_right1p5, flux_right2p07]
    return flux_values


def list_creator(spectrum, band):
    """
    creates a list of potential lines from a brute-force analysis of the band
    """
    wav, flux = read_spectrum(spectrum)
    wav_band, flux_band = band_selector(wav, flux, band)

    print(band + " band list:")
    short_flux = flux_band[2:-2]
    left_mask = ((short_flux < flux_band[:-4]) &
                 (short_flux < flux_band[1:-3]) &
                 (flux_band[:-4] > flux_band[1:-3]))

    right_mask = ((short_flux < flux_band[3:-1]) &
                  (short_flux < flux_band[4:]) &
                  (flux_band[4:] > flux_band[3:-1]))

    line_centers = wav_band[2:-2][left_mask * right_mask]  # find peaks using masking
    print("Line centers", line_centers * 1.0e4)
    print("In a spectrum with {} points".format(len(wav_band)),
          ", {} lines were found.".format(len(line_centers)))
    return line_centers


def silentremove(filename):
    """ Remove file without failing when it doesn't exist."""
    try:
        os.remove(filename)
    except OSError as e:  # this would be "except OSError, e:" before Python 2.6
        if e.errno != errno.ENOENT:  # errno.ENOENT = no such file or directory
            raise  # re-raise exception if a different error occured

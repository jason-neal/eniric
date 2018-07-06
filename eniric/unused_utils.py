# Old way things were done.


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
    """Calculate ratios between the different bands.

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
    print("In a spectrum with {0} points".format(len(wav_band)),
          ", {0} lines were found.".format(len(line_centers)))
    return line_centers

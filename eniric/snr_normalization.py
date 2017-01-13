
""" The photon flux can be scaled by multiplicative constant which affects the SNR of the spectra and the radial velocity precision.

For consistancy and valid comparision we normalize each spectra to acheive a consistant SNR at a specific location. """

# Normaize to SNR 100 in middle of J band 1.25 micron!

""" This is the original values to acheive a SNR of 100 in the middle of the J band at 1.25 micron """
def normalize_flux(flux_stellar, id_string):
    """Normalize flux to have SNR of 100 in middle of J band."""

    if("M0" in id_string):
        norm_constant = 1607

    elif("M3" in id_string):
        norm_constant = 1373

    elif("M6" in id_string):
        if("1.0" in id_string):
            norm_constant = 933
        elif("5.0" in id_string):
            norm_constant = 967
        else:
            norm_constant = 989

    elif("M9" in id_string):
        if("1.0" in id_string):
            norm_constant = 810
        elif("5.0" in id_string):
            norm_constant = 853
        else:
            norm_constant = 879
    else:
        print("Constant not defined. Aborting...")
        exit()

    return flux_stellar / ((norm_constant / 100.0)**2.0)


def determine_constant(wav, flux, SNR=100, band="J"):
    """ Determine the constant to acheive a SNR in the middle of band."

    wav: ndarray (micron)
    flux: ndarray (photons/s/cm**2)
    SNR: int,  default = 100
    band: str, default = "J"
    """
    band_mid = {"J": 1.25}
    # Option to specify own wavelength?
    index_reference = np.searchsorted(wav, band_mid[band])  # Searching for the closest index

    SN_estimate = np.sqrt(np.sum(flux_stellar[index_reference-1:index_reference+2]))

    print("\tSanity Check: The S/N for the {:s} reference model was of {:4.2f}.".format(id_string, SN_estimate))

    norm_constant = (SN_estimate / SNR)**2

    # Test it
    flux2 = flux / norm_constant
    SN_estimate2 = np.sqrt(np.sum(flux2[index_reference-1:index_reference+2]))
    print("\tSanity Check: The S/N for the {:s} normalized result was of {:4.2f}.".format(id_string, SN_estimate))

    return norm_constant

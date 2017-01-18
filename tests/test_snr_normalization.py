import pytest
import numpy as np
import eniric.utilities as utils

# Test using hypothesis
# from hypothesis import given, example
# import hypothesis.strategies as st
from eniric.snr_normalization import snr_constant_band, snr_constant_wav, sampling_index

def test_snr_normalization():
    """ Test SNR after normalizing function is the desired value.
    Testing on middle pf J band."""

    test_data = ("data/PHOENIX-ACES_spectra/Z-0.0/lte02800-4.50"
                "-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes_wave.dat")

    band = "J"
    band_mid = {"J": 1.25}
    wav, flux = utils.read_spectrum(test_data)

    for desired_snr in [50.0, 100.0, 150.0, 200.0]:

        index_reference = np.searchsorted(wav, band_mid[band])  # Searching for the closest index to 1.25
        SN_estimate = np.sqrt(np.sum(flux[index_reference-1:index_reference+2]))

        assert round(SN_estimate, 0) != desired_snr     # Assert SNR is not correct

        norm_const = snr_constant_band(wav, flux, SNR=desired_snr, band=band)

        new_flux = flux / norm_const
        new_SN_estimate = np.sqrt(np.sum(new_flux[index_reference-1:index_reference+2]))

        assert round(new_SN_estimate, 0) == desired_snr

        assert snr_constant_band(wav, flux, band="J", SNR=desired_snr) == snr_constant_wav(wav, flux, 1.25, SNR=desired_snr)


def test_band_snr_norm():
    """ Compared to wav snr norm"""
    # snr_constant_band
    test_data = ("data/PHOENIX-ACES_spectra/Z-0.0/lte02800-4.50"
                "-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes_wave.dat")
    wav, flux = utils.read_spectrum(test_data)
    assert snr_constant_band(wav, flux, band="J", SNR=100) == snr_constant_wav(wav, flux, wav_snr=1.25, SNR=100)



def test_sampling_index():

    # Some hard coded examples
    # odd number
    assert sampling_index(6, 1) == [6]
    assert np.all(sampling_index(100, 3) == [99, 100, 101])
    assert np.all(sampling_index(5, 5) == [3, 4, 5, 6, 7])
    # even number
    assert np.all(sampling_index(10, 4) == [8, 9, 10, 11])
    assert np.all(sampling_index(10, 2) == [9, 10])
    # number that is at end of array.

    # Check number of values correct.
    for sample in [1, 2, 3, 4, 5, 7, 10, 15]:
        assert len(sampling_index(20, sample)) == sample
        assert type(sampling_index(20, sample)[0]) == np.int64  # index values must be int


def test_sampling_index_array():
    """ Sampling index when array_length is given, or when index goes out of bounds."""
    assert np.all(sampling_index(100, 3, array_length=200) == [99, 100, 101])

    with pytest.raises(ValueError):
        sampling_index(3, 10)   # index will be < 0
    with pytest.raises(ValueError):
        sampling_index(3, 9, array_length=50)   # index will be < 0
    with pytest.raises(ValueError):
        sampling_index(46, 10, array_length=50)   # an index will be > (array_length - 1)

import numpy as np
import pytest
import os

import eniric.IOmodule as Io
import eniric.Qcalculator as Q
# Test using hypothesis
# from hypothesis import given, example
# import hypothesis.strategies as st
import eniric.snr_normalization as snrnorm
import eniric.utilities as utils
import eniric

file_error_to_catch = getattr(__builtins__, 'FileNotFoundError', IOError)


@pytest.mark.xfail(raises=file_error_to_catch)
def test_snr_normalization():
    """Test SNR after normalizing function is the desired value.

    Testing on middle of J band.
    """
    test_data = os.path.join(eniric.paths["phoenix_dat"],
        "Z-0.0/lte02800-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes_wave.dat")

    band = "J"
    band_mid = {"J": 1.25}
    wav, flux = utils.read_spectrum(test_data)

    for desired_snr in [50.0, 100.0, 150.0, 200.0]:

        index_reference = np.searchsorted(wav, [band_mid[band]])[0]  # Searching for the closest index to 1.25
        snr_estimate = np.sqrt(np.sum(flux[index_reference - 1:index_reference + 2]))

        assert round(snr_estimate, 0) != desired_snr     # Assert SNR is not correct

        norm_const = snrnorm.snr_constant_band(wav, flux, snr=desired_snr, band=band)

        new_flux = flux / norm_const
        new_snr_estimate = np.sqrt(np.sum(new_flux[index_reference - 1:index_reference + 2]))

        assert round(new_snr_estimate, 0) == desired_snr

        assert (snrnorm.snr_constant_band(wav, flux, band="J", snr=desired_snr) ==
                snrnorm.snr_constant_wav(wav, flux, 1.25, snr=desired_snr))


@pytest.mark.xfail(raises=file_error_to_catch)
def test_band_snr_norm():
    """Compared to wav snr norm."""
    # snr_constant_band
    test_data = os.path.join(
        eniric.paths["resampled"], "Spectrum_M0-PHOENIX-ACES_Jband_vsini1.0_R100k_res3.txt")
    wav, flux = Io.pdread_2col(test_data)

    assert (snrnorm.snr_constant_band(wav, flux, band="J", snr=100) ==
            snrnorm.snr_constant_wav(wav, flux, wav_ref=1.25, snr=100))

    assert (snrnorm.snr_constant_band(wav, flux, band="J", snr=100) !=
            snrnorm.snr_constant_wav(wav, flux, wav_ref=1.24, snr=100))


def test_sampling_index():
    """Some hard coded examples of sampling index."""
    # odd number
    assert snrnorm.sampling_index(6, 1) == [6]
    assert np.all(snrnorm.sampling_index(100, 3) == [99, 100, 101])
    assert np.all(snrnorm.sampling_index(5, 5) == [3, 4, 5, 6, 7])

    # even number
    assert np.all(snrnorm.sampling_index(10, 4) == [8, 9, 10, 11])
    assert np.all(snrnorm.sampling_index(10, 2) == [9, 10])
    # number that is at end of array.

    # Check number of values correct.
    for sample in [1, 2, 3, 4, 5, 7, 10, 15]:
        assert len(snrnorm.sampling_index(20, sample)) == sample
        assert type(snrnorm.sampling_index(20, sample)[0]) == np.int64  # index values must be int


def test_sampling_index_array():
    """Sampling index when array_length is given, or when index goes out of bounds."""
    assert np.all(snrnorm.sampling_index(100, 3, array_length=200) == [99, 100, 101])

    with pytest.raises(ValueError):
        snrnorm.sampling_index(3, 10)   # index will be < 0
    with pytest.raises(ValueError):
        snrnorm.sampling_index(3, 9, array_length=50)   # index will be < 0
    with pytest.raises(ValueError):
        snrnorm.sampling_index(46, 10, array_length=50)   # an index will be > (array_length - 1)


@pytest.mark.parametrize("bad_string",
    ["id-string", "M0-K-1.0-100", "M0-P-1.0-100k"])
def test_errors_in_snr_get_reference_spectrum(bad_string):
    """Testing Eorros in getting the reference spectrum."""
    with pytest.raises(ValueError):
        snrnorm.get_reference_spectrum(bad_string)


@pytest.mark.parametrize("bad_string", ["Alpha=", "smpl="])
def test_notimplemented_errors_in_snr_get_reference_spectrum(bad_string):
    """Testing getting the reference spectrum.

    Currently "Alpha=" in the id-string is not implemented.
    Currently "smpl=" in the id-string is not implemented.
    """
    with pytest.raises(NotImplementedError):
        snrnorm.get_reference_spectrum(bad_string)


@pytest.mark.xfail(raises=file_error_to_catch)
def test_valid_snr_get_reference_spectrum():
    """Testing getting the reference spectrum."""
    ref_band = "J"
    wav_ref, flux_ref = snrnorm.get_reference_spectrum("M0-K-1.0-100k", ref_band=ref_band)
    band_min, band_max = utils.band_limits(ref_band)

    # Test the wavelength is in the reference band wavelength range
    assert np.all(wav_ref <= band_max)
    assert np.all(wav_ref >= band_min)

    # test properties of output
    assert len(wav_ref) == len(flux_ref)
    assert isinstance(wav_ref, np.ndarray)
    assert isinstance(flux_ref, np.ndarray)


@pytest.mark.xfail()  # size is too big
def test_normalize_flux_new_verse_old():
    test_data = os.path.join(eniric.paths["resampled"], "Spectrum_M0-PHOENIX-ACES_Kband_vsini5.0_R100k_res3.txt")
    id_string = "M0-K-5.0-100k"
    wav, flux = utils.read_spectrum(test_data)
    new_norm = snrnorm.normalize_flux(flux, id_string, new=True)
    old_norm = snrnorm.normalize_flux(flux, id_string, new=False)

    rvprec_new = Q.RVprec_calc(wav, new_norm)
    rvprec_old = Q.RVprec_calc(wav, old_norm)

    assert abs(rvprec_new.value - rvprec_old.value) < 1


def test_old_does_does_not_handle_changed_band():
    test_data = os.path.join(eniric.paths["resampled"], "Spectrum_M0-PHOENIX-ACES_Kband_vsini5.0_R100k_res3.txt")
    id_string = "M0-K-5.0-100k"
    wav, flux = utils.read_spectrum(test_data)
    with pytest.raises(ValueError):
        snrnorm.normalize_flux(flux, id_string, new=False, ref_band="K")

    with pytest.raises(ValueError):
        snrnorm.normalize_flux(flux, id_string, new=False, snr=101)

@pytest.mark.parametrize("id_string", [
    "M0-1.0", "M3-1.0", "M6-1.0", "M9-1.0", "M0-5.0", "M3-5.0", "M6-5.0",
    "M9-5.0", "M0-10.0", "M3-10.0", "M6-10.0", "M9-10.0"])
def test_snr_old_norm_constant(id_string):
    norm_const = snrnorm.old_norm_constant(id_string)
    assert isinstance(norm_const, float)


@pytest.mark.parametrize("bad_string", [
    "M0-1", "M0-2.5", "M8-1.0", "M6-5", "M9-10", "T0-3.0", "", "AB-CDE"])
def test_snr_old_norm_constant_with_bad_id_str(bad_string):
    """Fixed to the set of values in first paper."""
    with pytest.raises(ValueError):
        snrnorm.old_norm_constant(bad_string)

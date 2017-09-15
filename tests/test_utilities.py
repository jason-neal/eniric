"""Test utilities for eniric."""

from __future__ import division, print_function

# Test using hypothesis
import hypothesis.strategies as st
import numpy as np
import pytest
from hypothesis import given, settings

import eniric.utilities as utils
import eniric

# For python2.X compatibility
file_error_to_catch = getattr(__builtins__, 'FileNotFoundError', IOError)


@pytest.mark.xfail(raises=file_error_to_catch)
def test_read_spectrum():
    """Test reading in a _wave_photon.dat is the same as a _wave.dat."""
    photon = "data/test_data/sample_lte03900-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes_wave_photon.dat"
    wave = "data/test_data/sample_lte03900-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes_wave.dat"
    wave_wav, wave_flux = utils.read_spectrum(wave)
    photon_wav, photon_flux = utils.read_spectrum(photon)

    assert np.allclose(photon_wav, wave_wav)
    assert np.allclose(photon_flux, wave_flux)


@pytest.mark.xfail(raises=file_error_to_catch)
def test_get_spectrum_name():
    """Test specifing file names with stellar parameters."""
    test = ("Z-0.0/lte02800-4.50"
            "-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes_wave.dat")

    assert utils.get_spectrum_name("M6", flux_type="wave") == test

    test_alpha = ("Z-0.0.Alpha=+0.20/"
                  "lte02600-6.00-0.0.Alpha=+0.20.PHOENIX-ACES-AGSS-COND-2011-HiRes_wave_photon.dat")
    assert utils.get_spectrum_name("M9", logg=6, alpha=0.2) == test_alpha

    test_pos_feh = ("Z+0.5/"
                    "lte03500-0.00+0.5.PHOENIX-ACES-AGSS-COND-2011-HiRes_wave_photon.dat")
    assert utils.get_spectrum_name("M3", logg=0, feh=0.5, alpha=0.0) == test_pos_feh

    test_photon = ("Z-0.0/lte02800-4.50"
                   "-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes_wave_photon.dat")
    assert utils.get_spectrum_name("M6") == test_photon

    # Catch Errors
    with pytest.raises(NotImplementedError):
        utils.get_spectrum_name("K1")       # Stellar type not added
    with pytest.raises(NotImplementedError):
        utils.get_spectrum_name("A8")       # Stellar type not added
    with pytest.raises(ValueError):
        utils.get_spectrum_name("MO")       # Miss spelled M0
    with pytest.raises(ValueError):
        utils.get_spectrum_name("X10")      # Not valid spectral type in [OBAFGKML]


@pytest.mark.parametrize("bad_alpha", [-0.3, 0.3, 1])
def test_spectrum_name_with_bad_alpha(bad_alpha):
    """Bad_alpha is outside range -0.2-0.2 for Mdwarf science case."""
    with pytest.raises(ValueError):
        utils.get_spectrum_name("M0", alpha=bad_alpha)


@pytest.mark.parametrize("alpha", [-0.2, 0.1, 0.2])
def test_spectrum_name_with_ok_alpha(alpha):
    name = utils.get_spectrum_name("M0", alpha=alpha)

    assert isinstance(name, str)
    assert str(alpha) in name
    assert "Alpha=" in name


@pytest.mark.xfail(raises=file_error_to_catch)
def test_org_name():
    """Test org flag of utils.get_spectrum_name, suposed to be temporary."""
    test_org = "lte03900-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes_wave.dat"
    assert utils.get_spectrum_name("M0", org=True) == test_org


@given(st.lists(st.floats(allow_nan=False, allow_infinity=False)),
       st.floats(allow_nan=False, allow_infinity=False),
       st.floats(allow_nan=False, allow_infinity=False),
       st.floats(allow_nan=False, allow_infinity=False))
def test_wav_selector(x, y, wav_min, wav_max):
    """Test some properties of wavelength selector."""
    y = [xi + y for xi in x]   # just to make y different
    x1, y1 = utils.wav_selector(x, y, wav_min, wav_max)

    # All values in selected should be less than the max and greater
    # than the min value.
    assert all(x1 >= wav_min)
    assert all(x1 <= wav_max)
    assert len(x1) == len(y1)
    assert isinstance(x1, np.ndarray)
    assert isinstance(y1, np.ndarray)


@pytest.mark.parametrize("band", [
    "VIS", "GAP", "z", "Y", "h", "J", "K", "CONT", "NIR"])
def test_band_limits(band):
    """Test getting limits out of band."""
    band_min, band_max = utils.band_limits(band)

    assert band_min < band_max
    assert band_min != band_max
    assert 0 < band_min < 3
    assert 0 < band_max < 3


@pytest.mark.parametrize("band", ["Z", "H", "J", "K"])
def test_band_selector(band):
    """Test band selector selects the wav and flux in the given band."""
    wav = np.linspace(0.5, 3, 100)
    flux = wav**2

    band_min, band_max = utils.band_limits(band)
    assert not np.all(wav > band_min)      # Assert wav goes outside band
    assert not np.all(wav < band_max)

    wav, flux = utils.band_selector(wav, flux, band)
    assert np.all(wav > band_min)
    assert np.all(wav < band_max)


@pytest.mark.parametrize("band,error",[
    ("X", ValueError),
    ("M0", ValueError),
    (1, AttributeError),
    (np.array(1), AttributeError),
    (["list", "of", "strings"], AttributeError),
])
def test_band_limits_raises_errors(band, error):
    """Test it raises the Value and Attribute Errors."""
    with pytest.raises(error):
        utils.band_limits(band)


@pytest.mark.parametrize("band,error",[
    ("X", ValueError),
    ("M0", ValueError),
    (1, AttributeError),
    (["list", "of", "strings"], AttributeError),
    (np.linspace(1,2,10), AttributeError)
])
def test_band_selector_raises_errors(band, error):
    """Test it raises the Value and Attribute Errors"""
    wav = np.linspace(0.5, 3, 100)
    flux = wav ** 2

    with pytest.raises(error):
        utils.band_selector(wav, flux, band)


@pytest.mark.parametrize("band", ["ALL", ""])
def test_band_selector_with_no_selection(band):
    """If band = "ALL" or ""."""
    wav = np.linspace(0.5, 3, 100)
    flux = wav ** 2
    wav2, flux2 = utils.band_selector(wav, flux, band)

    # No changes
    assert np.all(wav == wav2)
    assert np.all(flux2 == flux)


@settings(max_examples=100)
@given(st.lists(st.floats(min_value=1e-7, max_value=1e-5, allow_infinity=False,
       allow_nan=False), unique=True, min_size=3, max_size=25),
       st.floats(min_value=1e-2, max_value=200), st.floats(min_value=1e-4, max_value=1))
def test_rotational_kernal(delta_lambdas, vsini, epsilon):
    """Test that the new and original code produces the same output."""
    delta_lambdas = np.sort(np.asarray(delta_lambdas), kind='quicksort')
    delta_lambdas = np.append(np.flipud(delta_lambdas), np.insert(delta_lambdas, 0, 0))
    delta_lambda_l = np.max(delta_lambdas) * 2

    new_profile = utils.rotation_kernel(delta_lambdas, delta_lambda_l, vsini, epsilon)

    assert len(new_profile) == len(delta_lambdas)
    # other properties to test?

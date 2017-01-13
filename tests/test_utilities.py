"""test_eniric.py"""

import pytest
import numpy as np
from eniric.utilities import get_spectrum_name, wav_selector

# Test using hypothesis
from hypothesis import given, example
import hypothesis.strategies as st


def test_get_spectrum_name():
    """ Test specifing file names with stellar parameters."""
    test = ("PHOENIX-ACES_spectra/Z-0.0/lte02800-4.50"
            "-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes_wave.dat")
    assert get_spectrum_name("M6") == test

    test_alpha = ("PHOENIX-ACES_spectra/Z-0.0.Alpha=+0.20/"
                  "lte02600-6.00-0.0.Alpha=+0.20.PHOENIX-ACES-AGSS-COND-2011-HiRes_wave.dat")
    assert get_spectrum_name("M9", logg=6, alpha=0.2) == test_alpha

    test_pos_feh = ("PHOENIX-ACES_spectra/Z+0.5/"
                    "lte03500-0.00+0.5.PHOENIX-ACES-AGSS-COND-2011-HiRes_wave.dat")
    assert get_spectrum_name("M3", logg=0, feh=0.5, alpha=0.0) == test_pos_feh

    # Catch Errors
    with pytest.raises(NotImplementedError):
        get_spectrum_name("K1")       # Stellar type not added
    with pytest.raises(NotImplementedError):
        get_spectrum_name("A8")       # Stellar type not added
    with pytest.raises(ValueError):
        get_spectrum_name("MO")       # Miss spelled M0
    with pytest.raises(ValueError):
        get_spectrum_name("X10")      # Not valid spectral type in [OBAFGKML]


def test_org_name():
    """ Test org flag of get_spectrum_name, suposed to be temporary."""
    test_org = "PHOENIX-ACES_spectra/lte03900-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes_wave.dat"
    assert get_spectrum_name("M0", org=True) == test_org


@given(st.lists(st.floats()), st.floats(), st.floats(), st.floats())
def test_wav_selector(x, y, wav_min, wav_max):
    """Test some properties of wavelength selector"""
    y = [xi + y for xi in x]   # just to make y different
    x1, y1 = wav_selector(x, y, wav_min, wav_max)

    # All values in selected should be less than the max and greater
    # than the min value.
    assert all(x1 >= wav_min)
    assert all(x1 <= wav_max)
    assert len(x1) == len(y1)
    assert isinstance(x1, np.ndarray)
    assert isinstance(y1, np.ndarray)

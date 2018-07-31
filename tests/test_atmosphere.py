"""Test of atmosphere.py functions."""
import os

import numpy as np
import pytest
from hypothesis import assume, given, strategies as st

import eniric
from eniric.atmosphere import Atmosphere, consecutive_truths


def test_consecutive_truths():
    """Test consecutive truths lists count of consecutive ones."""
    array1 = np.array([1, 1, 1, 0, 1, 1, 0, 0, 1, 0], dtype=bool)  # [3,2 1]
    array2 = np.array([0, 0, 1, 1, 1, 1, 1, 0, 1, 1], dtype=bool)  # [5, 2]
    array3 = np.array([1, 0, 1, 0, 1, 0, 1, 0, 1, 0], dtype=bool)  # [1,1,1,1,1]
    array4 = np.zeros(5, dtype=bool)  # []
    array5 = np.ones(5, dtype=bool)  # [5]

    assert np.all(consecutive_truths(array1) == [3, 2, 1])
    assert np.all(consecutive_truths(array2) == [5, 2])
    assert np.all(consecutive_truths(array3) == [1, 1, 1, 1, 1])
    assert np.all(consecutive_truths(array5) == [5])
    assert np.all(consecutive_truths(array4) == [])

    # Check sum of trues equal sum of returned list
    assert np.sum(array1) == np.sum(consecutive_truths(array1))

    rand_array = np.asarray(
        np.floor(np.random.random(50) * 2), dtype=bool
    )  # random values
    assert np.sum(rand_array) == np.sum(consecutive_truths(rand_array))


@pytest.mark.xfail()
def test_barycenter_shift():
    """Test barycentric shift code."""
    raise False


@pytest.mark.xfail()
def test_old_barycenter_shift():
    """Test old barycentric shift code."""
    raise False


@pytest.mark.xfail()  # If missing the datafiles
def test_prepare_atmosphere():
    """Test that an atmosphere file is loaded.

    Test band is close to the band limits (accounting for the given offset. +- 100km/s).
    """
    raise False


# If I multiply the flux by the mask then the smallest flux should be 0.98
# If not then there would be an issue.


size = 100


@given(size=st.integers(min_value=10, max_value=50))
@given(
    st.lists(st.floats(min_value=1, max_value=3000), min_size=size, max_size=size),
    st.lists(st.floats(min_value=0, max_value=1), min_size=size, max_size=size),
    st.lists(st.booleans(), min_size=size, max_size=size),
)
def test_Amosphere_class(wave, transmission, mask):
    atmos = Atmosphere(wave, transmission, mask)
    assert np.all(atmos.wavelength == wave)
    assert np.all(atmos.transmission == transmission)
    assert np.all(atmos.wavelength == mask)
    assert atmos.mask.dtype == np.bool


def test_Amosphere_class_nomask(wave, transmission):
    atmos = Atmosphere(wave, transmission)
    assert np.all(atmos.wavelength == wave)
    assert np.all(atmos.transmission == transmission)
    assert np.all(atmos.mask == 1)
    assert len(atmos.mask) == len(atmos.transmssion)
    assert atmos.mask.dtype == np.bool


@pytest.fixture(
    params=[
        "Average_TAPAS_2014_H.txt",
        "Average_TAPAS_2014_K.txt",
        "Average_TAPAS_2014_J.txt",
    ]
)
def atm_model(request):
    """Get atmospheric model name to load."""
    return os.path.join(eniric.paths["atmmodel"], request.param)


def test_atmosphere_from_file(atm_model):
    atmos = Atmosphere._from_file(atm_model)
    assert len(atmos.self) == len(atmos.transmssion)
    assert len(atmos.transmssion[atmos.mask]) != len(
        atmos.transmssion
    )  # mask is not all ones


# @given(size=st.integers(min_value=10, max_value=100))
@given(st.lists(st.floats(min_value=0, max_value=1), min_size=size, max_size=size))
def test_atmosphere_from_file(atmmodel):
    atmos = Atmosphere._from_file(atmmodel)
    assert len(atmos.self) == len(atmos.transmssion)
    assert len(atmos.transmssion[atmos.mask]) != len(
        atmos.transmssion
    )  # mask is not all ones


# @given(size=st.integers(min_value=10, max_value=100))
@given(
    st.lists(st.floats(min_value=1, max_value=3000), min_size=size, max_size=size),
    st.lists(st.floats(min_value=0, max_value=1), min_size=size, max_size=size),
)
@given(percent=st.floats(min_value=0.5, max_value=99.9))
def test_atmosphere_masking(wav, trans, percent):
    cutoff = 1 - (percent / 100.0)
    assume(np.any(trans < cutoff))
    assume(np.any(trans > cutoff))
    atmos = Atmosphere(wavelength=wav, transmission=trans)
    org_mask = atmos.mask
    atmos.mask_transmission(percent)
    assert np.any(atmos.mask != org_mask)  # mask has changed
    assert np.all(atmos.transmission[atmos.mask] < cutoff)


@pytest.mark.xfail()
def test_values_within_rv_of_tell_line_are_masked():
    x, y, z = 0, 0, 0
    # Enable object
    atmos = Atmosphere(x, y, z)
    # mask pixels
    atmos.mask_transmission()
    # RV shift mask
    # ...

    for pixel, mask_value in zip(atmos.wavelength, atmos.mask):
        if mask_value != 0:
            # Already masked out
            pass
        else:
            # Find rv limits to this pixel.
            wl_lower, wl_upper = pixel * (1 - 30 / 3e5), pixel * (1 + 30 / 3e5)
        # Check mask is wider
        wl_mask = (atmos.wavlength >= wl_lower) * (atmos.wavelength < wl_upper)
        assert np.all(atmos.mask[wl_mask] == 1)

    assert False

# todo
# test atm masking
# test rv usage
# test input into rv work

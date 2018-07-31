"""Test of atmosphere.py functions."""
import os

import numpy as np
import pytest
from hypothesis import assume, given, strategies as st

import eniric
from eniric.atmosphere import Atmosphere, barycenter_shift, consecutive_truths


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


@pytest.mark.xfail()  # If missing the datafiles
def test_prepare_atmosphere():
    """Test that an atmosphere file is loaded.

    Test band is close to the band limits (accounting for the given offset. +- 100km/s).
    """
    raise False


# If I multiply the flux by the mask then the smallest flux should be 0.98
# If not then there would be an issue.


size = 10


@given(
    st.lists(st.floats(min_value=1, max_value=3000), min_size=size, max_size=size),
    st.lists(st.floats(min_value=0, max_value=1), min_size=size, max_size=size),
    st.lists(st.booleans(), min_size=size, max_size=size),
)
def test_atmosphere_class(wave, transmission, mask):
    atmos = Atmosphere(np.array(wave), np.array(transmission), np.array(mask))
    assert np.all(atmos.wl == wave)
    assert np.all(atmos.transmission == transmission)
    assert np.all(atmos.mask == mask)
    assert atmos.mask.dtype == np.bool


@given(
    st.lists(st.floats(min_value=1, max_value=3000), min_size=size, max_size=size),
    st.lists(st.floats(min_value=0, max_value=1), min_size=size, max_size=size),
    st.lists(st.booleans(), min_size=size, max_size=size),
)
def test_atmosphere_class_turns_lists_to_arrays(wave, transmission, mask):
    atmos = Atmosphere(wave, transmission, mask)
    assert isinstance(atmos.wl, np.ndarray)
    assert isinstance(atmos.transmission, np.ndarray)
    assert isinstance(atmos.mask, np.ndarray)
    assert atmos.mask.dtype == np.bool


@given(
    st.lists(st.floats(min_value=1, max_value=3000), min_size=size, max_size=size),
    st.lists(st.floats(min_value=0, max_value=1), min_size=size, max_size=size),
)
def test_Amosphere_class_nomask(wave, transmission):
    atmos = Atmosphere(wave, transmission)
    assert np.all(atmos.wl == wave)
    assert np.all(atmos.transmission == transmission)
    assert np.all(atmos.mask == 1)
    assert len(atmos.mask) == len(atmos.transmission)
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


@pytest.fixture(params=[1, 2, 20])
def atmosphere_fixture(request, atm_model):
    percent_cutoff = request.param
    atm = Atmosphere._from_file(atm_model)
    atm.mask_transmission(percent_cutoff)
    return atm


def test_atmosphere_from_file(atm_model):
    atmos = Atmosphere._from_file(atmmodel=atm_model)
    print(atmos)
    print(atmos.transmission)
    print(atmos.wl)
    assert len(atmos.wl) == len(atmos.transmission)
    assert len(atmos.transmission[atmos.mask]) != len(
        atmos.transmission
    )  # mask is not all ones


@given(
    trans=st.lists(st.floats(min_value=0, max_value=1), min_size=size, max_size=size),
    percent=st.floats(min_value=0.5, max_value=99),
)
def test_atmosphere_masking(trans, percent):
    cutoff = 1 - (percent / 100.0)

    assume(np.any(np.array(trans) < cutoff))
    assume(np.any(np.array(trans) > cutoff))

    atmos = Atmosphere(wavelength=np.arange(len(trans)), transmission=trans)
    org_mask = atmos.mask
    atmos.mask_transmission(percent)

    assert np.any(atmos.mask != org_mask)  # mask has changed
    assert np.all(atmos.transmission[atmos.mask] < cutoff)


@pytest.mark.xfail()
@given(rv=st.floats(min_value=-2000, max_value=2000))
def test_values_within_rv_of_tell_line_are_masked(atmosphere_fixture):
    # Enable object
    atmos = atmosphere_fixture
    org_mask = atmos.mask.copy()
    # RV shift mask
    atmos.bary_shift_mask()
    # wav_lower = self.wl - delta_lambdas
    # wav_upper = self.wl + delta_lambdas

    for pixel, mask_value, org_val in zip(atmos.wl, atmos.mask, org_mask):
        if mask_value != 0:
            # Already masked out
            pass
        else:
            # Find rv limits to this pixel.
            wl_lower, wl_upper = pixel * (1 - 30 / 3e5), pixel * (1 + 30 / 3e5)
        # Check mask is wider
        wl_mask = (atmos.wl >= wl_lower) * (atmos.wl < wl_upper)
        assert np.all(atmos.mask[wl_mask] == 1)

    assert False


@pytest.mark.xfail()
def test_atmos_barycenter_shift_mask():
    """Test barycentric shift code."""
    raise False


def test_barycenter_shift_verse_class(atmosphere_fixture):
    """Test barycentric shift code."""
    atmos = atmosphere_fixture

    mask30kms = barycenter_shift(
        atmos.wl, atmos.transmission, atmos.mask, consecutive_test=False
    )

    assert not np.allclose(mask30kms, atmos.mask)
    atmos.bary_shift_mask()

    assert np.all_close(atmos.mask)


# todo
# test atm masking
# test rv usage
# test input into rv work

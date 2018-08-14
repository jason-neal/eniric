"""Test of atmosphere.py functions."""
import os

import numpy as np
import pytest
from hypothesis import assume, given, strategies as st

import eniric
from eniric.atmosphere import Atmosphere, barycenter_shift, consecutive_truths
from eniric.broaden import resolution_convolution
from eniric.utilities import band_limits

size = 50  # Size of arrays if need consistent length


@pytest.mark.xfail()
def test_Atmosphere_funtional_test(short_atmosphere):
    """ Check Mask still works after shifting and masks out correct values.
    Load in telluric file.
    Mask the lines,
    Bary shift the lines.
    Check masks are still ok.
    """
    atm = short_atmosphere
    atm.mask = np.ones_like(atm.wl, dtype=bool)  # reset mask
    assert not np.all(atm.transmission[atm.mask] >= 0.98)
    atm.broaden(100000)
    atm.mask_transmission(2)
    assert np.all(atm.transmission[atm.mask] >= 0.98)
    atm.bary_shift_mask(consecutive_test=True)
    assert np.all(atm.transmission[atm.mask] >= 0.98)


@pytest.mark.parametrize(
    "wave, transmission, std, mask",
    [
        ([1, 2, 3, 4], [.5, .6, .7, 8], [0., 0.2, 0.1, 0.2], [0, 1, 1, 1]),
        ([2000, 2100], [0.97, 0.99], [0.1, 0.1], [False, True]),
        ([], [], [], []),
    ],
)
def test_atmosphere_class(wave, transmission, std, mask):
    atmos = Atmosphere(
        np.array(wave), np.array(transmission), np.array(mask), std=np.array(std)
    )
    assert np.all(atmos.wl == wave)
    assert np.all(atmos.transmission == transmission)
    assert np.all(atmos.mask == mask)
    assert np.all(atmos.std == std)
    assert atmos.mask.dtype == np.bool


@pytest.mark.parametrize(
    "wave, transmission, std, mask ", [([1, 2, 3], [4, 5, 6], [7, 8, 9], [10, 11, 12])]
)
def test_atmosphere_class_turns_lists_to_arrays(wave, transmission, std, mask):
    atmos = Atmosphere(wave, transmission, mask, std)
    assert isinstance(atmos.wl, np.ndarray)
    assert isinstance(atmos.transmission, np.ndarray)
    assert isinstance(atmos.std, np.ndarray)
    assert isinstance(atmos.mask, np.ndarray)
    assert atmos.mask.dtype == np.bool


@pytest.mark.parametrize(
    "wave, transmission", [([1, 2, 3], [0.4, 0.5, 0.6]), ([7, 8, 9], [.10, .11, .12])]
)
def test_atmosphere_class_nomask(wave, transmission):
    atmos = Atmosphere(wave, transmission)
    assert np.all(atmos.wl == wave)
    assert np.all(atmos.transmission == transmission)
    assert np.all(atmos.mask == 1)
    assert len(atmos.mask) == len(atmos.transmission)
    assert atmos.mask.dtype == np.bool


@pytest.fixture(
    params=[
        # "Average_TAPAS_2014_H.txt",
        "Average_TAPAS_2014_K.txt",
        "Average_TAPAS_2014_J.txt",
    ]
)
def atm_model(request):
    """Get atmospheric model name to load."""
    return os.path.join(eniric.paths["atmmodel"], request.param)


@pytest.fixture(params=[2.5, 4])
def atmosphere_fixture(request, atm_model):
    percent_cutoff = request.param
    atm = Atmosphere.from_file(atm_model)
    atm.mask_transmission(percent_cutoff)
    return atm


@pytest.fixture()
def short_atmosphere(atmosphere_fixture):
    # First 2000 data points only to speed up tests
    return atmosphere_fixture[:2000]


@pytest.fixture(params=[(0, 1500), (8000, 9000)])
def sliced_atmmodel_default_mask(request, atm_model):
    """To do own masking. Sliced in different places."""
    lower, upper = request.param  # slice limits
    atm = Atmosphere.from_file(atm_model)
    return atm[int(lower) : int(upper)]


def test_atmosphere_from_file(atm_model):
    atmos = Atmosphere.from_file(atmmodel=atm_model)
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
    assume(np.any(np.array(trans) >= cutoff))

    atmos = Atmosphere(wavelength=np.arange(len(trans)), transmission=trans)
    org_mask = atmos.mask
    atmos.mask_transmission(percent)

    assert np.any(atmos.mask != org_mask)
    assert np.all(atmos.transmission[atmos.mask] >= cutoff)
    assert np.all(atmos.transmission[~atmos.mask] < cutoff)


@pytest.mark.parametrize("rv", [-1471, 65])  # some fixed RV values
def test_values_within_the_rv_of_telluric_lines_are_masked(
    sliced_atmmodel_default_mask, rv
):
    # Enable object
    atmos = sliced_atmmodel_default_mask
    org_mask = atmos.mask.copy()
    # RV shift mask
    atmos.bary_shift_mask(rv=rv, consecutive_test=False)

    for pixel, mask_value, org_val in zip(atmos.wl, atmos.mask, org_mask):
        if mask_value == 0:
            # Already masked out
            pass
        else:
            # Find rv limits to this pixel.
            wl_lower, wl_upper = pixel * (1 - rv / 3e5), pixel * (1 + rv / 3e5)
            wl_mask = (atmos.wl >= wl_lower) * (atmos.wl < wl_upper)
            assert np.all(atmos.mask[wl_mask] == 1)


@pytest.mark.parametrize("consec_test", [True, False])
def test_atmos_barycenter_shift_mask(sliced_atmmodel_default_mask, consec_test):
    """Test barycentric shift code."""
    # Bary mask should have more pixels mask (at 0) so count will be lower
    atmos = sliced_atmmodel_default_mask
    org_mask = atmos.mask.copy()
    org_number_masked = np.sum(org_mask)
    org_len = len(org_mask)
    atmos.bary_shift_mask(consecutive_test=consec_test)
    new_number_masked = np.sum(atmos.mask)

    assert (new_number_masked < org_number_masked) or (
        (org_len == np.sum(org_number_masked)) or (org_number_masked == 0)
    )
    assert len(atmos.mask) == org_len
    assert np.all(
        (atmos.mask * org_mask) == atmos.mask
    )  # zeros should not turn into ones


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


@pytest.mark.parametrize("consec_test", [True, False])
def test_barycenter_shift_verse_class(short_atmosphere, consec_test):
    """Test barycentric shift code is equivalent inside class.

    TODO: Delete this test when barycenter_shift function is removed."""
    atmos = short_atmosphere
    mask30kms = barycenter_shift(atmos.wl, atmos.mask, consecutive_test=consec_test)

    assert not np.allclose(mask30kms, atmos.mask) or (
        len(mask30kms) == np.sum(mask30kms) or np.sum(mask30kms) == 0
    )
    atmos.bary_shift_mask(consecutive_test=consec_test)
    # They are now close
    assert np.allclose(mask30kms, atmos.mask)


@pytest.mark.xfail()
@pytest.mark.parametrize("resolution", [1000, 50000])
def test_atmos_broadening(atmosphere_fixture, resolution):
    # Test broadening transmission preforms instrumental convolution
    atm = atmosphere_fixture[:4000]
    atm_org = atm.copy()

    #    normalize: bool = True,

    new_trans = resolution_convolution(
        atm_org.wl,
        atm_org.wl,
        atm_org.transmission,
        R=resolution,
        fwhm_lim=5,
        num_procs=1,
        normalize=True,
    )

    atm.broaden(resolution=resolution, num_procs=1)

    assert not np.allclose(
        atm.transmission, atm_org.transmission
    )  # Should have changed
    assert np.allclose(atm.transmission, new_trans)  # Should be equal

    assert np.allclose(atm.wl, atm_org.wl)  # Should not have changed
    assert np.allclose(atm.mask, atm_org.mask)  # Should not have changed
    assert np.allclose(atm.std, atm_org.std)  # Should not have changed


@pytest.mark.xfail()
@pytest.mark.parametrize("resolution", [1000, 20000])
def test_atmos_broadening_reduces_number_of_masked_points(
    atmosphere_fixture, resolution
):
    """Test broadening transmission preforms instrumental convolution."""
    percent = 4
    cuttoff = 1 - percent / 100
    atm = atmosphere_fixture[:4000]
    atm.mask_transmission(percent)

    atm_org = atm.copy()

    atm.broaden(resolution=resolution, num_procs=1)
    assert np.allclose(atm.wl, atm_org.wl)  # Should not have changed
    assert np.allclose(atm.mask, atm_org.mask)  # Should not have changed
    assert np.allclose(atm.std, atm_org.std)  # Should not have changed

    # Transmission changes but number of points in mask may not
    assert not np.allclose(
        atm.transmission, atm_org.transmission
    )  # Should have changed
    atm.mask_transmission(percent)

    if np.allclose(atm.mask, atm_org.mask):  # Mask changes after convolution
        # If masks don't change check both conditions still the same
        assert np.all((atm.transmission > cuttoff) == (atm_org.transmission > cuttoff))


@pytest.mark.xfail()
@pytest.mark.parametrize("resolution", [5000000])
def test_atmos_broadening_at_higher_res_should_not_reduce_number_of_masked_points(
    atmosphere_fixture, resolution
):
    # TODO Fix up this test
    # Test broadening transmission preforms instrumental convolution
    atm = atmosphere_fixture[:5000]
    atm.mask_transmission(1)

    atm_org = atm.copy()

    atm.broaden(resolution=resolution, num_procs=1)
    assert np.allclose(atm.wl, atm_org.wl)  # Should not have changed
    assert np.allclose(atm.mask, atm_org.mask)  # Should not have changed
    assert np.allclose(atm.std, atm_org.std)  # Should not have changed
    assert not np.allclose(
        atm.transmission, atm_org.transmission
    )  # Should have changed
    atm.mask_transmission(1)

    assert np.allclose(atm.mask, atm_org.mask)  # Mask changes after convolution
    assert np.sum(atm.mask) >= np.sum(atm_org.mask)  # Less points are masked out


def test_Atmosphere_has_getitem():
    assert hasattr(Atmosphere, "__getitem__")


def test_Atmosphere_sliceable(short_atmosphere):
    atm = short_atmosphere
    assert len(atm.wl) != 100
    # Indexing on object.
    atm2 = atm[:100]

    assert len(atm2.wl) == 100
    assert len(atm2.transmission) == 100
    assert len(atm2.std) == 100
    assert len(atm2.mask) == 100


def test_Atmosphere_copyable(short_atmosphere):
    atm = short_atmosphere[10:50]

    atm2 = atm.copy()

    assert np.all(atm.wl == atm2.wl)
    assert not (atm.wl is atm2.wl)

    assert np.all(atm.transmission == atm2.transmission)
    assert not (atm.transmission is atm2.transmission)
    assert np.all(atm.std == atm2.std)
    assert not (atm.std is atm2.std)
    assert np.all(atm.mask == atm2.mask)
    assert not (atm.mask is atm2.mask)


def test_Atmosphere_copyable_attr():
    assert hasattr(Atmosphere, "copy")


def test_Atmosphere_band_select():
    """Small test to check band selection."""
    band = "K"  # "K": (2.07, 2.35),

    atm = Atmosphere(
        [2.0, 2.1, 2.2, 2.3, 2.4, 2.5], np.arange(6), std=None, mask=[1, 0, 0, 1, 1, 0]
    )
    atm2 = atm.band_select(band)

    assert np.allclose(atm2.wl, [2.1, 2.2, 2.3])
    assert np.allclose(atm2.mask, [0, 0, 1])


@pytest.mark.parametrize("band", ["TEST", "K"])
def test_from_band_is_a_constructor(band):
    """This is testing loading files using from_band method."""
    atm = Atmosphere.from_band(band)

    min_wl, max_wl = band_limits(band)
    assert isinstance(atm, Atmosphere)
    assert np.all(atm.wl <= max_wl * 1.01)
    assert np.all(atm.wl >= min_wl * 0.99)

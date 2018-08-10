from eniric.precisioncalc import RvPrecision
import numpy as np
import pytest
from eniric.Qcalculator import quality as qual
from astropy.constants import c

@pytest.mark.parametrize("l1, l2, l3", [(2, 2, 2), (60, 60, 60), (30, 30, None)])
def test_even_inputs_allowed(l1, l2, l3):
    """Sizes need to be the same. Initialization of RvPrecision should not fail."""
    if l3 is None:
        l3 = l2

    RvPrecision(np.random.rand(l1), np.random.rand(l2), np.random.rand(l3))
    assert True


@pytest.mark.parametrize("l1, l2, l3", [(1, 2, 3), (20, 60, 20), (5, 5, 1), (1, 10, 10), (1, 2, None)])
def test_uneven_inputs_errors(l1, l2, l3):
    """Sizes need to be the same."""
    if l3 is None:
        l3 = l2

    with pytest.raises(ValueError):
        RvPrecision(np.random.rand(l1), np.random.rand(l2), np.random.rand(l3))


@pytest.mark.parametrize("length", [10, 20, 40])
def test_mask_initialized_as_ones(length):
    """If mask is None."""
    rv = RvPrecision(np.random.randn(length), np.random.randn(length), mask=None)
    assert np.all(rv.mask == 1)
    assert len(rv.mask) == length


# @pytest.fixture()
# def wave_flux_mask():
#     # pull wave and flux from data fixture eventually
#
#     return wave, flux, mask


from eniric.Qcalculator import RVprec_calc, RVprec_calc_masked, RV_prec_calc_Trans


#############################
# Testing consistency these test will be removed later

from hypothesis import given, assume
from hypothesis import strategies as st

#@given(wav=st.lists(st.floats(min_value=0.5, max_value=5000), min_size=2, max_size=500),
#       flux=st.lists(st.floats(min_value=0, max_value=5000), min_size=2, max_size=500),)
@given(st.integers(min_value=2, max_value=100).flatmap(
    lambda n: st.lists(st.lists(st.floats(min_value=0.5, max_value=5000),
                                min_size=n, max_size=n), min_size=2, max_size=2)))
def test_rv1_prec_consistency(wav_flux):
    wavelength = np.asarray(wav_flux[0])
    wavelength.sort()
    assume(np.all(np.diff(wavelength) > 0.000005))
    flux = np.array(wav_flux[0])

    rv_prec = RvPrecision(wavelength, flux)
    old_rv = RVprec_calc(wavelength, flux)

    assert old_rv == rv_prec.rv_prec()

@given(st.integers(min_value=10, max_value=100).flatmap(
    lambda n: st.lists(st.lists(st.floats(min_value=0.5, max_value=5000),
                                min_size=n, max_size=n), min_size=3, max_size=3)))
#@given(wav=st.lists(st.floats(min_value=0.5, max_value=5000), min_size=2, max_size=500),
#       flux=st.lists(st.floats(min_value=0, max_value=5000), min_size=2, max_size=500),
#       mask= st.lists(st.integers(min_value=0, max_value=1), min_size=2, max_size=500))
def test_rv2_prec_clumped_consistency(wav_flux_mask):

    wav = np.array(wav_flux_mask[0])
    wav.sort()
    assume(np.all(np.diff(wav) > 0.000005))
    flux = np.array(wav_flux_mask[1])
    mask = wav_flux_mask[2]
    mask = np.array(mask) > np.mean(mask)
    assume(sum(mask) > len(mask)/4)

    rv_prec = RvPrecision(wav, flux, mask)
    print(mask, np.max(mask), np.min(mask))
    rv2 = rv_prec.rv_prec_clumped()
    rv2_old = RVprec_calc_masked(wav, flux, mask)
    rv2_ = rv_prec.rv_prec_mask()
    assert rv2 == rv2_old
    assert rv2 == rv2_


@given(st.integers(min_value=10, max_value=100).flatmap(
    lambda n: st.lists(st.lists(st.floats(min_value=0.5, max_value=5000),
                                min_size=n, max_size=n), min_size=3, max_size=3)))
def test_rv2_prec_consistency(wav_flux_mask):
    wav = np.array(wav_flux_mask[0])
    wav.sort()
    assume(np.all(np.diff(wav) > 0.000005))
    flux = np.array(wav_flux_mask[1])
    mask = wav_flux_mask[2]
    mask = np.array(mask) > np.mean(mask)

    assume(sum(mask) > len(mask)/4)

    rv_prec = RvPrecision(wav, flux, mask)
    print(mask, np.max(mask), np.min(mask))
    rv2 = rv_prec.rv_prec_mask()
    rv2_old = RVprec_calc_masked(wav, flux, mask)

    assert rv2 == rv2_old



@given(st.integers(min_value=10, max_value=100).flatmap(
    lambda n: st.lists(st.lists(st.floats(min_value=0.5, max_value=5000),
                                min_size=n, max_size=n), min_size=3, max_size=3)))
def test_rv3_prec_consistency(wav_flux_tell):

    wav = np.array(wav_flux_tell[0])
    wav.sort()
    assume(np.all(np.diff(wav) > 0.000005))
    flux = np.array(wav_flux_tell[1])
    tell = wav_flux_tell[2]

    telluric = np.array(tell/np.max(tell))
    assume(np.all(tell >= 0))

    rv_prec = RvPrecision(wav, flux, telluric=telluric)

    rv3 = rv_prec.rv_prec_trans()
    rv3_old = RV_prec_calc_Trans(wav, flux, telluric)
    assert rv3 == rv3_old


def test_pixel_mask_verse_clump():
    length = 500
    wavelength = np.random.randn(length)+1.1
    flux = np.random.randn(length) + 1
    mask = np.random.randint(0, 1, length)
    rv_prec = RvPrecision(wavelength, flux, mask)

    rv2_mask = rv_prec.rv_prec_mask()
    rv2_clump = rv_prec.rv_prec_clumped()
    print(rv2_mask, "masked")
    print(rv2_clump, "clumped")

    assert rv2_mask == rv2_clump


@pytest.mark.parametrize("percentage", [2, 5])
def test_mask_deep_telluric(resampled_data, percentage):
    """Test telluric masking."""
    id, wav, flux = resampled_data
    tell = flux / np.max(flux)

    # Check tell between 0 and 1, and there are some deep telluric lines.
    assert np.all((tell<= 1) | (tell >= 0))
    assert not np.all(tell >= (1 - percentage / 100))

    # Initialize class
    rv = RvPrecision(wav, flux, telluric=tell)

    # Check initial mask
    assert np.all(rv.mask == 1)
    assert not np.all(rv.telluric[rv.mask] >= (1 - percentage/100))

    rv.mask_deep_telluric(percentage)
    assert np.all(rv.telluric[rv.mask] >= (1 - percentage/100))


def test_quality_same(resampled_data):
    id, wav, flux = resampled_data
    tell = flux / np.max(flux)

    # Initialize class
    rv = RvPrecision(wav, flux)

    old_quality = qual(wav, flux)

    assert old_quality == rv.Q()

################################

def test_rv_precision_functional_test(resampled_data):
    # Load a spectra from file
    id_string, wav, flux = resampled_data

    # Obtain the telluric mask spectra
    from eniric_scripts.aces_precision import get_corresponding_atm

    wav_atm, flux_atm, mask_atm = get_corresponding_atm(
        wav, id_string.split("-")[1], bary=True
    )

    # Add to RvPrecision object
    rv_prec = RvPrecision(wav, flux, mask=mask_atm, telluric=flux_atm)

    # Calculate spectral quality
    quality = rv_prec.Q()
    print("quality", quality)
    print("quality old", qual(wav, flux))

    # Calculate the precision for full spectrum
    rv1 = rv_prec.rv_prec()
    print("rv1", rv1)
    # assert general relation between rv1 and quality
    assert rv1 == c / (quality * np.sqrt(rv_prec.total_flux()))

    # rv3
    rv3 = rv_prec.rv_prec_trans()
    rv3_old = RV_prec_calc_Trans(wav, flux, flux_atm)
    print("rv3", rv3)
    print("rv3_old", rv3_old)

    # Calculate the precision from masking
    rv2 = rv_prec.rv_prec_mask()
    print("rv2", rv2)
    print("rv2 old", RVprec_calc_masked(wav, flux, mask=mask_atm))

    # Check that it is the same as clumped masking
    rv2_clumped = rv_prec.rv_prec_clumped()

    assert (rv1 < rv3) and (rv3 < rv2)
    assert rv2 == rv2_clumped


#@pytest.fixture()
#TODO multiple tests hypothesis?

# TODo sort out nan/inf issues.
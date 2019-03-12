import numpy as np
import pytest
from astropy import units as u
from hypothesis import given, strategies as st

from eniric.atmosphere import Atmosphere
from eniric.legacy import RVprec_calc_masked, RVprec_calc_weights_masked, mask_clumping


def test_RV_prec_masked(test_spec):
    """Test same precision results between past pre-clumped version and mask clump version."""
    wav = test_spec[0]
    flux = test_spec[1]
    mask = test_spec[2]
    if mask is not None:
        mask = mask.round().astype(bool)
    else:
        mask = np.ones_like(wav)
    print(mask)

    # Pre clumping as done in Figueira et al. 2016
    wav_masked, flux_masked = mask_clumping(wav, flux, mask)
    rv_chunks = RVprec_calc_masked(wav_masked, flux_masked, mask=None)

    rv_masked = RVprec_calc_masked(wav, flux, mask)

    assert rv_masked.value == rv_chunks.value
    assert isinstance(rv_masked, u.Quantity)
    assert rv_masked.unit == u.m / u.s


@given(st.lists(st.booleans(), min_size=5, max_size=300))
def test_mask_clumping_of_mask(mask):
    """Masking mask show return all ones."""
    wav_clumped, flux_clumped = mask_clumping(mask, mask, mask)
    assert len(wav_clumped) == len(flux_clumped)
    for wav_i, flux_i in zip(wav_clumped, flux_clumped):
        assert np.all(wav_i)
        assert np.all(flux_i)
    # sum of masked_clumped should equal mask sum
    mask_sum = np.sum(mask)
    assert np.sum([np.sum(wav_i) for wav_i in wav_clumped]) == mask_sum
    assert np.sum([np.sum(flux_i) for flux_i in flux_clumped]) == mask_sum


def test_manual_clumping():
    """Test properties of clumping function using manually defined masked_arrays."""
    wav = np.arange(15)
    flux = np.arange(15, 30)
    mask = np.array([1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0])
    mask_bool = np.array([1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0], dtype=bool)

    wav_masked, flux_masked = mask_clumping(wav, flux, mask)
    wav_masked_bool, flux_masked_bool = mask_clumping(wav, flux, mask_bool)

    expected_wav = [np.arange(0, 4), np.arange(7, 10), np.arange(11, 14)]
    expected_flux = [np.arange(15, 19), np.arange(22, 25), np.arange(26, 29)]
    for i, __ in enumerate(wav_masked):
        assert np.allclose(wav_masked[i], expected_wav[i])
        assert np.allclose(flux_masked[i], expected_flux[i])
        assert np.allclose(wav_masked_bool[i], expected_wav[i])
        assert np.allclose(flux_masked_bool[i], expected_flux[i])

    assert len(expected_wav) == len(wav_masked)
    assert len(expected_flux) == len(flux_masked)


def test_legacy_RV_warns_nonfinite(grad_flag):
    """Some warning tests."""
    with pytest.warns(RuntimeWarning, match="divide by zero"):
        RVprec_calc_masked(
            np.array([1, 2, 3, 4]),
            np.array([1, 2, 3, 4]),
            np.array([0, 1, 0, 0]),
            grad=grad_flag,
        )


@pytest.mark.xfail(ModuleNotFoundError, reason="Issue with Starfish install.")
@pytest.mark.parametrize("band", ["H", "J", "K"])
def test_weights_clumping_grad(testing_spectrum, grad_flag, band):
    # Test masked clumping verse weight mask with new gradients
    # Test on an actual spectrum to check sizes of difference
    wav, flux = testing_spectrum
    atm = Atmosphere.from_band(band)
    mask = np.ones_like(wav)
    mask = atm.at(wav).mask

    clumped = RVprec_calc_masked(wav, flux, mask=mask, grad=grad_flag)
    weighted = RVprec_calc_weights_masked(wav, flux, mask=mask, grad=grad_flag)

    # They are not exactly the same by are within a specific percentage.
    ratio = (weighted.value - clumped.value) / clumped.value

    if grad_flag:
        # The difference for the new gradient is smaller.
        assert abs(ratio) < 0.008
    else:
        assert abs(ratio) < 0.04
    assert clumped.unit == weighted.unit

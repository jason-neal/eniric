"""Test Qcaclulator."""

import astropy.units as u
import numpy as np
import pytest
from astropy import constants as const

import eniric.Qcalculator as Q

# test RVprec_calc returns a single values
# test it returns a quantity in m/s
# test it can handle wavelength and flux alo being quantities.

m_per_s = u.meter / u.second
per_s_cm2 = (1 / u.second) / (u.centimeter ** 2)


def test_rvprev_calc():
    """Test that RVprec_calc can handle inputs as Quantities or unitless and returns a Quantity."""
    wav = np.arange(1, 101)
    flux = np.random.random(100)

    rv = Q.RVprec_calc(wav, flux)
    assert rv.unit == m_per_s
    assert not hasattr(rv.value, "__len__")  # assert value is a scalar
    assert isinstance(rv, u.Quantity)

    # Test it can handle quantities inputs
    rv1 = Q.RVprec_calc(wav * u.micron, flux)
    assert not hasattr(rv1.value, "__len__")  # assert value is a scalar
    assert rv1 == rv
    assert rv1.unit == m_per_s

    rv2 = Q.RVprec_calc(wav * u.micron, flux * per_s_cm2)
    assert not hasattr(rv2.value, "__len__")  # assert value is a scalar
    assert rv1 == rv2
    assert rv == rv2
    assert rv2.unit == m_per_s


def test_rvprev_calc_with_lists():
    """Test that it can handle list input also."""
    wav = list(np.arange(100))
    flux = list(np.random.random(100))

    rv = Q.RVprec_calc(wav, flux)
    assert not hasattr(rv.value, "__len__")  # assert value is a scalar
    assert isinstance(rv, u.Quantity)
    assert rv.unit == m_per_s


def test_sqrt_sum_wis():
    """Test that sqrt_sum_wis can handle inputs as Quantities or unitless.

    Returns a dimensionless unscaled Quantity.
    """
    wav = np.arange(1, 101)
    flux = np.random.random(100)

    sqrtsumwis = Q.sqrt_sum_wis(wav, flux)
    assert not isinstance(
        sqrtsumwis, u.Quantity
    )  # Doesn't turn into quantity if does not have to.
    assert not hasattr(sqrtsumwis, "__len__")  # assert value is a scalar

    sqrtsumwis2 = Q.sqrt_sum_wis(
        wav * u.micron, (flux / u.second) / (u.centimeter ** 2)
    )  # with some units
    assert not hasattr(sqrtsumwis2.value, "__len__")  # assert value is a scalar
    assert isinstance(sqrtsumwis2, u.Quantity)
    assert (
        sqrtsumwis2.unit == u.dimensionless_unscaled
    )  # unscaled and dimensionless quantity

    assert sqrtsumwis == sqrtsumwis2.value

    # Test relation to RVprec_calc
    assert Q.RVprec_calc(wav, flux) == const.c / sqrtsumwis


def test_RV_prec_calc_Trans():
    """Transmission should not have units."""
    wav = np.arange(1, 101)
    flux = np.random.random(100)
    trans = np.random.random(100)

    rv_trans = Q.RV_prec_calc_Trans(wav, flux, trans)
    assert not hasattr(rv_trans.value, "__len__")  # assert scalar
    assert rv_trans.unit == m_per_s

    # dimensionless_unscaled unit is ok
    rv_trans2 = Q.RV_prec_calc_Trans(wav, flux, trans * u.dimensionless_unscaled)
    assert not hasattr(rv_trans2.value, "__len__")  # assert  scalar
    assert rv_trans2.unit == m_per_s

    assert rv_trans == rv_trans2

    with pytest.raises(TypeError):
        # transmission mistakenly given as a flux unit
        Q.RV_prec_calc_Trans(wav, flux, (trans / u.s) / (u.centimeter ** 2))

    with pytest.raises(ValueError):
        Q.RV_prec_calc_Trans(wav, flux, trans + 1)

    with pytest.raises(ValueError):
        Q.RV_prec_calc_Trans(wav, flux, trans * -5)


def test_SQRTSumWisTrans():
    """Test square root sum of weights when including change of variance due to atmospheric transmission."""
    wav = np.arange(1, 101)
    flux = np.random.random(100)
    trans = np.random.random(100)

    swrtsum_trans = Q.sqrt_sum_wis_trans(wav, flux, trans)
    assert not isinstance(
        swrtsum_trans, u.Quantity
    )  # Doesn't turn into quantity if does not have to.
    assert not hasattr(swrtsum_trans, "__len__")  # assert scalar

    # dimensionless_unscaled unit is ok for transmission
    sqrtsum_trans2 = Q.sqrt_sum_wis_trans(wav, flux, trans * u.dimensionless_unscaled)
    assert not hasattr(sqrtsum_trans2.value, "__len__")  # assert scalar
    assert isinstance(sqrtsum_trans2, u.Quantity)
    assert (
        sqrtsum_trans2.unit == u.dimensionless_unscaled
    )  # unscaled and dimensionless quantity

    sqrtsum_trans3 = Q.sqrt_sum_wis_trans(wav * u.micron, flux, trans)
    assert not hasattr(sqrtsum_trans3.value, "__len__")  # assert value is a scalar
    assert isinstance(sqrtsum_trans3, u.Quantity)
    assert (
        sqrtsum_trans3.unit == u.dimensionless_unscaled
    )  # unscaled and dimensionless quantity

    with pytest.raises(TypeError):
        # transmission mistakenly given as a flux unit
        Q.sqrt_sum_wis_trans(wav, flux, (trans / u.s) / (u.centimeter ** 2))

    with pytest.raises(ValueError):
        Q.sqrt_sum_wis_trans(wav, flux, trans + 1)

    with pytest.raises(ValueError):
        Q.sqrt_sum_wis_trans(wav, flux, trans * -5)


def test_transmission_reduces_precision():
    """Check that a transmission vector reduces precision calculation."""
    wav = np.arange(100.)
    flux = np.random.random(100)
    transmission = np.random.random(100)

    # Value should be less then normal if trans <=1
    assert Q.RVprec_calc(wav, flux) < Q.RV_prec_calc_Trans(wav, flux, transmission)
    # Unitary transmission should give equivalent result.
    assert Q.RVprec_calc(wav, flux) == Q.RV_prec_calc_Trans(
        wav, flux, np.ones_like(wav)
    )


def test_RV_prec_masked():
    """Test same prections results between past pre-clumped version and mask version."""
    wav = np.arange(100)
    flux = np.random.random(100) * 10
    mask = np.asarray(np.floor(2 * np.random.random(100)), dtype=bool)

    # Pre clumping as in nIR_precion.py
    wav_chunks, flux_chunks = Q.bug_fixed_clumping_method(wav, flux, mask)
    rv_chunks = Q.RVprec_calc_masked(wav_chunks, flux_chunks, mask=None)

    rv_masked = Q.RVprec_calc_masked(wav, flux, mask)

    assert rv_masked.value == rv_chunks.value
    assert isinstance(rv_masked, u.Quantity)
    assert rv_masked.unit == u.m / u.s


def test_mask_clumping():
    """Test properties of clumping function using masked_arrays."""
    wav = np.arange(1, 101)
    flux = np.random.random(100) * 10
    mask = np.asarray(np.floor(2 * np.random.random(100)), dtype=bool)

    wav_chunks, flux_chunks = Q.bug_fixed_clumping_method(wav, flux, mask)
    wav_masked, flux_masked = Q.mask_clumping(wav, flux, mask)
    for i, __ in enumerate(wav_chunks):
        assert np.all(wav_chunks[i] == wav_masked[i])
        assert np.all(flux_chunks[i] == flux_masked[i])

    # All values of clumped mask should be True.
    mask_clumped, mask2_clumped = Q.mask_clumping(mask, mask, mask)
    assert len(mask_clumped) == len(mask2_clumped)
    for i, __ in enumerate(mask_clumped):
        assert np.all(mask_clumped[i])
        assert np.all(mask2_clumped[i])


def test_manual_clumping():
    """Test properties of clumping function using masked_arrays."""
    wav = np.arange(15)
    flux = np.arange(15, 30)
    mask = np.array([1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0])  # , dtype=bool
    mask_bool = np.array([1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0], dtype=bool)

    wav_chunks, flux_chunks = Q.bug_fixed_clumping_method(wav, flux, mask)
    wav_masked, flux_masked = Q.mask_clumping(wav, flux, mask)
    wav_masked_bool, flux_masked_bool = Q.mask_clumping(wav, flux, mask_bool)

    expected_wav = [np.arange(0, 4), np.arange(7, 10), np.arange(11, 14)]
    expected_flux = [np.arange(15, 19), np.arange(22, 25), np.arange(26, 29)]
    for i, chunk in enumerate(wav_chunks):
        assert np.all(wav_chunks[i] == expected_wav[i])
        assert np.all(flux_chunks[i] == expected_flux[i])
        assert np.all(wav_masked[i] == expected_wav[i])
        assert np.all(flux_masked[i] == expected_flux[i])
        assert np.all(wav_masked_bool[i] == expected_wav[i])
        assert np.all(flux_masked_bool[i] == expected_flux[i])

    assert len(wav_chunks) == len(wav_masked)
    assert len(flux_chunks) == len(flux_masked)
    assert len(expected_wav) == len(wav_masked)
    assert len(expected_flux) == len(flux_masked)

    # All values of clumped mask should be True.
    mask_clumped, mask2_clumped = Q.mask_clumping(mask, mask, mask)
    assert len(mask_clumped) == len(mask2_clumped)
    for i, __ in enumerate(mask_clumped):
        assert np.all(mask_clumped[i])
        assert np.all(mask2_clumped[i])


@pytest.mark.parametrize("wave_unit", [1, u.centimeter, u.nanometer])
@pytest.mark.parametrize("flux_unit", [1, per_s_cm2, 1. / u.second])
def test_sqrt_sum_wis_with_quantities(wave_unit, flux_unit):
    """Assert that wis returns dimensionless.

    Assert that wis returns with quantities is ok.
    """
    wav = np.arange(1, 101) * wave_unit
    flux = (np.random.randn(100) + 1) * flux_unit
    wis = Q.sqrt_sum_wis(wav, flux)

    if isinstance(wis, u.Quantity):
        assert wis.unit == u.dimensionless_unscaled
    else:
        assert True


@pytest.mark.parametrize("wave_unit", [1, u.centimeter, u.nanometer])
@pytest.mark.parametrize("flux_unit", [1, per_s_cm2, 1. / u.second])
def test_sqrt_sum_wis_trans_with_quantities(wave_unit, flux_unit):
    """Assert that wis returns dimensionless.

    Assert that wis returns with quantities is ok.
    """
    wav = np.arange(1, 101) * wave_unit
    flux = (np.random.randn(100) + 1) * flux_unit
    transmission = np.random.rand(len(wav))
    wis = Q.sqrt_sum_wis_trans(wav, flux, transmission)

    if isinstance(wis, u.Quantity):
        assert wis.unit == u.dimensionless_unscaled
    else:
        assert True


@pytest.mark.parametrize("wave_unit", [1, u.nanometer])
@pytest.mark.parametrize("flux_unit", [1, per_s_cm2])
@pytest.mark.parametrize("trans_unit", [m_per_s, per_s_cm2, u.meter])
def test_sqrt_sum_wis_trans_with_trans_unit_fails(wave_unit, flux_unit, trans_unit):
    """Assert a transmission with a unit fails with type error."""
    wav = np.arange(1, 101) * wave_unit
    flux = (np.random.randn(100) + 1) * flux_unit
    transmission = np.random.rand(len(wav))

    with pytest.raises(TypeError):
        Q.sqrt_sum_wis_trans(wav, flux, transmission * trans_unit)


@pytest.mark.parametrize("wave_unit", [1, u.centimeter, u.nanometer])
@pytest.mark.parametrize("flux_unit", [1, per_s_cm2, 1. / u.second])
def test_sqrt_sum_wis_transmission_outofbounds(wave_unit, flux_unit):
    """Transmission must be within 0-1."""
    wav = np.arange(1, 101) * wave_unit
    flux = (np.random.randn(100) + 1) * flux_unit
    transmission1 = np.random.randn(len(wav))
    transmission2 = np.random.rand(len(wav))

    transmission1[0] = 5  # Outside 0-1
    transmission2[-1] = -2  # Outside 0-1

    with pytest.raises(ValueError):
        Q.sqrt_sum_wis_trans(wav, flux, transmission1)  # Higher value
    with pytest.raises(ValueError):
        Q.sqrt_sum_wis_trans(wav, flux, transmission2)  # Lower value


@pytest.mark.parametrize("scale", [0.1, 1, 2, 100, 0.1, 0.5])
def test_quality_independent_of_flux_level(scale):
    """Q of a spectrum is independent of flux level."""
    wavelength = np.arange(100)
    flux = np.random.random(100)
    assert np.allclose(Q.quality(wavelength, flux), Q.quality(wavelength, flux * scale))


@pytest.mark.parametrize("wave_unit", [1, u.centimeter, u.nanometer])
@pytest.mark.parametrize("flux_unit", [1, 1 / u.second / u.micron ** 2, 1. / u.second])
def test_quality_independent_of_units(wave_unit, flux_unit):
    """Quality should be unitless, or dimensionless_unscaled...

    Needed to remove unit from flux in quality() to make quality unitless.
    """
    wave = np.arange(10) + 1
    flux = np.random.randn(10) + 1
    wave = wave * wave_unit
    flux = flux * flux_unit
    q = Q.quality(wave, flux)
    print(q)
    if isinstance(q, u.Quantity):
        assert q.unit == u.dimensionless_unscaled
    else:
        assert True

"""Test Qcalculator."""

import astropy.units as u
import hypothesis.strategies as st
import numpy as np
import pytest
from astropy import constants as const
from astropy.units import Quantity
from hypothesis import given

import eniric.Qcalculator as Q
from eniric.Qcalculator import mask_check, pixel_weights

m_per_s = u.meter / u.second
per_s_cm2 = (1 / u.second) / (u.centimeter ** 2)
c = const.c


# Define some fixtures for Qcalculator.
@pytest.fixture(
    params=[
        (np.arange(1, 101), np.random.random(100), None),
        (np.linspace(2.1, 2.5, 200), np.random.random(200), np.random.random(200)),
        (
            np.linspace(0.5, 1.5, 50),
            np.random.random(50),
            np.floor(2 * np.random.random(50)),
        ),
    ]
)
def test_spec(request):
    """Wave and flux, mask examples."""
    return request.param


# 3 situations each variable, no unit. a unit. or a dimensionless unscaled unit.
@pytest.fixture(params=[1, u.micron, u.dimensionless_unscaled])
def wav_unit(request):
    """Iterate some units on wavelength."""
    return request.param


@pytest.fixture(params=[1, per_s_cm2, u.dimensionless_unscaled])
def flux_unit(request):
    """Iterate some units on flux"""
    return request.param


@pytest.fixture(params=[1, u.dimensionless_unscaled])
def trans_unit(request):
    """Iterate some units on mask/transmission."""
    return request.param


@pytest.fixture(params=[True, False])
def grad_flag(request):
    """Gradient flag parameter."""
    return request.param


def test_rvprev_calc(test_spec, wav_unit, flux_unit, trans_unit):
    """Test that RVprec_calc can handle inputs as Quantities or unitless and returns a scalar Quantity."""
    wav = test_spec[0] * wav_unit
    flux = test_spec[1] * flux_unit
    mask = test_spec[2]
    if test_spec[2] is not None:
        mask *= trans_unit

    rv = Q.RVprec_calc(wav, flux, mask)
    assert rv.unit == m_per_s
    assert not hasattr(rv.value, "__len__")  # assert value is a scalar
    assert isinstance(rv, u.Quantity)


def test_rvprev_calc_with_lists(test_spec):
    """Test that it can handle list input also."""
    wav = list(test_spec[0])
    flux = list(test_spec[1])
    mask = test_spec[2]
    rv = Q.RVprec_calc(wav, flux, mask)
    assert not hasattr(rv.value, "__len__")  # assert value is a scalar
    assert isinstance(rv, u.Quantity)
    assert rv.unit == m_per_s


def test_sqrt_sum_wis_with_no_units(test_spec):
    """Test that sqrt_sum_wis can handle inputs as Quantities or unitless.
    Returns a dimensionless unscaled Quantity.
    """
    sqrtsumwis = Q.sqrt_sum_wis(test_spec[0], test_spec[1], test_spec[2])
    # Doesn't turn into quantity if does not have to.
    assert not isinstance(sqrtsumwis, u.Quantity)
    assert not hasattr(sqrtsumwis, "__len__")  # assert value is a scalar


def test_sqrt_sum_wis(test_spec, wav_unit, flux_unit, trans_unit):
    """Test that sqrt_sum_wis can handle inputs as Quantities or unitless.

    Returns a dimensionless unscaled Quantity.
    """
    wav = test_spec[0] * wav_unit
    flux = test_spec[1] * flux_unit
    mask = test_spec[2]
    if test_spec[2] is not None:
        mask *= trans_unit

    sqrtsumwis = Q.sqrt_sum_wis(wav, flux, mask)
    print("wav", wav, type(wav))
    print("flux", flux, type(flux))
    print("mask", mask, type(mask))
    print("sqrtsumwis", sqrtsumwis, type(sqrtsumwis))
    if (
        (isinstance(wav, Quantity))
        or (isinstance(flux, Quantity))
        or (isinstance(mask, Quantity) and (test_spec[2] is not None))
    ):
        assert isinstance(sqrtsumwis, u.Quantity)
        # Check is unscaled and dimensionless Quantity
        assert sqrtsumwis.unit == u.dimensionless_unscaled
        assert not hasattr(sqrtsumwis.value, "__len__")  # assert value is a scalar
    else:
        assert not hasattr(sqrtsumwis, "__len__")  # assert value is a scalar


def test_relation_of_rv_to_sqrtsumwis(test_spec, wav_unit, flux_unit, trans_unit):
    """Test relation of sqrtsumwis to RVprec_calc."""
    wav = test_spec[0] * wav_unit
    flux = test_spec[1] * flux_unit
    mask = test_spec[2]
    if test_spec[2] is not None:
        mask *= trans_unit
    assert np.all(
        Q.RVprec_calc(wav, flux, mask=mask) == c / Q.sqrt_sum_wis(wav, flux, mask=mask)
    )


def test_transmission_reduces_precision(test_spec):
    """Check that a transmission vector reduces precision calculation."""
    wav = test_spec[0]
    flux = test_spec[1]
    transmission = test_spec[2]

    # Value should be less then normal if trans <=1
    if transmission is not None:
        assert Q.RVprec_calc(wav, flux, mask=None) < Q.RVprec_calc(
            wav, flux, mask=transmission
        )
    # mask=None is the same as mask of all 1.
    assert Q.RVprec_calc(wav, flux, mask=None) == Q.RVprec_calc(
        wav, flux, mask=np.ones_like(wav)
    )


def test_improved_gradient_reduces_precision(test_spec):
    """Check that the gradient produces  larger RV error."""
    wav = test_spec[0]
    flux = test_spec[1]
    transmission = test_spec[2]

    # mask=None is the same as mask of all 1.
    assert Q.RVprec_calc(wav, flux, mask=transmission, grad=False) <= Q.RVprec_calc(
        wav, flux, mask=transmission, grad=True
    )


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

    # Pre clumping as in nIR_precision.py
    wav_masked, flux_masked = Q.mask_clumping(wav, flux, mask)
    rv_chunks = Q.RVprec_calc_masked(wav_masked, flux_masked, mask=None)

    rv_masked = Q.RVprec_calc_masked(wav, flux, mask)

    assert rv_masked.value == rv_chunks.value
    assert isinstance(rv_masked, u.Quantity)
    assert rv_masked.unit == u.m / u.s


@given(st.lists(st.booleans(), min_size=5, max_size=300))
def test_mask_clumping_of_mask(mask):
    """Masking mask show return all ones."""
    wav_clumped, flux_clumped = Q.mask_clumping(mask, mask, mask)
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

    wav_masked, flux_masked = Q.mask_clumping(wav, flux, mask)
    wav_masked_bool, flux_masked_bool = Q.mask_clumping(wav, flux, mask_bool)

    expected_wav = [np.arange(0, 4), np.arange(7, 10), np.arange(11, 14)]
    expected_flux = [np.arange(15, 19), np.arange(22, 25), np.arange(26, 29)]
    for i in range(len(wav_masked)):
        assert np.allclose(wav_masked[i], expected_wav[i])
        assert np.allclose(flux_masked[i], expected_flux[i])
        assert np.allclose(wav_masked_bool[i], expected_wav[i])
        assert np.allclose(flux_masked_bool[i], expected_flux[i])

    assert len(expected_wav) == len(wav_masked)
    assert len(expected_flux) == len(flux_masked)


@pytest.mark.parametrize("scale", [0.1, 1, 2, 100, 0.1, 0.5])
def test_quality_independent_of_flux_level(scale):
    """Q of a spectrum is independent of flux level."""
    wavelength = np.arange(100)
    flux = np.random.random(100)
    assert np.allclose(Q.quality(wavelength, flux), Q.quality(wavelength, flux * scale))


def test_quality_independent_of_units(test_spec, wav_unit, flux_unit):
    """Quality should be returned as unitless (not a quantity)."""
    wave = test_spec[0] * wav_unit
    flux = test_spec[1] * flux_unit
    q = Q.quality(wave, flux)

    assert not isinstance(q, Quantity)
    assert isinstance(q, float)
    assert not hasattr(q, "__len__")  # assert value is a scalar


@pytest.mark.parametrize(
    "wav_unit2, flux_unit2, quantity",
    ([(u.nanometer, 1 / u.second, True), (u.meter, u.erg, True), (1, 1, False)]),
)
@pytest.mark.parametrize(
    "wave, flux, expected",
    [
        ([1, 2, 3], [1, 2, 3], [1, 2, 3]),
        ([2.2, 2.3, 2.4], [.99, 0.97, 0.7], [0.195555556, 11.46621134, 59.9868571]),
    ],
)
def test_pixel_weights_gradient_with_fixed_values(
    wave, flux, expected, wav_unit2, flux_unit2, quantity
):
    """Pixel_weights with gradient formulation.
    Change detector with some simple numbers.
    """
    result = pixel_weights(wave * wav_unit2, flux * flux_unit2, grad=True)
    assert np.allclose(result, expected)
    # Checks quantities also work
    if quantity:
        assert isinstance(result, Quantity)
        assert result.unit == u.dimensionless_unscaled


@pytest.mark.parametrize(
    "wav_unit2, flux_unit2, quantity",
    [(u.nanometer, 1 / u.second, True), (u.meter, u.erg, True), (1, 1, False)],
)
@pytest.mark.parametrize(
    "wave, flux, expected",
    [
        ([1, 2, 3], [1, 2, 3], [1, 2]),
        ([2.2, 2.3, 2.4], [.99, 0.97, 0.7], [0.195555556, 39.75680412]),
    ],
)
def test_pixel_weights(wave, flux, expected, wav_unit2, flux_unit2, quantity):
    """Pixel_weights with finite difference formulation.
    Change detector with some simple numbers.
    """
    print(wave, flux, expected, wav_unit2, flux_unit2, quantity)
    result = pixel_weights(wave * wav_unit2, flux * flux_unit2, grad=False)
    assert np.allclose(result, expected)
    # Checks quantities also work
    if quantity:
        assert isinstance(result, Quantity)
        assert result.unit == u.dimensionless_unscaled


@pytest.mark.parametrize(
    "mask",
    [np.array([1, 0, -0.1]), np.array([1.1, 0, 1])],  # lower than 0  # Greater than 1
)
def test_mask_check(mask):
    with pytest.raises(ValueError):
        mask_check(mask)


@pytest.mark.parametrize(
    "mask", [np.array([1, 0, 1, 0, .2, 1, .4]) * u.m, np.array([1, 0, 1]) * u.s]
)
def test_mask_check_type_error(mask):
    with pytest.raises(TypeError):
        mask_check(mask)


@pytest.mark.parametrize("trans_unit2", [m_per_s, per_s_cm2, u.meter])
def test_sqrt_sum_wis_with_mask_with_unit_fails(
    test_spec, wav_unit, flux_unit, trans_unit2
):
    """Assert a transmission with a unit fails with type error."""
    wav = test_spec[0] * wav_unit
    flux = test_spec[1] * flux_unit
    transmission = np.random.rand(len(wav)) * trans_unit2

    with pytest.raises(TypeError):
        Q.sqrt_sum_wis(wav, flux, mask=transmission)

    with pytest.raises(TypeError):
        Q.RVprec_calc(wav, flux, mask=transmission)


def test_sqrt_sum_wis_transmission_outofbounds(test_spec, wav_unit, flux_unit):
    """Transmission must be within 0-1."""
    wav = test_spec[0] * wav_unit
    flux = test_spec[1] * flux_unit
    transmission1 = np.random.randn(len(wav))
    transmission2 = np.random.rand(len(wav))

    transmission1[0] = 5  # Outside 0-1
    transmission2[-1] = -2  # Outside 0-1

    # Higher value
    with pytest.raises(ValueError):
        Q.RVprec_calc(wav, flux, mask=transmission1)

    with pytest.raises(ValueError):
        Q.sqrt_sum_wis(wav, flux, mask=transmission1)

        # Lower value
    with pytest.raises(ValueError):
        Q.sqrt_sum_wis(wav, flux, mask=transmission2)

    with pytest.raises(ValueError):
        Q.sqrt_sum_wis(wav, flux, mask=transmission2)


def test_sqrtsumwis_warns_nonfinite(grad_flag):
    """Some warning tests."""
    with pytest.warns(RuntimeWarning, match="divide by zero"):
        Q.RVprec_calc_masked(
            np.array([1, 2, 3, 4]),
            np.array([1, 2, 3, 4]),
            np.array([0, 1, 0, 0]),
            grad=grad_flag,
        )

    with pytest.warns(RuntimeWarning, match="divide by zero"):
        Q.RVprec_calc_masked(
            np.array([1, 2, 3, 4]),
            np.array([1, 2, 3, 4]),
            np.array([0, 1, 0, 0]),
            grad=grad_flag,
        )

    with pytest.warns(UserWarning, match="This will cause infinite errors."):
        Q.sqrt_sum_wis(
            np.array([1, 2, 3]),
            np.array([1, 2, 3]),
            np.array([0, 0, 0]),
            grad=grad_flag,
        )  # All masked

    with pytest.warns(UserWarning, match="Weight sum is not finite"):
        Q.sqrt_sum_wis(
            np.array([2, 2, 2]),
            np.array([1, 2, 3]),
            np.array([1, 1, 1]),
            grad=grad_flag,
        )  # infinate gradient

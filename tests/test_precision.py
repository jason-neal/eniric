"""Test eniric.precision."""

import astropy.units as u
import numpy as np
import pytest
from astropy import constants as const
from astropy.units import Quantity

from eniric.atmosphere import Atmosphere
from eniric.precision import (
    incremental_quality,
    incremental_rv,
    mask_check,
    pixel_weights,
    quality,
    rv_precision,
    sqrt_sum_wis,
)
from eniric.snr_normalization import snr_constant_band
from eniric.utilities import (
    band_limits,
    load_aces_spectrum,
    wav_selector,
    weighted_error,
)

m_per_s = u.meter / u.second
per_s_cm2 = (1 / u.second) / (u.centimeter ** 2)
c = const.c

xfail = pytest.mark.xfail


def test_rvprev_calc(test_spec, wav_unit, flux_unit, trans_unit):
    """Test that rv_precision can handle inputs as Quantities or unitless and returns a scalar Quantity."""
    wav = test_spec[0] * wav_unit
    flux = test_spec[1] * flux_unit
    mask = test_spec[2]
    if test_spec[2] is not None:
        mask *= trans_unit

    rv = rv_precision(wav, flux, mask)
    assert rv.unit == m_per_s
    assert not hasattr(rv.value, "__len__")  # assert value is a scalar
    assert isinstance(rv, u.Quantity)


def test_rvprev_calc_with_lists(test_spec):
    """Test that it can handle list input also."""
    wav = list(test_spec[0])
    flux = list(test_spec[1])
    mask = test_spec[2]
    rv = rv_precision(wav, flux, mask)
    assert not hasattr(rv.value, "__len__")  # assert value is a scalar
    assert isinstance(rv, u.Quantity)
    assert rv.unit == m_per_s


def test_sqrt_sum_wis_with_no_units(test_spec):
    """Test that sqrt_sum_wis can handle inputs as Quantities or unitless.
    Returns a dimensionless unscaled Quantity.
    """
    sqrtsumwis = sqrt_sum_wis(test_spec[0], test_spec[1], test_spec[2])
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

    sqrtsumwis = sqrt_sum_wis(wav, flux, mask)
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
    """Test relation of sqrtsumwis to rv_precision."""
    wav = test_spec[0] * wav_unit
    flux = test_spec[1] * flux_unit
    mask = test_spec[2]
    if test_spec[2] is not None:
        mask *= trans_unit
        mask = mask ** 2
    assert np.all(
        rv_precision(wav, flux, mask=mask) == c / sqrt_sum_wis(wav, flux, mask=mask)
    )


def test_transmission_reduces_precision(test_spec):
    """Check that a transmission vector reduces precision calculation."""
    wav = test_spec[0]
    flux = test_spec[1]
    transmission = test_spec[2]

    # Value should be less then normal if trans <=1
    if transmission is not None:
        assert rv_precision(wav, flux, mask=None) < rv_precision(
            wav, flux, mask=transmission
        )
    # mask=None is the same as mask of all 1.
    assert rv_precision(wav, flux, mask=None) == rv_precision(
        wav, flux, mask=np.ones_like(wav)
    )


def test_both_gradients_possible(test_spec):
    wav = test_spec[0]
    flux = test_spec[1]
    transmission = test_spec[2]
    __ = rv_precision(wav, flux, mask=transmission, grad=False).value
    __ = rv_precision(wav, flux, mask=transmission, grad=True).value
    assert True


@pytest.mark.parametrize("scale", [0.1, 1, 2, 100, 0.1, 0.5])
def test_quality_independent_of_flux_level(scale):
    """Q of a spectrum is independent of flux level."""
    wavelength = np.arange(100)
    flux = np.random.random(100)
    assert np.allclose(quality(wavelength, flux), quality(wavelength, flux * scale))


def test_quality_independent_of_units(test_spec, wav_unit, flux_unit):
    """Quality should be returned as unitless (not a quantity)."""
    wave = test_spec[0] * wav_unit
    flux = test_spec[1] * flux_unit
    q = quality(wave, flux)

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
        (
            [2.2, 2.3, 2.4],
            [0.99, 0.97, 0.7],
            [0.195_555_556, 11.466_211_34, 59.986_857_1],
        ),
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
        ([2.2, 2.3, 2.4], [0.99, 0.97, 0.7], [0.195_555_556, 39.756_804_12]),
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
    "mask", [np.array([1, 0, 1, 0, 0.2, 1, 0.4]) * u.m, np.array([1, 0, 1]) * u.s]
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
        sqrt_sum_wis(wav, flux, mask=transmission ** 2)

    with pytest.raises(TypeError):
        rv_precision(wav, flux, mask=transmission ** 2)


def test_sqrt_sum_wis_transmission_outofbounds(test_spec, wav_unit, flux_unit):
    """Transmission must be within 0-1."""
    wav = test_spec[0] * wav_unit
    flux = test_spec[1] * flux_unit
    mask_1 = np.random.randn(len(wav))
    mask_2 = np.random.rand(len(wav))

    mask_1[0] = 5  # Outside 0-1
    mask_2[-1] = -2  # Outside 0-1

    # Higher value
    with pytest.raises(ValueError):
        rv_precision(wav, flux, mask=mask_1)

    with pytest.raises(ValueError):
        sqrt_sum_wis(wav, flux, mask=mask_1)

        # Lower value
    with pytest.raises(ValueError):
        sqrt_sum_wis(wav, flux, mask=mask_2)

    with pytest.raises(ValueError):
        sqrt_sum_wis(wav, flux, mask=mask_2)


def test_sqrtsumwis_warns_nonfinite(grad_flag):
    """Some warning tests."""
    with pytest.warns(UserWarning, match="This will cause infinite errors."):
        sqrt_sum_wis(
            np.array([1, 2, 3]),
            np.array([1, 2, 3]),
            np.array([0, 0, 0]),
            grad=grad_flag,
        )  # All masked

    with pytest.warns(UserWarning, match="Weight sum is not finite"):
        sqrt_sum_wis(
            np.array([2, 2, 2]),
            np.array([1, 2, 3]),
            np.array([1, 1, 1]),
            grad=grad_flag,
        )  # infinite gradient


@pytest.fixture(params=[1, 2, 5])
def increment_percent(request):
    return request.param


@pytest.fixture(params=["K", "J"])
def real_spec(request, use_test_config):
    band = request.param
    wav, flux = load_aces_spectrum([3900, 4.5, 0.0, 0])
    wav, flux = wav_selector(wav, flux, *band_limits(band))
    flux = flux / snr_constant_band(wav, flux, 100, band)
    atm = Atmosphere.from_band(band).at(wav)
    return wav, flux, atm.transmission


@xfail(raises=ModuleNotFoundError, reason="Issue with Starfish install.")
def test_increment_quality_gives_reasonable_length(real_spec, increment_percent):
    """The expected number of steps would be between the
    wavelength difference divided by the
    first and last point * the percent increment.
    """
    print(real_spec)
    wav, flux = real_spec[0], real_spec[1]
    x, q = incremental_quality(wav, flux, percent=increment_percent)
    d1 = wav[0] * increment_percent / 100
    d2 = wav[-1] * increment_percent / 100
    dlambda = wav[-1] - wav[0]
    len_q = len(q)

    assert len_q >= np.floor(dlambda / d1)
    assert len_q <= np.ceil(dlambda / d2 + 1)
    assert len(x) == len_q


@xfail(raises=ModuleNotFoundError, reason="Issue with Starfish install.")
def test_increments_rv_gives_reasonable_length(real_spec, increment_percent):
    """The expected number of steps would be between the
     wavelength difference divided by the
     first and last point * the percent increment.
     """
    wav, flux, mask = real_spec[0], real_spec[1], real_spec[2]
    x, rv = incremental_rv(wav, flux, mask=mask, percent=increment_percent)
    d1 = wav[0] * increment_percent / 100
    d2 = wav[-1] * increment_percent / 100
    dlambda = wav[-1] - wav[0]
    len_rv = len(rv)

    assert len_rv >= np.floor(dlambda / d1)
    assert len_rv <= np.ceil(dlambda / d2 + 1)
    assert len(x) == len_rv


@xfail(raises=ModuleNotFoundError, reason="Issue with Starfish install.")
@pytest.mark.parametrize("no_mask", [True, False])
def test_increments_rv_accumulate_same_as_full(real_spec, increment_percent, no_mask):
    """Assuming that the weighted rv from the steps should equal the rv from the band."""
    wav, flux, mask = real_spec[0], real_spec[1], real_spec[2]
    if no_mask:
        # Try with mask= None also.
        mask = None

    rv_full = rv_precision(wav, flux, mask=mask).value
    x, incremented_rv = incremental_rv(wav, flux, mask=mask, percent=increment_percent)
    incremented_weighted = weighted_error(incremented_rv)

    assert np.round(rv_full, 2) == np.round(incremented_weighted, 2)
    assert x[0] > wav[0]
    assert [-1] < wav[-1]

# Test Qcaclulator

import pytest
import numpy as np
import astropy.units as u
import eniric.Qcalculator as Q
import astropy.constants as const

# test RVprec_calc retuns a single values
# test it returns a quantity in m/s
# test it can handle wavelength and flux alo being quantities.


def test_RVprec_calc():
    """ Test that RVprec_calc can hande inputs as Quantities or unitless and returns a Quantity."""
    wav = np.arange(100)
    flux = np.random.random(100)

    rv = Q.RVprec_calc(wav, flux)
    assert rv.unit == u.meter / u.second
    assert not hasattr(rv.value, '__len__')       # assert value is a scalar
    assert isinstance(rv, u.Quantity)

    # Test it can handle quantities inputs
    rv1 = Q.RVprec_calc(wav * u.micron, flux)
    assert not hasattr(rv1.value, '__len__')       # assert value is a scalar
    assert rv1 == rv
    assert rv1.unit == u.meter / u.second

    rv2 = Q.RVprec_calc(wav * u.micron, (flux / u.second) / (u.centimeter**2))
    assert not hasattr(rv2.value, '__len__')       # assert value is a scalar
    assert rv1 == rv2
    assert rv == rv2
    assert rv2.unit == u.meter / u.second


def test_RVprec_calc_with_lists():
    """ Test that it can hande list input also."""
    wav = list(np.arange(100))
    flux = list(np.random.random(100))

    rv = Q.RVprec_calc(wav, flux)
    assert not hasattr(rv.value, '__len__')       # assert value is a scalar
    assert isinstance(rv, u.Quantity)
    assert rv. unit == u.meter / u.second


def test_SqrtSumWis():
    """ Test that SqrtSumWis can hande inputs as Quantities or unitless
    and returns a dimensionless unscaled Quantity. """
    wav = np.arange(100)
    flux = np.random.random(100)

    sqrtsumwis = Q.SqrtSumWis(wav, flux)
    assert not isinstance(sqrtsumwis, u.Quantity)  # Doesn't turn into quantity if does not have to.
    assert not hasattr(sqrtsumwis, '__len__')       # assert value is a scalar

    sqrtsumwis2 = Q.SqrtSumWis(wav * u.micron,  (flux / u.second) / (u.centimeter**2))   # with some units
    assert not hasattr(sqrtsumwis2.value, '__len__')   # assert value is a scalar
    assert isinstance(sqrtsumwis2, u.Quantity)
    assert sqrtsumwis2.unit == u.dimensionless_unscaled   # unscaled and dimentionless quantitiy

    assert sqrtsumwis == sqrtsumwis2.value

    # Test relation to RVprec_calc
    assert Q.RVprec_calc(wav, flux) == const.c / sqrtsumwis


def test_RV_prec_calc_Trans():

    """ Trans should not have units """
    wav = np.arange(100)
    flux = np.random.random(100)
    trans = np.random.random(100)

    rv_trans = Q.RV_prec_calc_Trans(wav, flux, trans)
    assert not hasattr(rv_trans.value, '__len__')  # assert scalar
    assert rv_trans.unit == u.meter / u.second

    # dimensionless_unscaled unit is ok
    rv_trans2 = Q.RV_prec_calc_Trans(wav, flux, trans * u.dimensionless_unscaled)
    assert not hasattr(rv_trans2.value, '__len__')  # assert  scalar
    assert rv_trans2.unit == u.meter / u.second

    assert rv_trans == rv_trans2

    with pytest.raises(TypeError):
        # transmission mistakenly given as a flux unit
        Q.RV_prec_calc_Trans(wav, flux, (trans/u.s)/(u.centimeter**2))

    with pytest.raises(ValueError):
        Q.RV_prec_calc_Trans(wav, flux, trans+1)

    with pytest.raises(ValueError):
        Q.RV_prec_calc_Trans(wav, flux, trans*-5)


def test_SQRTSumWisTrans():
    """ Test squareroot sum of weights when incuding change of variance due to atmospheric transmission."""
    wav = np.arange(100)
    flux = np.random.random(100)
    trans = np.random.random(100)

    swrtsum_trans = Q.SqrtSumWisTrans(wav, flux, trans)
    assert not isinstance(swrtsum_trans, u.Quantity)  # Doesn't turn into quantity if does not have to.
    assert not hasattr(swrtsum_trans, '__len__')  # assert scalar

    # dimensionless_unscaled unit is ok for transmission
    sqrtsum_trans2 = Q.SqrtSumWisTrans(wav, flux, trans * u.dimensionless_unscaled)
    assert not hasattr(sqrtsum_trans2.value, '__len__')   # assert scalar
    assert isinstance(sqrtsum_trans2, u.Quantity)
    assert sqrtsum_trans2.unit == u.dimensionless_unscaled  # unscaled and dimentionless quantitiy

    sqrtsum_trans3 = Q.SqrtSumWisTrans(wav * u.micron, flux, trans)
    assert not hasattr(sqrtsum_trans3.value, '__len__')  # assert value is a scalar
    assert isinstance(sqrtsum_trans3, u.Quantity)
    assert sqrtsum_trans3.unit == u.dimensionless_unscaled   # unscaled and dimentionless quantitiy

    with pytest.raises(TypeError):
        # transmission mistakenly given as a flux unit
        Q.SqrtSumWisTrans(wav, flux, (trans/u.s)/(u.centimeter**2))

    with pytest.raises(ValueError):
        Q.SqrtSumWisTrans(wav, flux, trans+1)

    with pytest.raises(ValueError):
        Q.SqrtSumWisTrans(wav, flux, trans*-5)


def test_transmission_reduces_precision():
    """Check that a transmission vector reduces precision calcualtion."""
    wav = np.arange(100)
    flux = np.random.random(100)
    transmission = np.random.random(100)
    # Value should be less then normal if trans <=1
    assert Q.RVprec_calc(wav, flux) < Q.RV_prec_calc_Trans(wav, flux, transmission)

    # Unitary transmission should give equivalent result.
    assert Q.RVprec_calc(wav, flux) == Q.RV_prec_calc_Trans(wav, flux, np.ones_like(wav))

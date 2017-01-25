# Test Qcaclulator

import pytest
import numpy as np
import eniric.Qcalculator as Q
import eniric.IOmodule as IO
import astropy.units as u
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
    """ Test that SqrtSumWis can hande inputs as Quantities or unitless and returns a dimensionless unscaled Quantity. """
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

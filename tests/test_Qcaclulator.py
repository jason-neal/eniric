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
    wav = np.arange(1, 101)
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
    wav = np.arange(1, 101)
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
    wav = np.arange(1, 101)
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
    wav = np.arange(1, 101)
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
    wav = np.arange(100.)
    flux = np.random.random(100)
    transmission = np.random.random(100)

    # Value should be less then normal if trans <=1
    assert Q.RVprec_calc(wav, flux) < Q.RV_prec_calc_Trans(wav, flux, transmission)
    # Unitary transmission should give equivalent result.
    assert Q.RVprec_calc(wav, flux) == Q.RV_prec_calc_Trans(wav, flux, np.ones_like(wav))


def test_RV_prec_masked():
    """ Test same prections results between past pre-clumped version and mask version
    """
    wav = np.arange(100)
    flux = np.random.random(100) * 10
    mask = np.asarray(np.floor(2*np.random.random(100)), dtype=bool)

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
    mask = np.asarray(np.floor(2*np.random.random(100)), dtype=bool)

    wav_chunks, flux_chunks = Q.bug_fixed_clumping_method(wav, flux, mask)
    wav_masked, flux_masked = Q.mask_clumping(wav, flux, mask)
    # assert len(wav_chunks) == len(wav_masked)
    # assert len(flux_chunks) == len(flux_masked)
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


def test_bugs_in_old_clumping_method():
    """ Test that it actually works on small tests"""
    val = np.arange(10)

    # Define masks and expected results from val
    mask1 = np.ones(10, dtype=bool)
    expected1 = [val]
    mask2 = np.array([1, 1, 1, 0, 1, 1, 0, 0, 0, 1], dtype=bool)
    expected2 = [np.arange(3), np.arange(4, 6), np.array([9])]

    # Failling examples with bugged code
    mask3 = np.zeros(10, dtype=bool)
    expected3 = []
    unexpected3 = [val]
    mask4 = np.array([0, 1, 0, 0, 0, 1, 1, 1, 1, 0], dtype=bool)
    expected4 = [np.array([1]), np.arange(5, 9)]
    unexpected4 = [np.arange(1), np.arange(2, 5), np.array([9])]

    x1, y1 = Q.bug_fixed_clumping_method(val, val, mask1)
    x1_bugged, y1_bugged = Q.bugged_clumping_method(val, val, mask1)   # This will work
    for i, __ in enumerate(x1):
        assert np.all(x1[i] == y1[i])
        assert np.all(x1[i] == x1_bugged[i])
        assert np.all(y1[i] == y1_bugged[i])
        assert np.all(x1[i] == expected1[i])

    x2, y2 = Q.bug_fixed_clumping_method(val, val, mask2)
    x2_bugged, y2_bugged = Q.bugged_clumping_method(val, val, mask2)  # This will work
    for i, __ in enumerate(x2):
        assert np.all(x2[i] == y2[i])
        assert np.all(x2[i] == x2_bugged[i])
        assert np.all(y2[i] == y2_bugged[i])
        assert np.all(x2[i] == expected2[i])

    # Failing examples where mask starts with 0.
    x3, y3 = Q.bug_fixed_clumping_method(val, val, mask3)
    x3_bugged, y3_bugged = Q.bugged_clumping_method(val, val, mask3)  # This will fail
    for i, __ in x3:
        assert np.all(x3[i] == y3[i])
        assert np.all(x3_bugged[i] == y3_bugged[i])
        assert not np.all(x3[i] == x3_bugged[i])
        assert not np.all(y3[i] == y3_bugged[i])
        assert np.all(x3[i] == expected3[i])
        assert not np.all(x3_bugged[i] == expected3[i])
        assert np.all(x3_bugged[i] == unexpected3[i])

    x4, y4 = Q.bug_fixed_clumping_method(val, val, mask4)
    x4_bugged, y4_bugged = Q.bugged_clumping_method(val, val, mask4)  # This will fail
    for i, __ in enumerate(x4):
        assert np.all(x4[i] == y4[i])
        assert np.all(x4_bugged[i] == y4_bugged[i])
        assert not np.all(x4[i] == x4_bugged[i])
        assert not np.all(y4[i] == y4_bugged[i])
        assert np.all(x4[i] == expected4[i])
        assert not np.all(x4_bugged[i] == expected4[i])
        assert np.all(x4_bugged[i] == unexpected4[i])

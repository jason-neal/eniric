"""Test utilities for eniric."""

import hypothesis.strategies as st
import numpy as np
import pytest
from astropy import constants as const
from hypothesis import given

import eniric.utilities as utils
from eniric.precision import quality
from eniric.utilities import (
    doppler_limits,
    doppler_shift_flux,
    doppler_shift_wav,
    load_aces_spectrum,
    load_btsettl_spectrum,
    mask_between,
    moving_average,
    rv_cumulative,
    rv_cumulative_full,
    weighted_error,
)

xfail = pytest.mark.xfail


@given(
    st.lists(st.floats(allow_nan=False, allow_infinity=False)),
    st.floats(allow_nan=False, allow_infinity=False),
    st.floats(allow_nan=False, allow_infinity=False),
    st.floats(allow_nan=False, allow_infinity=False),
)
def test_wav_selector(x, y, wav_min, wav_max):
    """Test some properties of wavelength selector."""
    y = [xi + y for xi in x]  # just to make y different
    x1, y1 = utils.wav_selector(x, y, wav_min, wav_max)

    assert all(x1 >= wav_min)
    assert all(x1 <= wav_max)
    assert len(x1) == len(y1)
    assert isinstance(x1, np.ndarray)
    assert isinstance(y1, np.ndarray)


@pytest.mark.parametrize("band", ["VIS", "GAP", "z", "Y", "h", "J", "K", "CONT", "NIR"])
def test_band_limits(band):
    """Test getting limits out of band."""
    band_min, band_max = utils.band_limits(band)

    assert band_min < band_max
    assert band_min != band_max
    assert 0 < band_min < 3
    assert 0 < band_max < 3


@pytest.mark.parametrize("band", ["Z", "H", "J", "K"])
def test_band_selector(band):
    """Test band selector selects the wav and flux in the given band."""
    wav = np.linspace(0.5, 3, 100)
    flux = wav ** 2

    band_min, band_max = utils.band_limits(band)
    assert not np.all(wav > band_min)  # Assert wav goes outside band
    assert not np.all(wav < band_max)

    wav, flux = utils.band_selector(wav, flux, band)
    assert np.all(wav > band_min)
    assert np.all(wav < band_max)


@pytest.mark.parametrize(
    "band,error",
    [
        ("X", ValueError),
        ("M0", ValueError),
        (1, AttributeError),
        (np.array(1), AttributeError),
        (["list", "of", "strings"], AttributeError),
    ],
)
def test_band_limits_raises_errors(band, error):
    """Test it raises the Value and Attribute Errors."""
    with pytest.raises(error):
        utils.band_limits(band)


@pytest.mark.parametrize(
    "band,error",
    [
        ("X", ValueError),
        ("M0", ValueError),
        (1, AttributeError),
        (["list", "of", "strings"], AttributeError),
        (np.linspace(1, 2, 10), AttributeError),
    ],
)
def test_band_selector_raises_errors(band, error):
    """Test it raises the Value and Attribute Errors"""
    wav = np.linspace(0.5, 3, 100)
    flux = wav ** 2

    with pytest.raises(error):
        utils.band_selector(wav, flux, band)


@pytest.mark.parametrize("band", ["ALL", ""])
def test_band_selector_with_no_selection(band):
    """If band = "ALL" or ""."""
    wav = np.linspace(0.5, 3, 100)
    flux = wav ** 2
    wav2, flux2 = utils.band_selector(wav, flux, band)

    # No changes
    assert np.all(wav == wav2)
    assert np.all(flux2 == flux)


@pytest.mark.parametrize("band", ["H", "J", "K", "VIS", "Z"])
def test_band_middle(band):
    """Test band middle is middle of band.
    Some repeated coding.
    """
    lower, upper = utils.band_limits(band)
    middle = utils.band_middle(band)

    assert lower < middle
    assert middle < upper
    assert (lower + upper) / 2 == middle


def test_band_midpoint_j():
    """Middle of J is 1.25 microns."""
    assert utils.band_middle("J") == 1.25


def test_silent_remove():
    """Test this doesn't raise and issue.
    Not really a good test."""
    utils.silent_remove("a_fake_filename_that_doesnt_exist.fake")
    assert True


####################################################
# Test Resolution conversions
@pytest.mark.parametrize(
    "resolution,result", [("60k", 60000), (80000, 80000), ("2000", 2000)]
)
def test_res2int(resolution, result):
    assert result == utils.res2int(resolution)


@pytest.mark.parametrize(
    "resolution,result",
    [
        (60000, "60k"),
        ("20000", "20k"),
        (np.float("20000"), "20k"),
        ("60k", "60k"),
        (80000, "80k"),
        (3000, "3k"),
    ],
)
def test_res2str(resolution, result):
    """Test single values in res2str"""
    assert result == utils.res2str(resolution)


@pytest.mark.parametrize(
    "resolutions,results",
    [
        ([60000, 80000, 100_000], [60000, 80000, 100_000]),
        (["60k", "80k", "100k"], [60000, 80000, 100_000]),
        (["6000", "8000", "100000"], [6000, 8000, 100_000]),
        (["10000", "20k", "300K"], [10000, 20000, 300_000]),
    ],
)
def test_resolutions2ints_lists(resolutions, results):
    """Test transformation of resolutions to integer resolutions."""
    assert results == utils.resolutions2ints(resolutions)


@pytest.mark.parametrize(
    "resolutions,results",
    [
        (["60k", "80k", "100k"], ["60k", "80k", "100k"]),
        ([60000, 80000, 100_000], ["60k", "80k", "100k"]),
        (["60000", "80K", "100k"], ["60k", "80k", "100k"]),
        ([np.float("60000"), np.int("2000")], ["60k", "2k"]),
    ],
)
def test_resolutions2strs_list(resolutions, results):
    """Test transformation of a list of resolutions to strings."""
    assert results == utils.resolutions2strs(resolutions)


@given(st.integers(min_value=1, max_value=1_000_000))
def test_res2int_doesnt_change_int(resolution):
    assert resolution == utils.res2int(resolution)


@pytest.mark.parametrize(
    "resolution", [1000, 10000, 200_000, "1k", "10k", "100k", "2000k"]
)
def test_compatibility_res2int_res2str(resolution):
    """Test res2int and rest2str reversible and do not change when operated twice.

    These are single values.
    """
    resolution = resolution
    res2str = utils.res2str
    res2int = utils.res2int

    assert res2str(resolution) == res2str(res2int(resolution))
    assert res2str(resolution) == res2str(res2str(resolution))

    assert res2int(resolution) == res2int(res2str(resolution))
    assert res2int(resolution) == res2int(res2int(resolution))


@pytest.mark.parametrize(
    "resolution",
    [
        [1000, 10000, 50000],
        [100_000, 2_000_000],
        ["1k"],
        ["10k", "20k"],
        ["100k", "2000k"],
    ],
)
def test_compatibility_resolutions2ints_resolutions2strs(resolution):
    """Test resolutions2ints and resolutions2strs reversible and do not change when operated twice.

    These are lists of values.
    """

    res2str = utils.resolutions2strs
    res2int = utils.resolutions2ints

    assert res2str(resolution) == res2str(res2int(resolution))
    assert res2str(resolution) == res2str(res2str(resolution))

    assert res2int(resolution) == res2int(res2str(resolution))
    assert res2int(resolution) == res2int(res2int(resolution))


@pytest.mark.parametrize(
    "resolutions",
    [[60000, "20000"], ["50k", "10k"], ["100000", "10000"], ["100000", "10000"]],
)
def test_res2int_fails_on_list(resolutions):
    with pytest.raises(TypeError):
        utils.res2int(resolutions)


@pytest.mark.parametrize("resolutions", [60000, "10k" "100000"])
def test_resolutions2ints_fails_on_single(resolutions):
    with pytest.raises(TypeError):
        utils.resolutions2ints(resolutions)


@pytest.mark.parametrize(
    "resolutions",
    [[60000, "20000"], ["50k", "10k"], ["100000", "10000"], ["100000", "10000"]],
)
def test_res2str_fails_on_list(resolutions):
    with pytest.raises(TypeError):
        utils.res2str(resolutions)


@pytest.mark.parametrize("resolutions", [60000, "10k" "100000"])
def test_resolutions2strs_fails_on_single(resolutions):
    with pytest.raises(TypeError):
        utils.resolutions2strs(resolutions)


###############################################


@given(
    st.lists(st.floats(allow_infinity=False, allow_nan=False), max_size=50),
    st.floats(1, 100),
    st.floats(1, 100),
)
def test_mask_between(x, x1, x2):
    """Correctly masks out values."""
    x = np.array(x)
    xmin = min(x1, x2)
    xmax = max(x1, x2)
    mask = mask_between(x, xmin, xmax)

    assert np.all(x[mask] >= xmin)
    assert np.all(x[mask] < xmax)
    # Dropped values are all outside range.
    assert np.all((x[~mask] >= xmax) | (x[~mask] < xmin))


@pytest.mark.parametrize(
    "input_, flag",
    [
        ([1, 2, 3, 4.5, 5], False),
        ([1, 2, 3, 4, 5], True),
        ([10, 2.5, 12, 4.5, 0.1], False),
        ([10, 2.5, 12.5, 4, 0.05], True),
    ],
)
def test_rv_cumulative(input_, flag):
    """Test it works on well formed input."""
    result = rv_cumulative(input_, single=flag)
    assert result[-2] == weighted_error(input_[:-1])
    assert result[-1] == weighted_error(input_)
    if flag:
        assert result[0] == input_[0]


@pytest.mark.parametrize(
    "input_", [[1, 2, 3, 4, 5], [10, 3, 12, 4, 50], [10, 66, 10.5, 4.5, 1]]
)
def test_rv_cumulative_full(input_):
    result = rv_cumulative_full(input_)
    assert len(result) == 9
    assert result[4] == weighted_error(input_)
    assert result[0] == input_[0]
    assert result[-1] == input_[-1]
    assert isinstance(result, np.ndarray)


@pytest.mark.parametrize("input_", [[1, 5], [3, 4, 5], [1, 1, 3, 7, 11, 16], []])
def test_rv_cumulative_full_errors_on_size(input_):
    """Test it errors when not given 5 values."""
    with pytest.raises(AssertionError):
        rv_cumulative_full(input_)


@pytest.mark.parametrize("window_size", [1, 5, 100])
def test_moving_average_size(window_size):
    x = np.arange(300)
    result = moving_average(x, window_size)
    assert isinstance(result, np.ndarray)
    assert len(result) == len(x)


@pytest.mark.parametrize(
    "input_, expected",
    [
        ([1, 1, 1], 1 / np.sqrt(3)),
        ([2, 2, 2], np.sqrt(4.0 / 3)),
        ([1], 1),
        ([1, 2, 3, 4, 5], 60 / np.sqrt(5269)),
    ],
)
def test_weighted_error(input_, expected):
    result = weighted_error(input_)
    assert np.allclose(result, expected)


@pytest.fixture(
    params=[
        np.array([1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8]),
        np.linspace(2.105, 2.115, 50),
        np.linspace(1100, 1105, 70),
        np.linspace(300, 305, 100),
    ]
)
def wavelength(request):
    wave = request.param
    return wave


@pytest.mark.parametrize("direction,multiplier", [(-1, 0), (1, 2)])
def test_doppler_shift_at_speed_of_light(wavelength, direction, multiplier):
    r"""It is test is only valid for the non-relativistic doppler shift.
     It is a test of the math but not physical (not valid at this speed)
     but we are checking the math of the equation works
     should result in a shift of \delta \lambda = +/- \lambda"""
    c_wave = doppler_shift_wav(wavelength, direction * const.c.to("km/s").value)
    assert np.all(c_wave == (wavelength * multiplier))


def doppler_shift_zero(wavelength):
    """Doppler shift of zero will give no change."""
    assert np.all(doppler_shift_wav(wavelength, vel=0.0) == wavelength)


@xfail(raises=ModuleNotFoundError, reason="Starfish is not installled correctly.")
@pytest.mark.parametrize("rv", [-100, 200])
def test_if_doppler_shift_changes_quality(testing_spectrum, rv):
    wavelength, flux_ = testing_spectrum[0], testing_spectrum[1]
    wmin_, wmax_ = 2.1, 2.3
    m1 = (wavelength < wmax_ + 0.2) * (wavelength > wmin_ - 0.2)
    # shorten spectra
    wavelength, flux_ = wavelength[m1], flux_[m1]

    mask = (wavelength < wmax_) * (wavelength > wmin_)
    wave = wavelength[mask]
    flux = flux_[mask]
    q1 = quality(wave, flux)

    wav2 = doppler_shift_wav(wavelength, rv)
    mask2 = (wav2 < wmax_) * (wav2 > wmin_)
    wav2 = wav2[mask2]
    flux2 = flux_[mask2]
    q2 = quality(wav2, flux2)
    assert q1 != q2

    flux3 = doppler_shift_flux(wavelength, flux_, rv, new_wav=wave)
    q3 = quality(wave, flux3)
    assert len(wave) == len(flux3)
    assert q1 != q3
    # There are differences due to interpolation
    assert q2 != q3


@pytest.mark.parametrize("rv", [-10, 50, 1000])
@pytest.mark.parametrize("wmin, wmax", [(2.1, 2.2), (1500, 1700)])
def test_doppler_limits(rv, wmin, wmax):
    """Doppler limits widens the wavelength range"""
    new_min, new_max = doppler_limits(rv, wmin, wmax)
    assert new_min < wmin
    assert new_max > wmax


@pytest.mark.parametrize("wmin, wmax", [(2.1, 2.2), (1500, 1700)])
def test_doppler_limits_rv_0(wmin, wmax):
    """RV of zero should have no effect to limits."""
    new_min, new_max = doppler_limits(0, wmin, wmax)
    assert new_min == wmin
    assert new_max == wmax


@xfail(raises=ModuleNotFoundError, reason="Issue with Starfish install.")
@pytest.mark.parametrize("photons", [True, False])
def test_load_btsettl_spectrum(photons, use_test_config):
    wav, flux = load_btsettl_spectrum(
        [3900, 4.5, 0, 0], photons=photons, air=False, wl_range=[21000, 22000]
    )
    assert len(wav) == len(flux)


@xfail(raises=ModuleNotFoundError, reason="Issue with Starfish install.")
@pytest.mark.parametrize("params", [[3900, 4.5, 1, 0], [3900, 4.5, 0, 1]])
def test_invalid_feh_alpha_load_btsettl_spectrum(params):
    # Invalid CIFIST parameters
    with pytest.raises(AssertionError):
        load_btsettl_spectrum(params, wl_range=[21000, 22000])


@xfail(raises=ModuleNotFoundError, reason="Issue with Starfish install.")
@pytest.mark.parametrize("photons", [True, False])
def test_load_aces_spectrum(photons, use_test_config):
    wav, flux = load_aces_spectrum(
        [3900, 4.5, 0, 0], photons=photons, air=False, wl_range=[21000, 22000]
    )
    assert len(wav) == len(flux)


@xfail(raises=ModuleNotFoundError, reason="Issue with Starfish install.")
@pytest.mark.parametrize(
    "params",
    [
        [20000, 4.5, 0, 0],
        [2200, 4.5, 0, 0],
        [3900, 4.7, 1, 0],
        [3900, 4.5, 0.2],
        [3900, 4.5, -1, 3],
    ],
)
def test_invalid_load_aces_spectrum(params):
    with pytest.raises(ValueError):
        load_aces_spectrum(params, wl_range=[21000, 22000])


@xfail(raises=ModuleNotFoundError, reason="Issue with Starfish install.")
@pytest.mark.parametrize("params", [[8000, 4.5, 0, 0], [3900, 0.5, 0, 0]])
def test_invalid_load_btsettl_spectrum(params):
    with pytest.raises(ValueError):
        load_btsettl_spectrum(params, wl_range=[21000, 22000])

"""Test utilities for eniric."""

import os

import hypothesis.strategies as st
import numpy as np
import pytest
from hypothesis import given, settings

import eniric
import eniric.utilities as utils
from eniric.utilities import mask_between


# @pytest.mark.xfail(raises=FileNotFoundError)
def test_read_spectrum():
    """Test reading in a _wave_photon.dat is the same as a _wave.dat."""
    photon = os.path.join(
        eniric.paths["test_data"],
        "sample_lte03900-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes_wave_photon.dat",
    )
    wave = os.path.join(
        eniric.paths["test_data"],
        "sample_lte03900-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes_wave.dat",
    )
    wave_wav, wave_flux = utils.read_spectrum(wave)
    photon_wav, photon_flux = utils.read_spectrum(photon)

    assert np.allclose(photon_wav, wave_wav)
    assert np.allclose(photon_flux, wave_flux)


# @pytest.mark.xfail(raises=FileNotFoundError)
def test_get_spectrum_name():
    """Test specifying file names with stellar parameters."""
    test = os.path.join("Z-0.0", "lte02800-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes_wave.dat")

    assert utils.get_spectrum_name("M6", flux_type="wave") == test

    test_alpha = os.path.join(
        "Z-0.0.Alpha=+0.20",
        "lte02600-6.00-0.0.Alpha=+0.20.PHOENIX-ACES-AGSS-COND-2011-HiRes_wave_photon.dat",
    )
    assert utils.get_spectrum_name("M9", logg=6, alpha=0.2) == test_alpha

    test_pos_feh = os.path.join(
        "Z+0.5", "lte03500-0.00+0.5.PHOENIX-ACES-AGSS-COND-2011-HiRes_wave_photon.dat"
    )
    assert utils.get_spectrum_name("M3", logg=0, feh=0.5, alpha=0.0) == test_pos_feh

    test_photon = os.path.join(
        "Z-0.0", "lte02800-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes_wave_photon.dat"
    )
    assert utils.get_spectrum_name("M6") == test_photon


# noinspection SpellCheckingInspection
@pytest.mark.parametrize("spec_type", ["MO", "ME", "M11", "X10", "Z3"])
def test_spectrum_name_value_error(spec_type):
    """Not valid spectral type in [OBAFGKML] or misspelled"""
    with pytest.raises(ValueError):
        utils.get_spectrum_name(spec_type)


@pytest.mark.parametrize("spec_type", ["O1", "B2", "A3", "F4", "G5", "K6", "M7", "L8"])
def test_notimplemented_spectrum_name(spec_type):
    with pytest.raises(NotImplementedError):
        utils.get_spectrum_name(spec_type)  # Stellar type not added (only M atm)


@pytest.mark.parametrize("bad_alpha", [-0.3, 0.3, 1])
def test_spectrum_name_with_bad_alpha(bad_alpha):
    """Bad_alpha is outside range -0.2-0.2 for M-dwarf science case."""
    with pytest.raises(ValueError):
        utils.get_spectrum_name("M0", alpha=bad_alpha)


@pytest.mark.parametrize("alpha", [-0.2, 0.1, 0.2])
def test_spectrum_name_with_ok_alpha(alpha):
    name = utils.get_spectrum_name("M0", alpha=alpha)

    assert isinstance(name, str)
    assert str(alpha) in name
    assert "Alpha=" in name


# @pytest.mark.xfail(raises=FileNotFoundError)
def test_org_name():
    """Test org flag of utils.get_spectrum_name, supposed to be temporary."""
    test_org = "lte03900-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes_wave.dat"
    assert utils.get_spectrum_name("M0", org=True) == test_org


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


##################################3


@settings(max_examples=100)
@given(
    st.lists(
        st.floats(
            min_value=1e-7, max_value=1e-5, allow_infinity=False, allow_nan=False
        ),
        unique=True,
        min_size=3,
        max_size=25,
    ),
    st.floats(min_value=1e-2, max_value=200),
    st.floats(min_value=1e-4, max_value=1),
)
def test_rotational_kernel(delta_lambdas, vsini, epsilon):
    """Test that the new and original code produces the same output."""
    delta_lambdas = np.sort(np.asarray(delta_lambdas), kind="quicksort")
    delta_lambdas = np.append(np.flipud(delta_lambdas), np.insert(delta_lambdas, 0, 0))
    delta_lambda_l = np.max(delta_lambdas) * 2

    new_profile = utils.rotation_kernel(delta_lambdas, delta_lambda_l, vsini, epsilon)

    assert len(new_profile) == len(delta_lambdas)
    # other properties to test?


@given(
    st.lists(
        st.floats(min_value=-100, max_value=100, allow_nan=False),
        min_size=1,
        unique=True,
    ),
    st.floats(min_value=-100, max_value=100, allow_nan=False),
    st.floats(min_value=0.001, max_value=100, allow_nan=False),
)
def test_unitary_gaussian(x, center, fwhm):
    """Just a quick simple test."""
    x = np.asarray(x)

    gaussian = utils.unitary_gaussian(x, center, fwhm)
    print(gaussian)
    # point at center should be the max
    assert len(gaussian) == len(x)
    assert np.allclose(np.max(gaussian), gaussian[np.argmin(abs(x - center))])


def test_unitary_gaussian_type_errors():
    """Testing for type errors."""
    x = np.arange(-10, 10)
    center = 0
    fwhm = 3

    with pytest.raises(TypeError):
        utils.unitary_gaussian(x, center, "fwhm")
    with pytest.raises(TypeError):
        utils.unitary_gaussian(x, "center", fwhm)
    with pytest.raises(TypeError):
        utils.unitary_gaussian(range(-10, 10), "center", fwhm)
    with pytest.raises(TypeError):
        utils.unitary_gaussian(1, "center", fwhm)


################################################

def test_silent_remove():
    """Test this doesn't raise and issue.
    Not really a good test."""
    utils.silent_remove("a_fake_filename_that_doesnt_exist.fake")
    assert True


####################################################
# Test Resolution conversions

@pytest.mark.parametrize(
    "resolution,result",
    [
        ("60k", 60000),
        (80000, 80000),
        ("2000", 2000),
    ],
)
def test_res2int(resolution, result):
    assert result == utils.res2int(resolution)


@pytest.mark.parametrize(
    "resolution,result", [(60000, "60k"),
                          ("20000", "20k"),
                          (np.float("20000"), "20k"),
                          ("60k", "60k"),
                          (80000, "80k"),
                          (3000, "3k")]
)
def test_res2str(resolution, result):
    """Test single values in res2str"""
    assert result == utils.res2str(resolution)


@pytest.mark.parametrize(
    "resolutions,results",
    [
        ([60000, 80000, 100000], [60000, 80000, 100000]),
        (["60k", "80k", "100k"], [60000, 80000, 100000]),
        (["6000", "8000", "100000"], [6000, 8000, 100000]),
        (["10000", "20k", "300K"], [10000, 20000, 300000]),
    ],
)
def test_resolutions2ints_lists(resolutions, results):
    """Test transformation of resolutions to integer resolutions."""
    assert results == utils.resolutions2ints(resolutions)


@pytest.mark.parametrize(
    "resolutions,results",
    [
        (["60k", "80k", "100k"], ["60k", "80k", "100k"]),
        ([60000, 80000, 100000], ["60k", "80k", "100k"]),
        (["60000", "80K", "100k"], ["60k", "80k", "100k"]),
        ([np.float("60000"), np.int("2000")], ["60k", "2k"]),
    ],
)
def test_resolutions2strs_list(resolutions, results):
    """Test transformation of a list of resolutions to strings."""
    assert results == utils.resolutions2strs(resolutions)


@given(st.integers(min_value=1, max_value=1000000))
def test_res2int_doesnt_change_int(resolution):
    assert resolution == utils.res2int(resolution)


@pytest.mark.parametrize(
    "resolution", [1000, 10000, 200000, "1k", "10k", "100k", "2000k"]
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
    [[1000, 10000, 50000], [100000, 2000000], ["1k"], ["10k", "20k"], ["100k", "2000k"]]
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


@pytest.mark.parametrize("resolutions", [
    [60000, "20000"],
    ["50k", "10k"],
    ["100000", "10000"],
    ["100000", "10000"],
])
def test_res2int_fails_on_list(resolutions):
    with pytest.raises(TypeError):
        utils.res2int(resolutions)


@pytest.mark.parametrize("resolutions", [
    60000,
    "10k"
    "100000",
])
def test_resolutions2ints_fails_on_single(resolutions):
    with pytest.raises(TypeError):
        utils.resolutions2ints(resolutions)


@pytest.mark.parametrize("resolutions", [
    [60000, "20000"],
    ["50k", "10k"],
    ["100000", "10000"],
    ["100000", "10000"],
])
def test_res2str_fails_on_list(resolutions):
    with pytest.raises(TypeError):
        utils.res2str(resolutions)


@pytest.mark.parametrize("resolutions", [
    60000,
    "10k"
    "100000",
])
def test_resolutions2strs_fails_on_single(resolutions):
    with pytest.raises(TypeError):
        utils.resolutions2strs(resolutions)


###############################################
@pytest.mark.parametrize(
    "filename",
    [
        os.path.join(
            eniric.paths["test_data"],
            "results",
            "Spectrum_M0-PHOENIX-ACES_Kband_vsini1.0_R100k.txt",
        ),
        os.path.join(
            eniric.paths["test_data"],
            "resampled",
            "Spectrum_M0-PHOENIX-ACES_Kband_vsini1.0_R100k_res3.0.txt",
        ),
    ],
)
def test_resampled_spectra_isnot_read_by_read_spectrum(filename):
    """Doesn't allow names with _vsini or _res in them."""
    with pytest.raises(ValueError, match="Using wrong function"):
        utils.read_spectrum(filename)


@given(st.lists(st.floats(allow_infinity=False, allow_nan=False), max_size=50), st.floats(1, 100),
       st.floats(1, 100))
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

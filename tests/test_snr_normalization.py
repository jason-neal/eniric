import numpy as np
import pytest

import eniric.snr_normalization as snrnorm
import eniric.utilities as utils
from scripts.phoenix_precision import convolve_and_resample

xfail = pytest.mark.xfail


@xfail(raises=ModuleNotFoundError, reason="Issue with Starfish install.")
@pytest.mark.parametrize("desired_snr", [100.0, 150.0])
@pytest.mark.parametrize("band", ["J", "Y", "VIS"])
def test_snr_normalization(desired_snr, band, testing_spectrum):
    """Test SNR after normalizing function is the desired value.

    Testing on middle of J band.
    """
    wav, flux = testing_spectrum
    band_mid = utils.band_middle(band)

    # Searching for the closest index to 1.25
    index_reference = np.searchsorted(wav, [band_mid])[0]
    snr_estimate = np.sqrt(np.sum(flux[index_reference - 1 : index_reference + 2]))

    assert round(snr_estimate, 1) != desired_snr  # Assert SNR is not correct

    norm_const = snrnorm.snr_constant_band(wav, flux, snr=desired_snr, band=band)

    new_flux = flux / norm_const
    new_snr_estimate = np.sqrt(
        np.sum(new_flux[index_reference - 1 : index_reference + 2])
    )

    assert round(new_snr_estimate, 0) == desired_snr


@xfail(raises=ModuleNotFoundError, reason="Issue with Starfish install.")
@pytest.mark.parametrize("desired_snr", [50.0, 200.0])
@pytest.mark.parametrize("band", ["H", "GAP", "K"])
def test_snr_normalization_constant(desired_snr, band, testing_spectrum):
    """Test snr_constant_band and snr_constant_wav produce same result."""
    wav, flux = testing_spectrum
    band_mid = utils.band_middle(band)

    assert snrnorm.snr_constant_band(
        wav, flux, band=band, snr=desired_snr
    ) == snrnorm.snr_constant_wav(wav, flux, band_mid, snr=desired_snr)


@xfail(raises=ModuleNotFoundError, reason="Issue with Starfish install.")
@pytest.mark.parametrize("sampling", [3, 4.1])
def test_band_snr_norm(testing_spectrum, sampling):
    """Compared to wav snr norm."""
    wav, flux = testing_spectrum
    wav, flux = convolve_and_resample(
        wav, flux, vsini=1, R=100_000, band="J", sampling=sampling
    )
    assert snrnorm.snr_constant_band(
        wav, flux, band="J", snr=100, sampling=sampling
    ) == snrnorm.snr_constant_wav(wav, flux, wav_ref=1.25, snr=100, sampling=sampling)

    assert snrnorm.snr_constant_band(
        wav, flux, band="J", snr=100, sampling=sampling, verbose=True
    ) != snrnorm.snr_constant_wav(
        wav, flux, wav_ref=1.24, snr=100, sampling=sampling, verbose=True
    )


def test_sampling_index():
    """Some hard coded examples of sampling index."""
    # odd number
    assert snrnorm.sampling_index(6, 1) == [6]
    assert np.all(snrnorm.sampling_index(100, 3) == [99, 100, 101])
    assert np.all(snrnorm.sampling_index(5, 5) == [3, 4, 5, 6, 7])

    # even number
    assert np.all(snrnorm.sampling_index(10, 4) == [8, 9, 10, 11])
    assert np.all(snrnorm.sampling_index(10, 2) == [9, 10])
    # number that is at end of array.


@pytest.mark.parametrize("sample", [1, 2, 3, 4, 5, 7, 10, 15])
def test_sampling_size_and_type(sample):
    # Check number of values correct.
    for sample in [1, 2, 3, 4, 5, 7, 10, 15]:
        assert len(snrnorm.sampling_index(20, sample)) == sample
        assert isinstance(
            snrnorm.sampling_index(20, sample)[0], (int, np.integer)
        )  # index values must be int


def test_sampling_index_array():
    """Sampling index when array_length is given, or when index goes out of bounds."""
    assert np.all(snrnorm.sampling_index(100, 3, array_length=200) == [99, 100, 101])

    with pytest.raises(ValueError):
        snrnorm.sampling_index(3, 10)  # index will be < 0
    with pytest.raises(ValueError):
        snrnorm.sampling_index(3, 9, array_length=50)  # index will be < 0
    with pytest.raises(ValueError):
        snrnorm.sampling_index(
            46, 10, array_length=50
        )  # an index will be > (array_length - 1)


@pytest.mark.parametrize("band", ["VIS", "Z", "NIR", "J", "Y"])
def test_snr_constant_band_returns_mid_value_const(band):
    size = 100
    np.random.seed(40)
    flux = 500 * np.random.rand(
        size
    )  # To give a random spectrum (but consistent between tests)
    lim = utils.band_limits(band)
    wav = np.linspace(lim[0], lim[1], size)

    band_const = snrnorm.snr_constant_band(wav, flux, band=band)
    wav_const = snrnorm.snr_constant_wav(wav, flux, wav_ref=utils.band_middle(band))

    assert isinstance(band_const, float)
    assert isinstance(wav_const, float)
    assert band_const == wav_const  # Since band calls wave at midpoint


@pytest.mark.parametrize("band", ["VIS", "K", "H"])
@pytest.mark.parametrize("verbose", [True, False])
def test_snr_normalization_logic(band, verbose):
    """Testing direct value.

    snr = sqrt(sum(3 pixels))
    const = (snr/ref_snr)**2
    if pixel value = 3 then  the normalization constant will be 1.
    """
    size = 100
    band = "K"
    lim = utils.band_limits(band)
    wav = np.linspace(lim[0], lim[1], size)
    flux = 3 * np.ones(size)
    band_const = snrnorm.snr_constant_band(wav, flux, snr=3, band=band, verbose=verbose)
    wav_const = snrnorm.snr_constant_wav(
        wav, flux, snr=3, wav_ref=(lim[0] + lim[1]) / 2, verbose=verbose
    )
    assert band_const == 1
    assert wav_const == 1


@pytest.mark.parametrize(
    "wav,band",
    [
        (np.linspace(0.8, 1, 50), "VIS"),  # "VIS": (0.38, 0.78)
        (np.linspace(2, 3, 50), "J"),  # "J": (1.17, 1.33)
        (np.linspace(2.0, 2.1, 50), "K"),  # "K": (2.07, 2.35)
        (np.linspace(2.25, 2.4, 50), "K"),  # "K": (2.07, 2.35)
    ],
)
def test_snr_constant_band_with_invalid_wavelength(wav, band):
    with pytest.raises(ValueError):
        snrnorm.snr_constant_band(wav, np.ones(50), band=band)


@pytest.mark.parametrize("wav_ref", [0.5, 4])
def test_snr_constant_wav_ref_outside_wav(wav_ref):
    """Wav-ref outside bounds of wav should raise ValueError"""
    wav = np.linspace(1, 3, 60)
    flux = np.random.randn(len(wav))

    with pytest.raises(ValueError):
        snrnorm.snr_constant_wav(wav, flux, wav_ref)

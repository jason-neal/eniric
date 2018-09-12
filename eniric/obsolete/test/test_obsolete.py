"""Move in obsolete testing functions."""
import os

import numpy as np
import pytest

import eniric
from eniric import (
    Qcalculator as Q,
    io_module as io,
    snr_normalization as snrnorm,
    utilities as utils,
)
from eniric.obsolete.nIR_run import main as nir_run
from eniric.obsolete.utilities import barycenter_shift

resampled_template = "Spectrum_{0}-PHOENIX-ACES_{1}band_vsini{2}_R{3}_res3.0.dat"
wave_photon_template = (
    "lte0{0}-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes_wave_photon.dat"
)


def test_get_spectrum_name():
    """Test specifying file names with stellar parameters."""
    test = os.path.join(
        "Z-0.0", "lte02800-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes_wave.dat"
    )

    assert eniric.obsolete.utilities.get_spectrum_name("M6", flux_type="wave") == test

    test_alpha = os.path.join(
        "Z-0.0.Alpha=+0.20",
        "lte02600-6.00-0.0.Alpha=+0.20.PHOENIX-ACES-AGSS-COND-2011-HiRes_wave_photon.dat",
    )
    assert (
        eniric.obsolete.utilities.get_spectrum_name("M9", logg=6, alpha=0.2)
        == test_alpha
    )

    test_pos_feh = os.path.join(
        "Z+0.5", "lte03500-0.00+0.5.PHOENIX-ACES-AGSS-COND-2011-HiRes_wave_photon.dat"
    )
    assert (
        eniric.obsolete.utilities.get_spectrum_name("M3", logg=0, feh=0.5, alpha=0.0)
        == test_pos_feh
    )

    test_photon = os.path.join(
        "Z-0.0", "lte02800-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes_wave_photon.dat"
    )
    assert eniric.obsolete.utilities.get_spectrum_name("M6") == test_photon


@pytest.mark.parametrize("spec_type", ["MO", "ME", "M11", "X10", "Z3"])
def test_spectrum_name_value_error(spec_type):
    """Not valid spectral type in [OBAFGKML] or misspelled"""
    with pytest.raises(ValueError):
        eniric.obsolete.utilities.get_spectrum_name(spec_type)


@pytest.mark.parametrize("spec_type", ["O1", "B2", "A3", "F4", "G5", "K6", "M7", "L8"])
def test_notimplemented_spectrum_name(spec_type):
    with pytest.raises(NotImplementedError):
        eniric.obsolete.utilities.get_spectrum_name(
            spec_type
        )  # Stellar type not added (only M atm)


@pytest.mark.parametrize("bad_alpha", [-0.3, 0.3, 1])
def test_spectrum_name_with_bad_alpha(bad_alpha):
    """Bad_alpha is outside range -0.2-0.2 for M-dwarf science case."""
    with pytest.raises(ValueError):
        eniric.obsolete.utilities.get_spectrum_name("M0", alpha=bad_alpha)


@pytest.mark.parametrize("alpha", [-0.2, 0.1, 0.2])
def test_spectrum_name_with_ok_alpha(alpha):
    name = eniric.obsolete.utilities.get_spectrum_name("M0", alpha=alpha)

    assert isinstance(name, str)
    assert str(alpha) in name
    assert "Alpha=" in name


def test_org_name():
    """Test org flag of utils.get_spectrum_name, supposed to be temporary."""
    test_org = "lte03900-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes_wave.dat"
    assert eniric.obsolete.utilities.get_spectrum_name("M0", org=True) == test_org


@pytest.mark.parametrize("noresample", [True, False])
@pytest.mark.parametrize("unnormalized", [True, False])
@pytest.mark.parametrize("org", [True, False])
def test_nir_run_raises_type_errors_on_non_lists(noresample, unnormalized, org):
    """Checking list inputs are needed.

    Checking over range of other boolean flags.
    """
    # Initialize parameters
    startype = ["M0", "M1"]
    vsini = [1]
    res = ["100k"]
    band = ["J", "K"]
    sample_rate = None

    # Check for TypeError on each parameter
    with pytest.raises(TypeError):
        nir_run(
            startype="M0",
            vsini=vsini,
            resolution=res,
            band=band,
            sample_rate=sample_rate,
            noresample=noresample,
            unnormalized=unnormalized,
            org=org,
        )

    with pytest.raises(TypeError):
        nir_run(
            startype=startype,
            vsini=1.5,
            resolution=res,
            band=band,
            sample_rate=sample_rate,
            noresample=noresample,
            unnormalized=unnormalized,
            org=org,
        )

    with pytest.raises(TypeError):
        nir_run(
            startype=startype,
            vsini=vsini,
            resolution="100k",
            band=band,
            sample_rate=sample_rate,
            noresample=noresample,
            unnormalized=unnormalized,
            org=org,
        )

    with pytest.raises(TypeError):
        nir_run(
            startype=startype,
            vsini=vsini,
            resolution=res,
            band="K",
            sample_rate=sample_rate,
            noresample=noresample,
            unnormalized=unnormalized,
            org=org,
        )

    with pytest.raises(TypeError):
        nir_run(
            startype=startype,
            vsini=vsini,
            resolution=res,
            band=band,
            sample_rate=5,
            noresample=noresample,
            unnormalized=unnormalized,
            org=org,
        )


@pytest.mark.parametrize("bad_string", ["id-string", "M0-K-1.0-100", "M0-P-1.0-100k"])
def test_errors_in_snr_get_reference_spectrum(bad_string):
    """Testing Errors in getting the reference spectrum."""
    with pytest.raises(ValueError):
        eniric.obsolete.snr_norm.get_reference_spectrum(bad_string)


@pytest.mark.parametrize("bad_string", ["Alpha=", "smpl="])
def test_notimplemented_errors_in_snr_get_reference_spectrum(bad_string):
    """Testing getting the reference spectrum.

    Currently "Alpha=" in the id-string is not implemented.
    Currently "smpl=" in the id-string is not implemented.
    """
    with pytest.raises(NotImplementedError):
        eniric.obsolete.snr_norm.get_reference_spectrum(bad_string)


def test_valid_snr_get_reference_spectrum():
    """Testing getting the reference spectrum."""
    ref_band = "J"
    wav_ref, flux_ref = eniric.obsolete.snr_norm.get_reference_spectrum(
        "M0-K-1.0-100k", ref_band=ref_band
    )
    band_min, band_max = utils.band_limits(ref_band)

    # Test the wavelength is in the reference band wavelength range
    assert np.all(wav_ref <= band_max)
    assert np.all(wav_ref >= band_min)

    # test properties of output
    assert len(wav_ref) == len(flux_ref)
    assert isinstance(wav_ref, np.ndarray)
    assert isinstance(flux_ref, np.ndarray)


def test_get_reference_spectrum_in_nonexistent_file():
    """Testing getting the reference spectrum."""
    with pytest.raises(FileNotFoundError):
        eniric.obsolete.snr_norm.get_reference_spectrum("M1-K-1.0-100k", ref_band="J")


def test_normalize_flux_new_verse_old(resampled_data):
    """Test only small differences due to new normalization."""
    id_string, wav, flux = resampled_data

    print("wav in max =", wav[0], wav[-1])
    new_norm = eniric.obsolete.snr_norm.normalize_flux(flux, id_string, new=True)
    old_norm = eniric.obsolete.snr_norm.normalize_flux(flux, id_string, new=False)

    print("new norm", new_norm)
    print("old_norm", old_norm)

    rvprec_new = Q.rv_precision(wav, new_norm)
    rvprec_old = Q.rv_precision(wav, old_norm)

    print("new rv=", rvprec_new, "old rv=", rvprec_old)
    assert np.abs(rvprec_new.value - rvprec_old.value) < 0.4


@pytest.mark.parametrize(
    "id_string",
    [
        "M0-1.0",
        "M3-1.0",
        "M6-1.0",
        "M9-1.0",
        "M0-5.0",
        "M3-5.0",
        "M6-5.0",
        "M9-5.0",
        "M0-10.0",
        "M3-10.0",
        "M6-10.0",
        "M9-10.0",
    ],
)
def test_snr_old_norm_constant(id_string):
    norm_const = eniric.obsolete.snr_norm.old_norm_constant(id_string)
    assert isinstance(norm_const, float)


@pytest.mark.parametrize(
    "bad_string", ["M0-1", "M0-2.5", "M8-1.0", "M6-5", "M9-10", "T0-3.0", "", "AB-CDE"]
)
def test_snr_old_norm_constant_with_bad_id_str(bad_string):
    """Fixed to the set of values in first paper."""
    with pytest.raises(ValueError):
        eniric.obsolete.snr_norm.old_norm_constant(bad_string)


@pytest.mark.parametrize(
    "id_string",
    [
        "M0-BAD-1.0-100k",
        "M9-A-5.0-50k",
        "MO-J-1.0-100k",
        "N0-J-1.0-100k",
        "M2--1.0-100k",
        "M0-J-2-100k",
        "M9-Z-5.0",
        "M0-J-1.0-100",
        "M0-J-1.0-1k",
        "M2-100k",
        "M0",
    ],
)
def test_decompose_bad_id_strings_give_errors(id_string):
    with pytest.raises(ValueError):
        eniric.obsolete.snr_norm.decompose_id_string(id_string)


@pytest.mark.parametrize(
    "id_string,expected",
    [
        ("M0-H-1.0-100k", ("M0", "H", "1.0", "100k")),
        ("M9-K-5.0-50k", ("M9", "K", "5.0", "50k")),
        ("M9-J-5.0-30k", ("M9", "J", "5.0", "30k")),
        ("M6-TEST-10.0-80k", ("M6", "TEST", "10.0", "80k")),
    ],
)
def test_decompose_id_string(id_string, expected):
    decomposed = eniric.obsolete.snr_norm.decompose_id_string(id_string)

    assert decomposed == expected
    assert len(decomposed) == 4


def test_old_normalization_does_not_handle_changed_band(resampled_data):
    id_string, wav, flux = resampled_data
    with pytest.raises(ValueError):
        eniric.obsolete.snr_norm.normalize_flux(
            flux, id_string, new=False, ref_band="K"
        )


def test_old_normalization_does_not_handle_changed_snr(resampled_data):
    id_string, wav, flux = resampled_data
    with pytest.raises(ValueError):
        eniric.obsolete.snr_norm.normalize_flux(flux, id_string, new=False, snr=101)


def test_get_ref_spectrum_with_ref_band_self(resampled_data):
    """Checks for upper or lower "self"."""
    id_string, wav, flux = resampled_data

    wav_ref, flux_ref = eniric.obsolete.snr_norm.get_reference_spectrum(
        id_string, ref_band="self"
    )

    # Reference is the same values
    assert np.allclose(wav, wav_ref)
    assert np.allclose(flux, flux_ref)


@pytest.mark.parametrize("ref_band", ["self", "SELF", "self", "SeLF"])
def test_get_self_band_can_be_any_case(resampled_data, ref_band):
    """Checks for upper or lower "self"."""

    id_string, wav, flux = resampled_data
    wav_ref, flux_ref = eniric.obsolete.snr_norm.get_reference_spectrum(
        id_string, ref_band=ref_band
    )

    # Reference is the same values
    assert np.allclose(wav, wav_ref)
    assert np.allclose(flux, flux_ref)


@pytest.mark.parametrize("consec_test", [True, False])
def test_barycenter_shift_verse_class(short_atmosphere, consec_test):
    """Test barycentric shift code is equivalent inside class.

    TODO: Delete this test when barycenter_shift function is removed."""
    atmos = short_atmosphere
    mask30kms = barycenter_shift(atmos.wl, atmos.mask, consecutive_test=consec_test)

    assert not np.allclose(mask30kms, atmos.mask) or (
        len(mask30kms) == np.sum(mask30kms) or np.sum(mask30kms) == 0
    )
    atmos.bary_shift_mask(consecutive_test=consec_test)
    # They are now close
    assert np.allclose(mask30kms, atmos.mask)


def test_read_spectrum():
    """Test reading in a _wave_photon.dat is the same as a _wave.dat."""
    photon = os.path.join(
        eniric.paths["test_data"],
        "obsolete",
        "sample_lte03900-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes_wave_photon.dat",
    )
    wave = os.path.join(
        eniric.paths["test_data"],
        "obsolete",
        "sample_lte03900-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes_wave.dat",
    )
    wave_wav, wave_flux = eniric.obsolete.utilities.read_spectrum(wave)
    photon_wav, photon_flux = eniric.obsolete.utilities.read_spectrum(photon)

    assert np.allclose(photon_wav, wave_wav)
    assert np.allclose(photon_flux, wave_flux)


@pytest.mark.parametrize(
    "filename",
    [
        os.path.join(
            eniric.paths["test_data"],
            "results",
            "Spectrum_M0-PHOENIX-ACES_Kband_vsini1.0_R100k.dat",
        ),
        os.path.join(
            eniric.paths["test_data"],
            "resampled",
            "Spectrum_M0-PHOENIX-ACES_Kband_vsini1.0_R100k_res3.0.dat",
        ),
    ],
)
def test_resampled_spectra_isnot_read_by_read_spectrum(filename):
    """Doesn't allow names with _vsini or _res in them."""
    with pytest.raises(ValueError, match="Using wrong function"):
        eniric.obsolete.utilities.read_spectrum(filename)


def test_band_snr_norm_test_data():
    """Compared to wav snr norm."""
    # snr_constant_band
    star, band, vel, res = "M0", "J", 1.0, "100k"
    test_data = os.path.join(
        eniric.paths["resampled"], resampled_template.format(star, band, vel, res)
    )
    wav, flux = io.pdread_2col(test_data)

    assert snrnorm.snr_constant_band(
        wav, flux, band="J", snr=100
    ) == snrnorm.snr_constant_wav(wav, flux, wav_ref=1.25, snr=100)

    assert snrnorm.snr_constant_band(
        wav, flux, band="J", snr=100
    ) != snrnorm.snr_constant_wav(wav, flux, wav_ref=1.24, snr=100)

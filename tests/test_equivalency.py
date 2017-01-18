
import numpy as np
import pytest
from hypothesis import given
from eniric.Qcalculator import RVprec_calc
from eniric.IOmodule import pdread_2col, pdread_3col, read_2col, read_3col

from eniric.original_code.Qcalculator import RVprec_calc as old_RVprec_calc
from eniric.original_code.IOmodule import read_2col as old_read_2col
from eniric.original_code.IOmodule import read_3col as old_read_3col

from eniric.IOmodule import pdwrite_2col, pdwrite_3col, write_e_2col, write_e_3col

import eniric.original_code.nIRanalysis as oldnIR
# import eniric.nIRanalysis as nIR
import eniric.utilities as eniric_utils
# To test the equivalence of code to check if it does the same thing:

# For python2.X compatibility
file_error_to_catch = getattr(__builtins__,'FileNotFoundError', IOError)



def test_convolution_indexing():
    """ To test equivalence of code to repalce for speed"""
    wav_extended = np.arange(100, 200)
    flux_extended = np.random.random(size=wav_extended.size)
    wav_val = 145
    R = 8
    FWHM = wav_val / R
    FWHM_lim = 5
    # Replace this code
    indexes = [i for i in range(len(wav_extended))
               if ((wav_val - FWHM_lim*FWHM) < wav_extended[i] <
               (wav_val + FWHM_lim*FWHM))]

    old_flux_2convolve = flux_extended[indexes[0]:indexes[-1]+1]
    old_wav_2convolve = wav_extended[indexes[0]:indexes[-1]+1]

    # Mask of wavelength range within 5 FWHM of wav
    index_mask = ((wav_extended > (wav_val - FWHM_lim*FWHM)) &
                  (wav_extended < (wav_val + FWHM_lim*FWHM)))

    new_flux_2convolve = flux_extended[index_mask]
    new_wav_2convolve = wav_extended[index_mask]

    assert np.all(old_flux_2convolve == new_flux_2convolve)
    assert np.all(old_wav_2convolve == new_wav_2convolve)


def test_rotational_convolution_indexing():
    wav_ext_rotation = np.arange(100, 200)
    flux_ext_rotation = np.random.random(size=wav_ext_rotation.size)
    wav = 145
    delta_lambda_L = 35

    # Old Code
    indexes = [i for i in range(len(wav_ext_rotation))
               if ((wav - delta_lambda_L) < wav_ext_rotation[i] <
               (wav + delta_lambda_L))]
    old_flux_2convolve = flux_ext_rotation[indexes[0]:indexes[-1]+1]
    old_wav_2convolve = wav_ext_rotation[indexes[0]:indexes[-1]+1]

    # New code
    index_mask = ((wav_ext_rotation > (wav - delta_lambda_L)) &
                  (wav_ext_rotation < (wav + delta_lambda_L)))

    new_flux_2convolve = flux_ext_rotation[index_mask]
    new_wav_2convolve = wav_ext_rotation[index_mask]

    assert np.all(old_flux_2convolve == new_flux_2convolve)
    assert np.all(old_wav_2convolve == new_wav_2convolve)


@pytest.mark.xfail(raises=file_error_to_catch)
def test_result_files_the_same():
    """ Test the result files are equal """

    print("Reading the data...")
    new_spectrum = "data/results/Spectrum_M0-PHOENIX-ACES_Yband_vsini1.0_R100k_unnormalized.txt"
    old_spectrum = "data/original_code/results/Spectrum_M0-PHOENIX-ACES_Yband_vsini1_R100k.txt"

    new_wavelength, new_ts, new_flux = pdread_3col(new_spectrum)
    old_wavelength, old_ts, old_spectrum = old_read_3col(old_spectrum)

    assert np.allclose(new_wavelength, np.array(old_wavelength, dtype="float64"))
    assert np.allclose(new_ts, np.array(old_ts, dtype="float64"))
    assert np.allclose(new_flux, np.array(old_spectrum, dtype="float64"))


@pytest.mark.xfail(raises=file_error_to_catch)    # Data file may not exist
def test_resampled_files_the_same():
    """ Test the resampled files are the same"""
    band = "GAP"
    new_spectrum = "data/resampled/Spectrum_M0-PHOENIX-ACES_{0}band_vsini1.0_R100k_unnormalized_res3.txt".format(band)
    old_spectrum = "data/original_code/resampled/original_code/results/Spectrum_M0-PHOENIX-ACES_{0}band_vsini1_R100k_res3.txt".format(band)
    new_wavelength, new_flux = pdread_2col(new_spectrum, noheader=True)
    old_wavelength, old_flux = old_read_2col(old_spectrum)
    assert np.allclose(old_wavelength, new_wavelength)
    assert np.allclose(old_flux, new_flux)


@pytest.mark.xfail(raises=file_error_to_catch)   # Data file may not exist
def test_resampled_RVprec_equal():
    """ Test quality of new and old spectra"""
    band = "GAP"
    new_spectrum = "data/resampled/Spectrum_M0-PHOENIX-ACES_{0}band_vsini1.0_R100k_unnormalized_res3.txt".format(band)
    old_spectrum = "data/original_code/resampled/original_code/results/Spectrum_M0-PHOENIX-ACES_{0}band_vsini1_R100k_res3.txt".format(band)

    new_wavelength, new_flux = pdread_2col(new_spectrum, noheader=True)

    old_wavelength, old_flux = old_read_2col(old_spectrum)

    # Scale with correct value

    new_RVprec = RVprec_calc(new_wavelength, new_flux)
    old_RVprec = old_RVprec_calc(old_wavelength, old_flux)
    assert np.allclose(new_RVprec.value, old_RVprec)
    assert new_RVprec.unit == "m / s"  # Check unit of precision


@pytest.mark.xfail(raises=file_error_to_catch)   # Data file may not exist
def test_list_creator():
    """ Test new masking in list creator is equivalent"""
    # test a couple of single bands only for speed

    # for band in ["H", "J", "K"]:
    for band in ["K"]:
        spectrum = "data/test_data/{0}_band_test_sample_lte03900-PHOENIX-ACES.dat".format(band)
        assert np.allclose(np.array(oldnIR.list_creator(spectrum, band)),
                           eniric_utils.list_creator(spectrum, band))


def test_pdread_2col():
    """ Test reading 2cols with pandas"""
    spectrum_1 = "data/test_data/Sample_input_phoenix.dat"
    spectrum_2 = "data/test_data/Sample_resampled_spectrum_res3.txt"

    wav_1_pd, flux_1_pd = pdread_2col(spectrum_1)
    wav_1, flux_1 = read_2col(spectrum_1)
    assert np.allclose(wav_1_pd, np.array(wav_1))
    assert np.allclose(flux_1_pd, np.array(flux_1))

    wav_2_pd, flux_2_pd = pdread_2col(spectrum_2, noheader=True)
    wav_2, flux_2 = read_2col(spectrum_2)
    assert np.allclose(wav_2_pd, np.array(wav_2))
    assert np.allclose(flux_2_pd, np.array(flux_2))


def test_pdread_3col():
    """ Test reading 3 cols with pandas.

    Use small sample file to reduce time for test.
    """
    filename = "data/test_data/Sample_results_spectrum.txt"

    wav_1_pd, theoretical_1_pd, flux_1_pd = pdread_3col(filename, noheader=True)
    wav_1, theoretical_1, flux_1 = read_3col(filename)
    assert np.allclose(wav_1_pd, np.array(wav_1))
    assert np.allclose(theoretical_1_pd, np.array(theoretical_1))
    assert np.allclose(flux_1_pd, np.array(flux_1))


def test_pdwriter():
    """ Check pd_writer same write_col with with exponential flag."""
    filedir = "data/test_data/"
    data = np.random.randn(3, 100) * 10000000
    pd2col_name = filedir + "pd2col_test.txt"
    pd3col_name = filedir + "pd3col_test.txt"
    twocol_name = filedir + "2col_test.txt"
    threecol_name = filedir + "3col_test.txt"

    # write files
    pdwrite_2col(pd2col_name, data[0], data[1])
    pdwrite_3col(pd3col_name, data[0], data[1],  data[2])
    write_e_2col(twocol_name, data[0], data[1])
    write_e_3col(threecol_name, data[0], data[1], data[2])

    # re-read files
    a = pdread_2col(pd2col_name)
    b = pdread_2col(twocol_name)
    c = pdread_3col(pd3col_name)
    d = pdread_3col(threecol_name)

    # check results the same
    assert np.allclose(a[0], b[0])
    assert np.allclose(a[1], b[1])
    assert np.allclose(c[0], d[0])
    assert np.allclose(c[1], d[1])
    assert np.allclose(c[2], d[2])

    # clean-up
    eniric_utils.silentremove(pd2col_name)
    eniric_utils.silentremove(pd3col_name)
    eniric_utils.silentremove(twocol_name)
    eniric_utils.silentremove(threecol_name)


def test_prepared_dat_files():
    """ Test that the flux inthe new prepared .dat files matches the original.
    This insures that all any conversions/scaling has been taken care of."""
    pass

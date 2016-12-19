
import numpy as np
import pandas as pd
import pytest
from hypothesis import given
from eniric.Qcalculator import RVprec_calc
from eniric.IOmodule import read_2col, read_3col

from eniric.original_code.Qcalculator import RVprec_calc as old_RVprec_calc
from eniric.original_code.IOmodule import read_2col as old_read_2col
from eniric.original_code.IOmodule import read_3col as old_read_3col

# To test the equivalence of code to check if it does the same thing:


def test_convolution_indexing():
    """ To test equivalence of code to repalce for speed"""
    wav_extended = np.arange(100, 200)
    flux_extended = np.random.random(size=wav_extended.size)
    wav_val = 145
    R = 8
    FWHM = wav_val/R
    FWHM_lim = 5
    # Replace this code
    indexes = [i for i in range(len(wav_extended))
               if ((wav_val - FWHM_lim*FWHM) < wav_extended[i] <
               (wav_val + FWHM_lim*FWHM))]

    old_flux_2convolve = flux_extended[indexes[0]:indexes[-1]+1]
    old_wav_2convolve = wav_extended[indexes[0]:indexes[-1]+1]

    # With this code

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


def test_result_files_the_same():
    """ Test the result files are equal """

    print("Reading the data...")
    new_spectrum = "data/results/Spectrum_M0-PHOENIX-ACES_Yband_vsini1_R100k.txt"
    old_spectrum = "data/original_code/results/Spectrum_M0-PHOENIX-ACES_Yband_vsini1_R100k.txt"

    data = pd.read_table(new_spectrum, delim_whitespace=True, names=["w", "ts", "s"], dtype=np.float64)
    new_wavelength, new_flux = data["w"].values, data['s'].values
    new_ts = data["ts"].values

    old_wavelength, old_ts, old_spectrum = old_read_3col(old_spectrum)

    assert np.allclose(new_wavelength, np.array(old_wavelength, dtype="float64"))
    assert np.allclose(new_ts, np.array(old_ts, dtype="float64"))
    assert np.allclose(new_flux, np.array(old_spectrum, dtype="float64"))


def test_resampled_files_the_same():
    """ Test the resampled files are the same"""
    new_spectrum = "data/resampled/results/Spectrum_M0-PHOENIX-ACES_Yband_vsini1_R100k_res3.txt"
    old_spectrum = "data/original_code/resampled/original_code/results/Spectrum_M0-PHOENIX-ACES_Yband_vsini1_R100k_res3.txt"
    new_wavelength, new_flux = read_2col(new_spectrum)
    old_wavelength, old_flux = old_read_2col(old_spectrum)
    assert np.all(old_wavelength == new_wavelength)
    assert np.all(old_flux == new_flux)


def test_resampled_RVprec_equal():
    """ Test quality of new and old spectra"""
    new_spectrum = "data/resampled/results/Spectrum_M0-PHOENIX-ACES_Yband_vsini1_R100k_res3.txt"
    old_spectrum = "data/original_code/resampled/original_code/results/Spectrum_M0-PHOENIX-ACES_Yband_vsini1_R100k_res3.txt"
    new_RVprec = RVprec_calc(spectrum_file=new_spectrum)
    old_RVprec = old_RVprec_calc(spectrum_file=old_spectrum)
    assert new_RVprec == old_RVprec


def test_RVprec_equal():
    """ Test quality of new and old spectra"""
    new_spectrum = "data/results/Spectrum_M0-PHOENIX-ACES_Yband_vsini1_R100k.txt"
    old_spectrum = "data/original_code/results/Spectrum_M0-PHOENIX-ACES_Yband_vsini1_R100k.txt"
    new_RVprec = RVprec_calc(spectrum_file=new_spectrum)
    old_RVprec = old_RVprec_calc(spectrum_file=old_spectrum)

    assert new_RVprec == old_RVprec

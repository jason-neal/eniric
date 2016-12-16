
import numpy as np
import pytest
from hypothesis import given
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
    indexes = [i for i in range(len(wav_extended)) if ((wav_val - FWHM_lim*FWHM) < wav_extended[i] < (wav_val + FWHM_lim*FWHM))]

    old_flux_2convolve = flux_extended[indexes[0]:indexes[-1]+1]
    old_wav_2convolve = wav_extended[indexes[0]:indexes[-1]+1]

    # With this code

    # Mask of wavelength range within 5 FWHM of wav
    index_mask = ((wav_extended > (wav_val - FWHM_lim*FWHM)) &
                  (wav_extended < (wav_val + FWHM_lim*FWHM)))

    new_flux_2convolve = flux_extended[index_mask]
    new_wav_2convolve = wav_extended[index_mask]

    assert  np.all(old_flux_2convolve == new_flux_2convolve)
    assert  np.all(old_wav_2convolve == new_wav_2convolve)

def test_rotational_convolution_indexing():
    wav_ext_rotation = np.arange(100, 200)
    flux_ext_rotation = np.random.random(size=wav_ext_rotation.size)
    wav = 145
    delta_lambda_L = 35

    # Old Code
    indexes = [i for i in range(len(wav_ext_rotation)) if ((wav - delta_lambda_L) < wav_ext_rotation[i] < (wav + delta_lambda_L))]
    old_flux_2convolve = flux_ext_rotation[indexes[0]:indexes[-1]+1]
    old_wav_2convolve = wav_ext_rotation[indexes[0]:indexes[-1]+1]

    # New code
    index_mask = ((wav_ext_rotation > (wav - delta_lambda_L)) &
          (wav_ext_rotation < (wav + delta_lambda_L)))

    new_flux_2convolve = flux_ext_rotation[index_mask]
    new_wav_2convolve = wav_ext_rotation[index_mask]

    assert np.all(old_flux_2convolve == new_flux_2convolve)
    assert np.all(old_wav_2convolve == new_wav_2convolve)

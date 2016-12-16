
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

    flux_2convolve = flux_extended[index_mask]
    wav_2convolve = wav_extended[index_mask]

    assert  np.all(old_flux_2convolve == flux_2convolve)
    assert  np.all(old_wav_2convolve == wav_2convolve)

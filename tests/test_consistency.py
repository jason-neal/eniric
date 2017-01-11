
# Consistency checks between original convolution code and the modified version

# Use pytest to run the same code and assert they are equal or close with np.allclose

# Decemeber 2016

from __future__ import division, print_function
import numpy as np
import pytest
from hypothesis import given, example
import hypothesis.strategies as st

import eniric.nIRanalysis as nir
import eniric.utilities as eniric_utils

import eniric.original_code.nIRanalysis as nir_org


@given(st.lists(st.floats(min_value=1e-7, max_value=1e-5, allow_infinity=False, allow_nan=False), unique=True, min_size=3, max_size=25), st.floats(min_value=1e-2, max_value=200), st.floats(min_value=1e-4, max_value=1))
def test_rotational_kernal(delta_lambdas, vsini, epsilon):
    """ Test that the new and original code produces the same output."""
    delta_lambdas = np.sort(np.asarray(delta_lambdas), kind='quicksort')
    delta_lambdas = np.append(np.flipud(delta_lambdas), np.insert(delta_lambdas, 0, 0))
    delta_lambda_L = np.max(delta_lambdas) * 2

    org_profile = nir_org.rotation_kernel(delta_lambdas, delta_lambda_L, vsini, epsilon)
    new_profile = eniric_utils.rotation_kernel(delta_lambdas, delta_lambda_L, vsini, epsilon)

    assert np.allclose(org_profile, new_profile)


@given(st.lists(st.lists(st.floats(min_value=1e-4, max_value=1, allow_infinity=False, allow_nan=False), min_size=2, max_size=2), min_size=1), st.floats(allow_infinity=False, allow_nan=False), st.floats(allow_infinity=False, allow_nan=False))
def test_wav_selector(wav_flux, wav_min, wav_max):
    """ Test that the wavelength selection code is equilivelent and works"""
    wav, flux = np.asarray(wav_flux).T[0], np.asarray(wav_flux).T[1]
    new_wav, new_flux = eniric_utils.wav_selector(wav, flux, wav_min, wav_max)
    org_wav, org_flux = nir_org.wav_selector(wav, flux, wav_min, wav_max)

    assert np.allclose(new_wav, org_wav)
    assert np.allclose(new_flux, org_flux)
    assert np.all(new_wav >= wav_min)
    assert np.all(new_wav <= wav_max)
    assert np.all(org_wav >= wav_min)
    assert np.all(org_wav <= wav_max)

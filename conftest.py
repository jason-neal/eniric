import os
from os.path import join

import astropy.units as u
import numpy as np
import pandas as pd
import pytest

import eniric.io_module as io
from eniric import config
from eniric.atmosphere import Atmosphere
from eniric.utilities import load_aces_spectrum


@pytest.fixture
def published_data():
    name = "data/precision/precision_figueira_2016.dat"
    df = pd.read_csv(name, sep="\t")
    return df


@pytest.fixture(
    params=[
        ("M0", "Z", 1, "60k"),
        ("M0", "J", 1, "100k"),
        ("M3", "Y", 5, "80k"),
        ("M6", "H", 10, "100k"),
        ("M9", "K", 1, "80k"),
    ]
)
def model_parameters(request):
    # Tuple of "SpType, band, vsini, R"
    return request.param


# Define some fixtures for eniric.precision.
per_s_cm2 = (1 / u.second) / (u.centimeter ** 2)


@pytest.fixture(
    params=[
        (np.arange(1, 101), np.random.random(100), None),
        (np.linspace(2.1, 2.5, 200), np.random.random(200), np.random.random(200)),
        (
            np.linspace(0.5, 1.5, 50),
            np.random.random(50),
            np.floor(2 * np.random.random(50)),
        ),
    ]
)
def test_spec(request):
    """Wave and flux, mask examples."""
    wav, flux, mask = request.param

    # Modify mask so that there are at least 2 consecutive zeros in "random" mask.
    if mask is not None:
        new_mask = np.ones_like(mask)
        new_mask[:-1] = new_mask[:-1] + mask[1:]  # Left add
        new_mask[1:] = new_mask[1:] + mask[:-1]  # Right add
        new_mask = new_mask > 2
        return wav, flux, new_mask.astype(int)
    else:
        return wav, flux, mask


# 3 situations each variable, no unit. a unit. or a dimensionless unscaled unit.
@pytest.fixture(params=[1, u.micron, u.dimensionless_unscaled])
def wav_unit(request):
    """Iterate some units on wavelength."""
    return request.param


@pytest.fixture(params=[1, per_s_cm2, u.dimensionless_unscaled])
def flux_unit(request):
    """Iterate some units on flux"""
    return request.param


@pytest.fixture(params=[1, u.dimensionless_unscaled])
def trans_unit(request):
    """Iterate some units on mask/transmission."""
    return request.param


@pytest.fixture(params=[True, False])
def grad_flag(request):
    """Gradient flag parameter."""
    return request.param


@pytest.fixture(
    params=[
        # "Average_TAPAS_2014_H.dat",
        "Average_TAPAS_2014_K.dat",
        "Average_TAPAS_2014_J.dat",
    ]
)
def atm_model(request):
    """Get atmospheric model name to load."""
    return join(config.pathdir, config.paths["atmmodel"], request.param)


@pytest.fixture(params=[2.5, 4])
def atmosphere_fixture(request, atm_model):
    percent_cutoff = request.param
    atm = Atmosphere.from_file(atm_model)
    atm.mask_transmission(percent_cutoff)
    return atm


@pytest.fixture()
def short_atmosphere(atmosphere_fixture):
    # First 2000 data points only to speed up tests
    return atmosphere_fixture[:2000]


@pytest.fixture(params=[(0, 1500), (8000, 9000)])
def sliced_atmmodel_default_mask(request, atm_model):
    """To do own masking. Sliced in different places."""
    lower, upper = request.param  # slice limits
    atm = Atmosphere.from_file(atm_model)
    return atm[int(lower) : int(upper)]


@pytest.fixture(params=[[3900, 4.5, 0, 0], [2600, 4.5, 0, 0]])
def testing_spectrum(request, use_test_config):
    wav, flux = load_aces_spectrum(request.param)
    return wav, flux


@pytest.fixture(scope="function")
def use_test_config():
    """Change configuration used to the test_config file."""
    original = config._path
    base_dir = os.path.dirname(__file__)
    test_filename = os.path.join(base_dir, "tests", "data", "test_config.yaml")
    config.change_file(test_filename)
    yield None
    config.change_file(original)

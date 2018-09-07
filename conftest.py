import os

import astropy.units as u
import numpy as np
import pandas as pd
import pytest

import eniric
import eniric.io_module as io
from eniric import utilities as utils
from eniric.atmosphere import Atmosphere

resampled_template = "Spectrum_{0}-PHOENIX-ACES_{1}band_vsini{2}_R{3}_res3.0.txt"


@pytest.fixture
def published_data():
    name = "data/precision/precision_figueira_2016.txt"
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


@pytest.fixture(
    params=[
        ("M0", "Z", 1.0, "60k"),
        ("M0", "K", 5.0, "60k"),
        ("M3", "Y", 5.0, "80k"),
        ("M6", "J", 10.0, "100k"),
        ("M6", "H", 10.0, "100k"),
        ("M9", "K", 1.0, "80k"),
    ]
)
def resampled_data(request):
    """Load a resampled spectra.

    Returns id-string, wavelength and flux.

    Fixture so that data files only get loaded once here.
    """
    star, band, vel, res = request.param
    id_string = "{0:s}-{1:s}-{2:.1f}-{3:s}".format(star, band, float(vel), res)

    test_data = os.path.join(
        eniric.paths["resampled"], resampled_template.format(star, band, vel, res)
    )
    wav, flux = io.pdread_2col(test_data)
    return id_string, wav, flux


# Define some fixtures for Qcalculator.
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
        # "Average_TAPAS_2014_H.txt",
        "Average_TAPAS_2014_K.txt",
        "Average_TAPAS_2014_J.txt",
    ]
)
def atm_model(request):
    """Get atmospheric model name to load."""
    return os.path.join(eniric.paths["atmmodel"], request.param)


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
def testing_spectrum(request):
    wav, flux = utils.load_aces_spectrum(request.param)
    return wav, flux

import os

import pandas as pd
import pytest

import eniric
from eniric.IOmodule import pdread_2col

resampled_template = "Spectrum_{0}-PHOENIX-ACES_{1}band_vsini{2}_R{3}_res3.0.txt"


@pytest.fixture
def published_data():
    name = "data/precision/precision_data_paper2015.txt"
    df = pd.read_csv(name, sep="\t")
    return df


@pytest.fixture(params=[
    ("M0", "Z", 1, "60k"),
    ("M0", "Z", 1, "100k"),
    ("M0", "Z", 10, "60k"),
    ("M0", "Z", 10, "100k"),
    ("M0", "H", 1, "60k"),
    ("M0", "H", 1, "100k"),
    ("M0", "H", 10, "60k"),
    ("M0", "H", 10, "100k"),
    ("M3", "Y", 10, "100k"),
    ("M3", "K", 10, "100k"),
    ("M0", "K", 1, "60k"),
    ("M0", "K", 1, "100k"),
    ("M0", "K", 5, "60k"),
    ("M0", "K", 5, "100k"),
    ("M6", "H", 1, "80k"),
    ("M6", "K", 1, "80k"),
    ("M9", "H", 1, "80k"),
    ("M9", "K", 1, "80k"),
])
def model_parameters(request):
    # Tuple of "SpType, band, vsini, R"
    return request.param


@pytest.fixture(params=[
    ("M0", "Z", 1.0, "60k"),
    ("M3", "Y", 10.0, "100k"),
    ("M6", "J", 10.0, "100k"),
    ("M9", "H", 1.0, "80k"),
    ("M0", "K", 5.0, "60k")
])
def resampled_data(request):
    """Load a resampled spectra.

    Returns id-string, wavelength and flux.

    Fixture so that data files only get loaded once here.
    """
    star, band, vel, res = request.param
    id_string = "{0:s}-{1:s}-{2:.1f}-{3:s}".format(star, band, float(vel), res)

    test_data = os.path.join(eniric.paths["resampled"],
                             resampled_template.format(star, band, vel, res))
    wav, flux = pdread_2col(test_data)
    return id_string, wav, flux

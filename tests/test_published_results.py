"""To test if the new code produces the same precision values on the published results."""

import numpy as np
import pandas as pd
import pytest

from eniric_scripts.nIR_precision import calculate_prec


@pytest.fixture
def published_data():
    name = "data/precision/precision_data_paper2015.txt"
    df = pd.read_csv(name, sep="\t")
    return df



@pytest.mark.parametrize("SpType,band,vsini,R,expected", [
    ("M0", "Y", 10.0, "100k", [20.7, 21.4, 20.9]),
    ("M9", "K", 5.0, "60k", [9.0, 10.1, 9.6]),
    ("M6", "H", 1.0, "80k", [4.1, 4.4, 4.2])
])
def test_old_calc_precision(SpType, band, vsini, R, expected):
    id_string = "{0:s}-{1:s}-{2:.1f}-{3:s}".format(SpType, band, float(vsini), R)
    results = calculate_prec([SpType], [band], [float(vsini)], [R], [3],
                             plot_atm=False, plot_ste=False,
                             plot_flux=False, paper_plots=False, rv_offset=0.0,
                             use_unshifted=False, snr=100, ref_band="J", new=False)

    # precision 1
    assert expected[0] == round(results[id_string][0].value, 1)

    # precision 3
    assert expected[2] == round(results[id_string][2].value, 1)


@pytest.mark.parametrize("SpType,band,vsini,R", [
    ("M0", "Z", 1, "60k"),
    ("M0", "H", 1, "60k"),
    ("M0", "Y", 10, "100k"),
    ("M0", "K", 5, "60k"),
    ("M0", "K", 5, "100k"),
    ("M6", "H", 1, "80k"),
    ("M9", "K", 5, "60k"),
    ("M9", "H", 1, "100k"),
    ("M6", "J", 10, "100k"),
    ("M3", "Y", 5, "80k")
])
def test_published_precision_with_old_normalization(SpType, band, vsini, R, published_data):
    # Testing calculate_prec with old normalization values are exactly the published values."""
    id_string = "{0:s}-{1:s}-{2:.1f}-{3:s}".format(SpType, band, float(vsini), R)
    results = calculate_prec([SpType], [band], [float(vsini)], [R], [3],
                             plot_atm=False, plot_ste=False,
                             plot_flux=False, paper_plots=False, rv_offset=0.0,
                             use_unshifted=False, snr=100, ref_band="J", new=False)
    print(published_data.head())

    # Precision 1
    published = published_data[published_data.Simulation == id_string]
    assert published["RV_Cond_1[m/s]"].values == round(results[id_string][0].value, 1)

    # precision 2 has changed
    assert published["RV_Cond_1[m/s]"].values != round(results[id_string][2].value, 1)

    # precision 3
    assert published["RV_Cond_3[m/s]"].values == round(results[id_string][2].value, 1)


@pytest.mark.parametrize("SpType,band,vsini,R", [
    ("M0", "Z", 1, "60k"),
    ("M0", "H", 1, "60k"),
    ("M0", "Y", 10, "100k"),
    ("M0", "K", 5, "60k"),
    ("M0", "K", 5, "100k"),
    ("M6", "H", 1, "80k"),
    ("M9", "K", 5, "60k"),
    ("M9", "H", 1, "100k"),
    ("M6", "J", 10, "100k"),
    ("M3", "Y", 5, "80k")
])
def test_published_precision_with_new_normalization(SpType, band, vsini, R, published_data):
    # Testing calculate_prec with new normalization values are similar."""
    id_string = "{0:s}-{1:s}-{2:.1f}-{3:s}".format(SpType, band, float(vsini), R)
    results = calculate_prec([SpType], [band], [float(vsini)], [R], [3],
                             plot_atm=False, plot_ste=False,
                             plot_flux=False, paper_plots=False, rv_offset=0.0,
                             use_unshifted=False, snr=100, ref_band="J", new=True)
    published = published_data[published_data.Simulation == id_string]
    print(published.head())

    # Precision 1
    published_prec_1 = published["RV_Cond_1[m/s]"].values
    result_1 = round(results[id_string][0].value, 1)
    assert np.abs(result_1 - published_prec_1) < 0.4
    assert (np.abs(result_1 - published_prec_1) / published_prec_1) < 0.025

    # precision 2 has changed
    published_prec_2 = published["RV_Cond_2[m/s]"].values
    result_2 = round(results[id_string][1].value, 1)
    assert np.abs(result_2 - published_prec_2) > 1
    # assert published["RV_Cond_1[m/s]"].values != round(results[id_string][2].value, 1)

    # precision 3
    published_prec_3 = published["RV_Cond_3[m/s]"].values
    result_3 = round(results[id_string][2].value, 1)
    assert np.abs(result_3 - published_prec_3) < 0.4
    assert (np.abs(result_3 - published_prec_3) / published_prec_3) < 0.025

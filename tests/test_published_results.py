"""To test if the new code produces the same precision values on the published results."""

import numpy as np
import pytest

from eniric_scripts.nIR_precision import calculate_prec


@pytest.mark.parametrize(
    "SpType,band,vsini,R,expected",
    [
        ("M3", "K", 10.0, "100k", [26.1, 29.6, 27.9]),
        ("M6", "H", 1.0, "80k", [4.1, 4.4, 4.2]),
        ("M0", "K", 5.0, "60k", [20.7, 23.9, 22.1]),
    ],
)
def test_old_calc_precision(SpType, band, vsini, R, expected):
    """Testing some direct examples from Figueira et al. 2016"""
    id_string = "{0:s}-{1:s}-{2:.1f}-{3:s}".format(SpType, band, float(vsini), R)
    results = calculate_prec(
        [SpType],
        [band],
        [float(vsini)],
        [R],
        [3],
        plot_atm=False,
        plot_ste=False,
        plot_flux=False,
        paper_plots=False,
        rv_offset=0.0,
        use_unshifted=False,
        snr=100,
        ref_band="J",
        new=False,
        grad=False,  # Old values without new gradient
    )

    for ii in [0, 2]:  # Condition 1 # Condition 3
        print("Condition # {}".format(ii + 1))
        assert expected[ii] == round(results[id_string][ii].value, 1)


def test_published_precision_with_old_normalization(model_parameters, published_data):
    # Testing calculate_prec with old normalization values are exactly the published values."""
    SpType, band, vsini, R = model_parameters
    id_string = "{0:s}-{1:s}-{2:.1f}-{3:s}".format(SpType, band, float(vsini), R)
    results = calculate_prec(
        [SpType],
        [band],
        [float(vsini)],
        [R],
        [3],
        plot_atm=False,
        plot_ste=False,
        plot_flux=False,
        paper_plots=False,
        rv_offset=0.0,
        use_unshifted=False,
        snr=100,
        ref_band="J",
        new=False,
        grad=False,  # Old values without new gradient
    )
    print(published_data.head())

    # Precision 1
    published = published_data[published_data.Simulation == id_string]
    assert published["RV_Cond_1[m/s]"].values == round(results[id_string][0].value, 1)

    # precision 2 has changed
    assert published["RV_Cond_2[m/s]"].values != round(results[id_string][1].value, 1)

    # precision 3
    assert published["RV_Cond_3[m/s]"].values == round(results[id_string][2].value, 1)


def test_published_precision_with_new_normalization(model_parameters, published_data):
    # Testing calculate_prec with new normalization values are similar."""
    SpType, band, vsini, R = model_parameters
    id_string = "{0:s}-{1:s}-{2:.1f}-{3:s}".format(SpType, band, float(vsini), R)
    results = calculate_prec(
        [SpType],
        [band],
        [float(vsini)],
        [R],
        [3],
        plot_atm=False,
        plot_ste=False,
        plot_flux=False,
        paper_plots=False,
        rv_offset=0.0,
        use_unshifted=False,
        snr=100,
        ref_band="J",
        new=True,
        grad=False,  # Old values without new gradient
    )
    published = published_data[published_data.Simulation == id_string]
    print(published.head())

    # Precision 1
    published_prec_1 = published["RV_Cond_1[m/s]"].values
    result_1 = round(results[id_string][0].value, 1)
    assert (np.abs(result_1 - published_prec_1) / published_prec_1) < 0.05

    # precision 2 has changed
    published_prec_2 = published["RV_Cond_2[m/s]"].values
    result_2 = round(results[id_string][1].value, 1)
    assert np.abs(result_2 - published_prec_2) > 1
    # assert published["RV_Cond_1[m/s]"].values != round(results[id_string][2].value, 1)

    # precision 3
    published_prec_3 = published["RV_Cond_3[m/s]"].values
    result_3 = round(results[id_string][2].value, 1)
    assert (np.abs(result_3 - published_prec_3) / published_prec_3) < 0.05


"""To test if the new code produces the same precision values on the published results."""

from __future__ import division, print_function

import numpy as np
import pytest

import eniric.IOmodule as io
import eniric.Qcalculator as Q
from bin.prec_1 import calc_prec1
from bin.nIR_precision import calculate_prec

# For python2.X compatibility
file_error_to_catch = getattr(__builtins__, 'FileNotFoundError', IOError)

path = "data/Published_Results/resampled/"


@pytest.mark.xfail(raises=file_error_to_catch)   # Data file may not exist
def test_precision_1():
    """New precision 1 test that works."""
    published_results = {1: 3.8, 5: 9.1, 10: 20.7}
    for vsini in [1, 5, 10]:
        # name = "Spectrum_M0-PHOENIX-ACES_Yband_vsini{0}.0_R100k_res3.txt".format(vsini)
        __, p1 = calc_prec1("M0", "Y", vsini, "100k", 3)

        assert abs(np.round(p1, 2).value - published_results[vsini]) < 0.5


@pytest.mark.parametrize("SpType,band,vsini,R,expected", [
    ("M0", "Y", 10.0, "100k", [20.7, 21.4, 20.9]),
    ("M9", "K", 5.0, "60k", [9.0, 10.1, 9.6]),
    ("M6", "H", 1.0, "80k", [4.1, 4.4, 4.2])
])
def test_old_calc_precision(SpType, band, vsini, R, expected):

    id = "{0:s}-{1:s}-{2:.1f}-{3:s}".format(SpType, band, float(vsini), R)
    results = calculate_prec([SpType], [band], [float(vsini)], [R], [3],
                             plot_atm=False, plot_ste=False,
                             plot_flux=False, paper_plots=False, rv_offset=0.0,
                             use_unshifted=False, snr=100, ref_band="J", new=False)

    # precision 1
    assert expected[0] == round(results[id][0].value, 1)

    # precision 3
    assert expected[2] == round(results[id][2].value, 1)


"""To test if the new code produces the same precision values on the published results."""

from __future__ import division, print_function

import numpy as np
import pytest

import eniric.IOmodule as io
import eniric.Qcalculator as Q
from bin.prec_1 import calc_prec1

# For python2.X compatibility
file_error_to_catch = getattr(__builtins__, 'FileNotFoundError', IOError)

path = "data/Published_Results/resampled/"


@pytest.mark.xfail(raises=file_error_to_catch)   # Data file may not exist
def test_presicion_1():
    """New precision 1 test that works."""
    published_results = {1: 3.8, 5: 9.1, 10: 20.7}
    path = "data/resampled/"
    for vsini in [1, 5, 10]:
        # name = "Spectrum_M0-PHOENIX-ACES_Yband_vsini{0}.0_R100k_res3.txt".format(vsini)
        __, p1 = calc_prec1("M0", "Y", vsini, "100k", 3, resampled_dir=path)

        # assert np.round(p1, 1).value == published_results[vsini]
        assert np.round(100 * p1, 1).value == published_results[vsini]  # With incorect normalization

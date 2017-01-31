
""" To test the equivalence of old and newcode to check if it does the same thing"""

from __future__ import division, print_function
import os
import numpy as np
import pytest
from hypothesis import given
import hypothesis.strategies as st

import eniric.IOmodule as IO
import eniric.utilities as utils

# For python2.X compatibility
file_error_to_catch = getattr(__builtins__, 'FileNotFoundError', IOError)


def test_pdread_2col():
    """ Test reading 2cols with pandas"""
    spectrum_1 = "data/test_data/Sample_input_phoenix.dat"
    spectrum_2 = "data/test_data/Sample_resampled_spectrum_res3.txt"

    wav_1_pd, flux_1_pd = IO.pdread_2col(spectrum_1)
    wav_1, flux_1 = IO.read_2col(spectrum_1)
    assert np.allclose(wav_1_pd, np.array(wav_1))
    assert np.allclose(flux_1_pd, np.array(flux_1))

    wav_2_pd, flux_2_pd = IO.pdread_2col(spectrum_2, noheader=True)
    wav_2, flux_2 = IO.read_2col(spectrum_2)
    assert np.allclose(wav_2_pd, np.array(wav_2))
    assert np.allclose(flux_2_pd, np.array(flux_2))


def test_pdread_3col():
    """ Test reading 3 cols with pandas.

    Use small sample file to reduce time for test.
    """
    filename = "data/test_data/Sample_results_spectrum.txt"

    wav_1_pd, theoretical_1_pd, flux_1_pd = IO.pdread_3col(filename, noheader=True)
    wav_1, theoretical_1, flux_1 = IO.read_3col(filename)
    assert np.allclose(wav_1_pd, np.array(wav_1))
    assert np.allclose(theoretical_1_pd, np.array(theoretical_1))
    assert np.allclose(flux_1_pd, np.array(flux_1))


def test_pdwriter():
    """ Check pd_writer same write_col with with exponential flag."""
    filedir = "data/test_data/"
    data = np.random.randn(3, 100) * 1e7
    pd2col_name = filedir + "pd2col_test.txt"
    pd3col_name = filedir + "pd3col_test.txt"
    twocol_name = filedir + "2col_test.txt"
    threecol_name = filedir + "3col_test.txt"

    # write files
    IO.pdwrite_2col(pd2col_name, data[0], data[1])
    IO.pdwrite_3col(pd3col_name, data[0], data[1],  data[2])
    IO.write_e_2col(twocol_name, data[0], data[1])
    IO.write_e_3col(threecol_name, data[0], data[1], data[2])

    # re-read files
    a = IO.pdread_2col(pd2col_name)
    b = IO.pdread_2col(twocol_name)
    c = IO.pdread_3col(pd3col_name)
    d = IO.pdread_3col(threecol_name)

    # check results the same
    assert np.allclose(a[0], b[0])
    assert np.allclose(a[1], b[1])
    assert np.allclose(c[0], d[0])
    assert np.allclose(c[1], d[1])
    assert np.allclose(c[2], d[2])

    # clean-up
    utils.silentremove(pd2col_name)
    utils.silentremove(pd3col_name)
    utils.silentremove(twocol_name)
    utils.silentremove(threecol_name)


def test_prepared_dat_files():
    """ Test that the flux inthe new prepared .dat files matches the original.
    This insures that all any conversions/scaling has been taken care of."""
    pass


def test_pdwrire_cols():
    """Test writer that can take variable column numbers"""
    filedir = "data/test_data"
    pd_multicol_name = os.path.join(filedir, "pd_multicol_test.txt")

    data1 = range(5)
    data2 = range(5, 10)
    bad_data = range(6)   # Different length

    # 0 means successful write
    assert 0 == IO.pdwrite_cols(pd_multicol_name, data1, data2, data1)
    assert 0 == IO.pdwrite_cols(pd_multicol_name, data1)
    assert 0 == IO.pdwrite_cols(pd_multicol_name, data1, data2, data1, data2,
                             header=["headers", "for", "column", "labels"])
    assert 0 == IO.pdwrite_cols(pd_multicol_name, data1, data2, sep=",", index=True)

    # test uneven dats lengths
    with pytest.raises(ValueError):
        IO.pdwrite_cols(pd_multicol_name, data1, data2, bad_data)

    # test bad header
    with pytest.raises(ValueError):
        IO.pdwrite_cols(pd_multicol_name, data1, data2, bad_data, header=["too", "many", "values"])

    with pytest.raises(TypeError):
        IO.pdwrite_cols(pd_multicol_name, data1, bad="keyword")

    # clean-up
    utils.silentremove(pd_multicol_name)

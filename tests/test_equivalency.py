
import os
import numpy as np
import pytest
from hypothesis import given

from eniric.Qcalculator import RVprec_calc as rvprec_calc
import eniric.IOmodule as IO

from eniric.IOmodule import pdread_2col, pdread_3col, read_2col, read_3col

from eniric.IOmodule import pdwrite_2col, pdwrite_3col, write_e_2col, write_e_3col

# import eniric.nIRanalysis as nIR
import eniric.utilities as eniric_utils
# To test the equivalence of code to check if it does the same thing:

# For python2.X compatibility
file_error_to_catch = getattr(__builtins__,'FileNotFoundError', IOError)


def test_pdread_2col():
    """ Test reading 2cols with pandas"""
    spectrum_1 = "data/test_data/Sample_input_phoenix.dat"
    spectrum_2 = "data/test_data/Sample_resampled_spectrum_res3.txt"

    wav_1_pd, flux_1_pd = pdread_2col(spectrum_1)
    wav_1, flux_1 = read_2col(spectrum_1)
    assert np.allclose(wav_1_pd, np.array(wav_1))
    assert np.allclose(flux_1_pd, np.array(flux_1))

    wav_2_pd, flux_2_pd = pdread_2col(spectrum_2, noheader=True)
    wav_2, flux_2 = read_2col(spectrum_2)
    assert np.allclose(wav_2_pd, np.array(wav_2))
    assert np.allclose(flux_2_pd, np.array(flux_2))


def test_pdread_3col():
    """ Test reading 3 cols with pandas.

    Use small sample file to reduce time for test.
    """
    filename = "data/test_data/Sample_results_spectrum.txt"

    wav_1_pd, theoretical_1_pd, flux_1_pd = pdread_3col(filename, noheader=True)
    wav_1, theoretical_1, flux_1 = read_3col(filename)
    assert np.allclose(wav_1_pd, np.array(wav_1))
    assert np.allclose(theoretical_1_pd, np.array(theoretical_1))
    assert np.allclose(flux_1_pd, np.array(flux_1))


def test_pdwriter():
    """ Check pd_writer same write_col with with exponential flag."""
    filedir = "data/test_data/"
    data = np.random.randn(3, 100) * 10000000
    pd2col_name = filedir + "pd2col_test.txt"
    pd3col_name = filedir + "pd3col_test.txt"
    twocol_name = filedir + "2col_test.txt"
    threecol_name = filedir + "3col_test.txt"

    # write files
    pdwrite_2col(pd2col_name, data[0], data[1])
    pdwrite_3col(pd3col_name, data[0], data[1],  data[2])
    write_e_2col(twocol_name, data[0], data[1])
    write_e_3col(threecol_name, data[0], data[1], data[2])

    # re-read files
    a = pdread_2col(pd2col_name)
    b = pdread_2col(twocol_name)
    c = pdread_3col(pd3col_name)
    d = pdread_3col(threecol_name)

    # check results the same
    assert np.allclose(a[0], b[0])
    assert np.allclose(a[1], b[1])
    assert np.allclose(c[0], d[0])
    assert np.allclose(c[1], d[1])
    assert np.allclose(c[2], d[2])

    # clean-up
    eniric_utils.silentremove(pd2col_name)
    eniric_utils.silentremove(pd3col_name)
    eniric_utils.silentremove(twocol_name)
    eniric_utils.silentremove(threecol_name)


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
    eniric_utils.silentremove(pd_multicol_name)

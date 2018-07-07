"""To test the equivalence of old and new code to check if it does the same thing."""
import os

import numpy as np
import pytest

import eniric.IOmodule as io
import eniric.utilities as utils


def test_pdread_2col():
    """Test reading 2cols with pandas."""
    spectrum_1 = "data/test_data/Sample_input_phoenix.dat"
    spectrum_2 = "data/test_data/Sample_resampled_spectrum_res3.txt"

    wav_1_pd, flux_1_pd = io.pdread_2col(spectrum_1)
    wav_1, flux_1 = io.read_2col(spectrum_1)
    assert np.allclose(wav_1_pd, np.array(wav_1))
    assert np.allclose(flux_1_pd, np.array(flux_1))

    wav_2_pd, flux_2_pd = io.pdread_2col(spectrum_2, noheader=True)
    wav_2, flux_2 = io.read_2col(spectrum_2)
    assert np.allclose(wav_2_pd, np.array(wav_2))
    assert np.allclose(flux_2_pd, np.array(flux_2))


def test_pdread_3col():
    """Test reading 3 cols with pandas.

    Use small sample file to reduce time for test.
    """
    filename = "data/test_data/Sample_results_spectrum.txt"

    wav_1_pd, theoretical_1_pd, flux_1_pd = io.pdread_3col(filename, noheader=True)
    wav_1, theoretical_1, flux_1 = io.read_3col(filename)
    assert np.allclose(wav_1_pd, np.array(wav_1))
    assert np.allclose(theoretical_1_pd, np.array(theoretical_1))
    assert np.allclose(flux_1_pd, np.array(flux_1))


def test_pdwriter():
    """Check pd_writer same write_col with with exponential flag."""
    filedir = "data/test_data/"
    data = np.random.randn(3, 100) * 1e7
    pd2col_name = filedir + "pd2col_test.txt"
    pd3col_name = filedir + "pd3col_test.txt"
    twocol_name = filedir + "2col_test.txt"
    threecol_name = filedir + "3col_test.txt"

    # write files
    io.pdwrite_2col(pd2col_name, data[0], data[1])
    io.pdwrite_3col(pd3col_name, data[0], data[1], data[2])
    io.write_e_2col(twocol_name, data[0], data[1])
    io.write_e_3col(threecol_name, data[0], data[1], data[2])

    # re-read files
    a = io.pdread_2col(pd2col_name)
    b = io.pdread_2col(twocol_name)
    c = io.pdread_3col(pd3col_name)
    d = io.pdread_3col(threecol_name)

    # check results the same
    assert np.allclose(a[0], b[0])
    assert np.allclose(a[1], b[1])
    assert np.allclose(c[0], d[0])
    assert np.allclose(c[1], d[1])
    assert np.allclose(c[2], d[2])

    # clean-up
    utils.silent_remove(pd2col_name)
    utils.silent_remove(pd3col_name)
    utils.silent_remove(twocol_name)
    utils.silent_remove(threecol_name)


def test_prepared_dat_files():
    """Test that the flux inthe new prepared .dat files matches the original.

    This insures that all any conversions/scaling has been taken care of.
    """
    pass


def test_pdwrire_cols():
    """Test writer that can take variable column numbers."""
    filedir = "data/test_data"
    pd_multicol_name = os.path.join(filedir, "pd_multicol_test.txt")

    data1 = range(5)
    data2 = range(5, 10)
    bad_data = range(6)  # Different length

    # 0 means successful write
    assert 0 == io.pdwrite_cols(pd_multicol_name, data1, data2, data1)
    assert 0 == io.pdwrite_cols(pd_multicol_name, data1)
    assert 0 == io.pdwrite_cols(
        pd_multicol_name,
        data1,
        data2,
        data1,
        data2,
        header=["headers", "for", "column", "labels"],
    )
    assert 0 == io.pdwrite_cols(pd_multicol_name, data1, data2, sep=",", index=True)

    # test uneven data lengths
    with pytest.raises(ValueError):
        io.pdwrite_cols(pd_multicol_name, data1, data2, bad_data)

    # test bad header
    with pytest.raises(ValueError):
        io.pdwrite_cols(
            pd_multicol_name, data1, data2, bad_data, header=["too", "many", "values"]
        )

    with pytest.raises(TypeError):
        io.pdwrite_cols(pd_multicol_name, data1, bad="keyword")

    # clean-up
    utils.silent_remove(pd_multicol_name)

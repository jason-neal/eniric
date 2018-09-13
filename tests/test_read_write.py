"""To test the equivalence of old and new code to check if it does the same thing."""

import numpy as np
import pytest

import eniric.io_module as io
import eniric.utilities as utils


def test_write_read_2col(tmpdir):
    """Check pd_writer same write_col with with exponential flag."""
    data = np.random.randn(2, 100) * 1e7

    tmp_file = tmpdir.join("file_col2.dat")
    tmp_file_name = str(tmp_file)
    # write files
    io.pdwrite_2col(tmp_file_name, data[0], data[1])

    # re-read files
    a = io.pdread_2col(tmp_file_name)

    # check results the same
    assert np.allclose(a[0], data[0])
    assert np.allclose(a[1], data[1])

    # clean-up
    utils.silent_remove(tmp_file_name)


def test_write_read_3col(tmpdir):
    """Check pd_writer same write_col with with exponential flag."""
    data = np.random.randn(3, 100) * 1e7

    tmp_file = tmpdir.join("file_col3.dat")
    tmp_file_name = str(tmp_file)

    # write files
    io.pdwrite_3col(tmp_file_name, data[0], data[1], data[2])

    # re-read files
    a = io.pdread_3col(tmp_file_name)

    # check results the same
    assert np.allclose(a[0], data[0])
    assert np.allclose(a[1], data[1])
    assert np.allclose(a[2], data[2])

    # clean-up
    utils.silent_remove(tmp_file_name)


def test_pdwrire_cols(tmpdir):
    """Test writer that can take variable column numbers."""
    tmpfile = tmpdir.join("pd_multicol_test.dat")
    pd_multicol_name = str(tmpfile)

    data1 = range(5)
    data2 = range(5, 10)

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
    # clean-up
    utils.silent_remove(pd_multicol_name)


def test_write_col_errors_different_length(tmpdir):
    tmpfile = tmpdir.join("pd_multicol_test.dat")
    pd_multicol_name = str(tmpfile)

    data1 = range(5)
    data2 = range(5, 10)
    bad_data = range(6)  # Different length

    # test uneven data lengths
    with pytest.raises(ValueError):
        io.pdwrite_cols(pd_multicol_name, data1, data2, bad_data)


def test_write_col_errors(tmpdir):
    tmpfile = tmpdir.join("pd_multicol_test.dat")
    pd_multicol_name = str(tmpfile)

    data1 = range(5)
    data2 = range(5, 10)
    bad_data = range(6)  # Different length

    # test bad header
    with pytest.raises(ValueError):
        io.pdwrite_cols(
            pd_multicol_name, data1, data2, bad_data, header=["too", "many", "values"]
        )


def test_write_col_errors_bad_keyword(tmpdir):
    tmpfile = tmpdir.join("pd_multicol_test.dat")
    pd_multicol_name = str(tmpfile)

    data1 = range(5)

    with pytest.raises(TypeError):
        io.pdwrite_cols(pd_multicol_name, data1, bad="keyword")

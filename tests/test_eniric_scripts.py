"""Tests for functions and bits of code from the scripts in 'scripts'."""

import argparse

import pytest

from scripts.split_atmmodel import check_positive


def test_check_positive():
    """Test that positive string values are returned as floats and negative values as errors."""
    assert check_positive("1") == 1.0
    assert isinstance(check_positive("20"), float)


@pytest.mark.parametrize(
    "val,error",
    [("-1", argparse.ArgumentTypeError), (-1.0, ValueError), (10, ValueError)],
)  # Not a float
def test_check_positive_errors(val, error):
    with pytest.raises(error):
        check_positive(val)

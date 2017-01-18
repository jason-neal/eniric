
""" Tests for functions and bits of code from the scripts in bin."""
import argparse
import pytest
from bin.split_atmmodel import check_positive


def test_check_positive():
    """ Test that positive values are returned as floats
    and negative values as errors."""

    assert check_positive("1") == 1.0
    assert check_positive(10) == 10.0

    with pytest.raises(argparse.ArgumentTypeError):
        check_positive("-1")
    with pytest.raises(argparse.ArgumentTypeError):
        check_positive(-1)
    with pytest.raises(argparse.ArgumentTypeError):
        check_positive(-0.1)

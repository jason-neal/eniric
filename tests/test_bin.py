"""Tests for functions and bits of code from the scripts in bin."""

from __future__ import division, print_function

import argparse

import pytest

from bin.split_atmmodel import check_positive
import eniric


def test_check_positive():
    """Test that positive string values are returned as floats and negative values as errors."""
    assert check_positive("1") == 1.0
    assert isinstance(check_positive("20"), float)

    with pytest.raises(ValueError):
        check_positive(10) == 10.0   # input must be float
    with pytest.raises(argparse.ArgumentTypeError):
        check_positive("-1")
    with pytest.raises(ValueError):
        check_positive(-1.0)

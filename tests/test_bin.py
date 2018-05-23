"""Tests for functions and bits of code from the scripts in bin."""

import argparse

import pytest

from bin.nIR_run import main as nir_run
from bin.split_atmmodel import check_positive


def test_check_positive():
    """Test that positive string values are returned as floats and negative values as errors."""
    assert check_positive("1") == 1.0
    assert isinstance(check_positive("20"), float)


@pytest.mark.parametrize("val,error", [
    ("-1", argparse.ArgumentTypeError),
    (-1.0, ValueError),
    (10, ValueError)])  # Not a float
def test_check_positive_errors(val, error):
    with pytest.raises(error):
        check_positive(val)


@pytest.mark.parametrize("noresample", [True, False])
@pytest.mark.parametrize("unnormalized", [True, False])
@pytest.mark.parametrize("org", [True, False])
def test_nir_run_raises_type_errors_on_non_lists(noresample, unnormalized, org):
    """Checking list inputs are needed.

    Checking over range of other boolean flags.
    """
    # Initialize parameters
    startype = ["M0", "M1"]
    vsini = [1]
    res = ["100k"]
    band = ["J", "K"]
    sample_rate = None

    # Check for TypeError on each parameter
    with pytest.raises(TypeError):
        nir_run(startype="M0", vsini=vsini, resolution=res, band=band, sample_rate=sample_rate,
                noresample=noresample, unnormalized=unnormalized, org=org)

    with pytest.raises(TypeError):
        nir_run(startype=startype, vsini=1.5, resolution=res, band=band, sample_rate=sample_rate,
                noresample=noresample, unnormalized=unnormalized, org=org)

    with pytest.raises(TypeError):
        nir_run(startype=startype, vsini=vsini, resolution="100k", band=band, sample_rate=sample_rate,
                noresample=noresample, unnormalized=unnormalized, org=org)

    with pytest.raises(TypeError):
        nir_run(startype=startype, vsini=vsini, resolution=res, band="K", sample_rate=sample_rate,
                noresample=noresample, unnormalized=unnormalized, org=org)

    with pytest.raises(TypeError):
        nir_run(startype=startype, vsini=vsini, resolution=res, band=band, sample_rate=5,
                noresample=noresample, unnormalized=unnormalized, org=org)

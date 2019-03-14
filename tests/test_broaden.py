import numpy as np
import pytest
from hypothesis import given, settings, strategies as st
from joblib import Parallel

from eniric.broaden import (
    convolution,
    resolution_convolution,
    rotation_kernel,
    rotational_convolution,
    unitary_gaussian,
)


@settings(max_examples=100)
@given(
    st.lists(
        st.floats(
            min_value=1e-7, max_value=1e-5, allow_infinity=False, allow_nan=False
        ),
        unique=True,
        min_size=3,
        max_size=25,
    ),
    st.floats(min_value=1e-2, max_value=200),
    st.floats(min_value=1e-4, max_value=1),
)
def test_rotational_kernel(delta_lambdas, vsini, epsilon):
    """Test that the new and original code produces the same output."""
    delta_lambdas = np.sort(np.asarray(delta_lambdas), kind="quicksort")
    delta_lambdas = np.append(np.flipud(delta_lambdas), np.insert(delta_lambdas, 0, 0))
    delta_lambda_l = np.max(delta_lambdas) * 2

    new_profile = rotation_kernel(delta_lambdas, delta_lambda_l, vsini, epsilon)

    assert len(new_profile) == len(delta_lambdas)


@given(
    st.lists(
        st.floats(min_value=-100, max_value=100, allow_nan=False),
        min_size=1,
        unique=True,
    ),
    st.floats(min_value=-100, max_value=100, allow_nan=False),
    st.floats(min_value=0.001, max_value=100, allow_nan=False),
)
def test_unitary_gaussian(x, center, fwhm):
    """Just a quick simple test."""
    x = np.asarray(x)

    gaussian = unitary_gaussian(x, center, fwhm)
    print(gaussian)
    # point at center should be the max
    assert len(gaussian) == len(x)
    assert np.allclose(np.max(gaussian), gaussian[np.argmin(abs(x - center))])


def test_unitary_gaussian_type_errors():
    """Testing for type errors."""
    x = np.arange(-10, 10)
    center = 0
    fwhm = 3

    with pytest.raises(TypeError):
        unitary_gaussian(x, center, "fwhm")
    with pytest.raises(TypeError):
        unitary_gaussian(x, "center", fwhm)
    with pytest.raises(TypeError):
        unitary_gaussian(range(-10, 10), "center", fwhm)
    with pytest.raises(TypeError):
        unitary_gaussian(1, "center", fwhm)


@pytest.mark.parametrize("num_proc", [1])
def test_convolution_can_accept_int(num_proc):
    n = 5
    x = np.linspace(2.0, 2.3, n)
    y = np.random.randn(n)
    convolution(x, y, vsini=1, R=1000, band="K", num_procs=num_proc)


@pytest.mark.parametrize("num_proc", [1])
def test_rot_convolution_can_accept_int(num_proc):
    n = 5
    x = np.linspace(2.0, 2.3, n)
    y = np.random.randn(n)
    rotational_convolution(x, x, y, vsini=1, num_procs=num_proc)


@pytest.mark.parametrize("num_proc", [1])
def test_res_convolution_can_accept_int(num_proc):
    n = 5
    x = np.linspace(2.0, 2.3, n)
    y = np.random.randn(n)
    resolution_convolution(x, x, y, R=1000, num_procs=num_proc)


@pytest.mark.parametrize("num_proc", [1, 2])
def test_convolution_can_accept_worker_pool(num_proc):
    n = 20
    x = np.linspace(2.0, 2.3, n)
    y = np.random.randn(n)
    with Parallel(num_proc) as parallel:
        convolution(
            x, y, vsini=1, R=100_000, band="K", num_procs=parallel, verbose=False
        )


@pytest.mark.parametrize("num_proc", [1, 2])
def test_rot_convolution_can_accept_worker_pool(num_proc):
    n = 20
    x = np.linspace(2.0, 2.3, n)
    y = np.random.randn(n)
    with Parallel(num_proc) as parallel:
        rotational_convolution(x, x, y, vsini=1, num_procs=parallel, verbose=False)


@pytest.mark.parametrize("num_proc", [1, 2])
def test_res_convolution_can_accept_worker_pool(num_proc):
    n = 20
    x = np.linspace(2.0, 2.3, n)
    y = np.random.randn(n)
    with Parallel(num_proc) as parallel:
        resolution_convolution(x, x, y, R=100_000, num_procs=parallel, verbose=False)


@pytest.mark.parametrize("num_proc", [3.14, "str"])
def test_rot_convolution_with_bad_num_proc(num_proc):
    n = 20
    x = np.linspace(2.0, 2.3, n)
    y = np.random.randn(n)
    with pytest.raises(
        TypeError, match="num_proc must be an int or joblib.parallel.Parallel"
    ):
        rotational_convolution(x, x, y, vsini=1, num_procs=num_proc, verbose=False)


@pytest.mark.parametrize("num_proc", [3.14, "str"])
def test_res_convolution_with_bad_num_proc(num_proc):
    n = 20
    x = np.linspace(2.0, 2.3, n)
    y = np.random.randn(n)
    with pytest.raises(
        TypeError, match="num_proc must be an int or joblib.parallel.Parallel"
    ):
        resolution_convolution(x, x, y, R=100_000, num_procs=num_proc, verbose=False)


def test_convolution_can_accept_None():
    n = 20
    x = np.linspace(2.0, 2.3, n)
    y = np.random.randn(n)
    convolution(x, y, vsini=1, R=100_000, band="K", num_procs=None)


def test_rot_convolution_can_accept_None():
    n = 20
    x = np.linspace(2.0, 2.3, n)
    y = np.random.randn(n)
    rotational_convolution(x, x, y, vsini=1, num_procs=None)


def test_res_convolution_can_accept_None():
    n = 20
    x = np.linspace(2.0, 2.3, n)
    y = np.random.randn(n)
    resolution_convolution(x, x, y, R=100_000, num_procs=None)


@pytest.mark.parametrize("zero_vsini", [0, 0.00, np.int(0.0), np.float64(0.0)])
def test_zero_rotation(zero_vsini):
    n = 20
    x = np.linspace(2.0, 2.3, n)
    y = np.random.randn(n)
    z = rotational_convolution(x, x, y, vsini=zero_vsini, num_procs=1)

    assert np.all(y == z)


@pytest.mark.parametrize("zero_vsini", [0, 0.00, np.int(0.0), np.float64(0.0)])
def test_zero_rotation_with_different_x(zero_vsini):
    n = 20
    x = np.linspace(1000, 2300, n)
    y = np.random.randn(n)
    z = rotational_convolution(x[4:16], x, y, vsini=zero_vsini, num_procs=1)

    assert np.all(y[4:16] == z)

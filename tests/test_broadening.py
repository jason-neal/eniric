import numpy as np
import pytest
from hypothesis import given, settings, strategies as st

from eniric.broaden import rotation_kernel, unitary_gaussian


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
    # other properties to test?


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

from datetime import datetime

import numpy as np
import pytest

from eniric.resample import log_resample, old_resample


@pytest.mark.parametrize("sampling", [1, 3, 5])
@pytest.mark.parametrize("resolution", [10, 1000, 50000, 100000])
@pytest.mark.parametrize("wave", [
    np.arange(2000, 2100),
    np.arange(100, 110, 0.1),
    np.arange(2.10, 2.15, 0.01)
])
def test_resamplers_equal(wave, sampling, resolution):
    grid_1 = old_resample(wave, sampling, resolution)
    grid_2 = log_resample(wave, sampling, resolution)
    print("len(grid_1)", len(grid_1))
    print("len(grid_1)", len(grid_2))
    assert np.allclose(grid_1, grid_2)


@pytest.mark.parametrize("sampling", [1.5, 3, 5])
@pytest.mark.parametrize("resolution", [10000, 50000, 100000])
@pytest.mark.parametrize("wave", [
    np.arange(2000, 2100),
    np.arange(100, 110, 0.1),
    np.arange(2.10, 2.15, 0.01)
])
def test_log_resamplers_faster(wave, sampling, resolution):
    """Test that the log resampling is faster."""
    t1 = datetime.now()
    old_resample(wave, sampling, resolution)
    t2 = datetime.now()
    log_resample(wave, sampling, resolution)
    t3 = datetime.now()
    print((t3 - t2))
    print((t2 - t1))
    assert (t3 - t2) < ((t2 - t1))  # * 0.1)  # 10 x Faster

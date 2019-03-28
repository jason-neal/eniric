import numpy as np
import pytest

from eniric.resample import log_resample, wl_logspace


@pytest.mark.parametrize("wav_start, wav_stop", [(2.1, 2.2), (1657.11, 1695)])
@pytest.mark.parametrize("sampling", [2, 5])
@pytest.mark.parametrize("resolution", [20000, 100_001])
def test_log_resample(wav_start, wav_stop, sampling, resolution):
    wav = np.linspace(wav_start, wav_stop, 50)

    resampled = log_resample(wav, sampling, resolution)

    assert resampled[0] == wav[0]
    assert resampled[-1] >= wav[-1]
    assert np.allclose(
        resampled[1:] / resampled[:-1],
        np.ones_like(resampled[:-1]) * (1 + 1 / (sampling * resolution)),
    )


@pytest.fixture(params=[(1, 1.1, 1.01), (1, 10, 1.25)])
def logspace_params(request):
    return request.param


def test_wl_logspace_endpoint(logspace_params):
    start, stop, base = logspace_params
    resample1 = wl_logspace(start, stop, base=base, end_point=False)
    resample2 = wl_logspace(start, stop, base=base, end_point=True)

    assert len(resample1) < len(resample2)

    assert resample1[-1] < stop
    assert resample2[-1] >= stop
    assert resample2[-1] < stop * base ** 2  # Not big


def test_wl_logspace_increases_by_bnd(logspace_params):
    start, stop, base = logspace_params
    resample1 = wl_logspace(start, stop, base=base, end_point=False)
    assert np.allclose(
        resample1[1:] / resample1[:-1], base * np.ones_like(resample1[:-1])
    )

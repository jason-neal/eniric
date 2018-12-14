import os

import pytest

from eniric import config

default_config = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'config.yaml')


class TestConfig:

    def test_default_filename(self):
        assert config.filename == default_config

    @pytest.mark.parametrize('key, value', [
        ('cache', {"location": ".joblib"}),
        ('Comments', 'Mid M dwarfs using emulator.\n'),
        ('chunk_ID', 0),
        ('spectrum_ID', 0),
        ('instrument_ID', 0)
    ])
    def test_base_keys(self, key, value):
        assert config[key] == value

    def test_name(self):
        assert config.name == 'default'

    def test_outdir(self):
        assert config.outdir == 'output/'

    def test_grid(self):
        assert isinstance(config.grid, dict)

    @pytest.mark.parametrize('key, value', [
        ('phoenix_raw', 'data/test_data/phoenix-raw'),
        ('btsettl_raw', 'data/test_data/btsettl-raw'),
        ('atmmodel', ['data/atmmodel']),
        ('precision', 'data/precsion'),
        ('test_data', "data/test_data"),
        ('wl_range', [6300, 6360]),
        ('buffer', 50)
    ])
    def test_paths_keys(self, key, value):
        assert config.paths[key] == value

    def test_paths(self):
        assert isinstance(config.paths, dict)

    def test_cache(self):
        assert isinstance(config.cache, dict)
        assert config.cache["location"] == ".joblib"

    def test_atmmodel(self):
        assert isinstance(config.atmmodel, dict)
        assert config.atmmodel["base"] == "Average_TAPAS_2014"

    def test_bands(self):
        assert isinstance(config.bands, dict)
        assert isinstance(config.bands["all"], list)
        assert "K" in config.bands["all"]
        assert "J" in config.bands["all"]

    def test_custom_bands(self):
        assert isinstance(config.custom_bands, dict)
        for key, value in config.custom_bands:
            assert isinstance(value, list)

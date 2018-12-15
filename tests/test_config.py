import os

import pytest
import yaml

from eniric import config

default_config = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'config.yaml')


class TestConfig:

    @pytest.mark.xfail()
    def test_default_filename(self):
        assert config.filename == default_config

    @pytest.mark.parametrize('key, value', [
        ('phoenix_raw', 'data/test_data/phoenix-raw'),
        ('btsettl_raw', 'data/test_data/btsettl-raw'),
        ('atmmodel', 'data/atmmodel'),
        ('precision', 'precision'),
        ('test_data', "data/test_data"),
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
        for value in config.custom_bands.values():
            assert isinstance(value, list)

    def test_lazy_load(self):
        previous = config.cache["location"]
        with open(config.filename, 'r+') as f:
            base = yaml.safe_load(f)
            base['cache'].update({'location': 'test_output'})
            yaml.dump(base, f)
        assert config.cache["location"] != previous
        assert config.cache["location"] == 'test_output'
        with open(config.filename, 'r+') as f:
            base = yaml.safe_load(f)
            base['cache'].update({"location" : previous})
            yaml.dump(base, f)
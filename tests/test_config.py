import os

import pytest
import yaml

from eniric import DEFAULT_CONFIG_FILE, config
from eniric._config import Config


class TestConfig:
    @pytest.fixture
    def test_config(self):
        """Config file for testing."""
        base_dir = os.path.dirname(__file__)
        filename = os.path.join(base_dir, "test_config.yaml")
        yield Config(filename)

    def test_default_filename(self):
        default_config = os.path.join(
            os.path.dirname(os.path.dirname(__file__)), "eniric", "config.yaml"
        )
        assert DEFAULT_CONFIG_FILE == default_config

    def test_base_dots(self):
        assert config.paths == config["paths"]
        assert config.cache == config["cache"]
        assert config.atmmodel == config["atmmodel"]
        assert config.bands == config["bands"]
        assert config.custom_bands == config["custom_bands"]

    @pytest.mark.parametrize(
        "key, value",
        [
            ("phoenix_raw", os.path.join(*["..", "data", "phoenix-raw"])),
            ("btsettl_raw", os.path.join(*["..", "data", "btsettl-raw"])),
            ("atmmodel", os.path.join(*["..", "data", "atmmodel"])),
            ("precision", "precision"),
            ("test_data", os.path.join(*["..", "data", "test_data"])),
        ],
    )
    def test_default_paths_keys(self, test_config, key, value):
        assert config.paths[key] == value

    @pytest.mark.parametrize(
        "key, value",
        [
            ("phoenix_raw", os.path.join(*["tests", "data", "phoenix-raw"])),
            ("btsettl_raw", os.path.join(*["tests", "data", "btsettl-raw"])),
            ("atmmodel", os.path.join(*["data", "atmmodel"])),
            ("precision", "precision"),
        ],
    )
    def test_paths_keys(self, test_config, key, value):
        assert test_config.paths[key] == value

    def test_paths(self):
        assert isinstance(config.paths, dict)

    def test_cache(self,):
        assert isinstance(config.cache, dict)
        assert config.cache["location"] == ".joblib"

    def test_atmmodel(self, test_config):
        assert isinstance(config.atmmodel, dict)
        assert config.atmmodel["base"] == "Average_TAPAS_2014"

    def test_default_bands(self):
        assert isinstance(config.bands, dict)
        assert isinstance(config.bands["all"], list)
        assert config.bands["all"] == [
            "VIS",
            "GAP",
            "Z",
            "Y",
            "J",
            "H",
            "K",
            "CONT",
            "NIR",
            "TEST",
        ]

    def test_bands(self, test_config):
        assert test_config.bands["all"] == ["K", "H", "J", "Y", "Z", "TEST"]

    def test_custom_bands(self):
        assert isinstance(config.custom_bands, dict)
        for value in config.custom_bands.values():
            assert isinstance(value, list)
            assert len(value) == 2

    def test_change_file(self):
        previous = config.name
        base_dir = os.path.dirname(__file__)
        filename = os.path.join(base_dir, "test_config.yaml")
        config.change_file(filename)
        assert config.name != previous
        assert config._path == filename

    def test_set_attr_fail_on_default(self):
        with pytest.raises(RuntimeError):
            config.name = "Stephen King"

    def test_set_base_attr(self, test_config):
        previous = test_config.name
        test_config.name = "new name"
        assert test_config.name == "new name"
        test_config.name = previous
        assert test_config.name == previous

    def test_set_non_base_attr(self, test_config):
        old_path = test_config.paths["btsettl_raw"]
        test_config.paths["btsettl_raw"] = "testpath_btsettl_raw"
        assert test_config.paths["btsettl_raw"] == "testpath_btsettl_raw"
        test_config.paths["btsettl_raw"] = old_path
        assert test_config.paths["btsettl_raw"] == old_path

    def test_copy_config(self, tmpdir):
        assert not os.path.exists(tmpdir.join("config.yaml"))
        config.copy_file(tmpdir)
        assert os.path.exists(tmpdir.join("config.yaml"))

    def test_lazy_load(self, test_config):
        previous = test_config.cache["location"]
        base = test_config._config
        base["cache"].update({"location": "test_output"})
        with open(test_config._path, "w") as f:
            yaml.safe_dump(base, f)

        assert test_config.cache["location"] == "test_output"
        test_config.cache["location"] = previous
        assert test_config.cache["location"] == previous

    def test_pathdir(self):
        assert config.pathdir == os.path.split(config._path)[0]
        with pytest.raises(AttributeError):
            config.pathdir = 5

    def test_pathdir_getter(self):
        assert config.pathdir == config.get_pathdir()

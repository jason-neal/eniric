import os

import pytest
import yaml

from eniric import DEFAULT_CONFIG_FILE, config
from eniric._config import Config

base_dir = os.path.dirname(__file__)
test_filename = os.path.join(base_dir, "data", "test_config.yaml")


class TestConfig:
    @pytest.fixture
    def test_config(self):
        """Config file for testing."""
        filename = test_filename
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
        "key, values",
        [
            ("phoenix_raw", ["..", "data", "phoenix-raw"]),
            ("btsettl_raw", ["..", "data", "btsettl-raw"]),
            ("atmmodel", ["..", "data", "atmmodel"]),
            ("precision_results", ["..", "data", "precision"]),
        ],
    )
    def test_default_paths_keys(self, key, values):
        assert config.paths[key] == os.path.join(*values)

    @pytest.mark.parametrize(
        "key, values",
        [
            ("phoenix_raw", ["phoenix-raw"]),
            ("btsettl_raw", ["btsettl-raw"]),
            ("atmmodel", ["..", "..", "data", "atmmodel"]),
            ("precision_results", ["..", "..", "data", "precision"]),
        ],
    )
    def test_paths_keys(self, test_config, key, values):
        assert test_config.paths[key] == os.path.join(*values)

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
        previous = config._path
        config.change_file(test_filename)
        assert config._path != previous
        assert config._path == test_filename
        config.change_file(previous)
        assert config._path == previous

    def test_set_attr_fail_on_default(self):
        with pytest.raises(RuntimeError):
            config.name = "Stephen King"

    def test_set_base_attr(self, test_config):
        previous_name = test_config.name
        test_config.name = "new name"
        assert test_config.name == "new name"
        test_config.name = previous_name
        assert test_config.name == previous_name

    def test_set_non_base_attr(self, test_config):
        old_path = test_config.paths["btsettl_raw"]
        test_config.paths["btsettl_raw"] = "testpath_btsettl_raw"
        assert test_config.paths["btsettl_raw"] == "testpath_btsettl_raw"
        test_config.paths["btsettl_raw"] = old_path
        assert test_config.paths["btsettl_raw"] == old_path

    @pytest.mark.parametrize("switch", [True, False])
    def test_copy_config(self, tmpdir, switch):
        previous_file = config._path
        assert not os.path.exists(tmpdir.join("config.yaml"))
        config.copy_file(tmpdir, switch=switch)
        if switch:
            assert str(tmpdir) in config.pathdir
        else:
            assert str(tmpdir) not in config.pathdir
        assert os.path.exists(tmpdir.join("config.yaml"))

        config.change_file(previous_file)  # Restore config

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

    def test_update_config_with_None(self):
        with pytest.raises(
            RuntimeError, match="The default file is not allowed to be overwritten."
        ):
            # Default config protection
            config.update(d=None)

    def test_update_test_config_with_None(self, test_config):
        previous_config = test_config
        test_config.update(d=None)
        assert previous_config == test_config

    def test_update_test_config_with_dict(self, test_config):
        temp_name = "test name"
        temp_atmmodel = {"base": "test_Average_TAPAS"}
        previous_name = test_config.name
        previous_atmmodel = test_config.atmmodel
        assert previous_name != temp_name
        assert previous_atmmodel != temp_atmmodel
        test_config.update(d={"name": temp_name, "atmmodel": temp_atmmodel})
        assert test_config.name == temp_name
        assert test_config.atmmodel == temp_atmmodel
        test_config.update(d={"name": previous_name, "atmmodel": previous_atmmodel})
        assert test_config.atmmodel == previous_atmmodel
        assert test_config.name == previous_name

    def test_update_with_kwargs(self, test_config):
        temp_name = "test name"
        temp_atmmodel = {"base": "test_Average_TAPAS"}
        previous_name = test_config.name
        previous_atmmodel = test_config.atmmodel
        assert previous_name != temp_name
        assert previous_atmmodel != temp_atmmodel
        test_config.update(d=None, name=temp_name, atmmodel=temp_atmmodel)
        assert test_config.name == temp_name
        assert test_config.atmmodel == temp_atmmodel
        test_config.update(d={"name": previous_name}, atmmodel=previous_atmmodel)
        assert test_config.atmmodel == previous_atmmodel
        assert test_config.name == previous_name

    @pytest.mark.parametrize("key", ["name", "cache", "bands"])
    def test_delitem(self, test_config, key):
        # Ability to delete from config.
        __ = test_config[key]
        del test_config[key]
        with pytest.raises(KeyError):
            __ = test_config[key]

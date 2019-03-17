# Inspired by Starfish

import os
import shutil

import oyaml as yaml

from . import DEFAULT_CONFIG_FILE


class Config(object):
    def __init__(self, path):
        """Creates a persistent Config object from the give file.

        This class is not meant to be instantiated by the User,
        but rather to be used via ``eniric.config``.

        Parameters
        ----------
        path: str or path-like
            The filename for creating the Config. Must be YAML.

        """
        self._path = path
        self._protect_rewrites = True

        with open(self._path, "r") as fd:
            self._config = yaml.safe_load(fd)

        self._protect_rewrites = os.path.abspath(path) == DEFAULT_CONFIG_FILE

    def get_pathdir(self):
        """Directory of the config file."""
        return os.path.dirname(self._path)

    pathdir = property(get_pathdir)

    def __contains__(self, key):
        return key in self._config

    def __delitem__(self, key):
        del self._config[key]

    def __eq__(self, other):
        return self._config.__eq__(other)

    def __getitem__(self, key):
        return self._config[key]

    def __setitem__(self, key, value):
        ret = self._config.__setitem__(key, value)
        return ret

    def __getattr__(self, key):
        if key in self:
            if key == "paths":
                paths = self["paths"]
                for k, value in paths.items():
                    if isinstance(value, list):
                        paths[k] = os.path.join(*value)
                return paths
            elif key == "cache":
                cache = self["cache"]
                if (cache["location"] is None) or (cache["location"] == "None"):
                    cache["location"] = None  # Disables caching
                elif isinstance(cache["location"], list):
                    cache["location"] = os.path.join(*cache["location"])
                return cache
            return self[key]
        else:
            return super().__getattribute__(key)

    def __setattr__(self, key, value):
        if key not in ["_path", "_protect_rewrites", "_config", "pathdir"]:
            if key in self:
                self.__setitem__(key, value)
                self._rewrite()

        super().__setattr__(key, value)

    def _rewrite(self):
        if self._protect_rewrites:
            raise RuntimeError(
                "The default file is not allowed to be overwritten. Please copy a file using "
                "config.copy_file(<path>) for your use."
            )
        with open(self._path, "w") as fd:
            yaml.safe_dump(self._config, fd)

    def update(self, d=None, **kwargs):
        """Update the config values and save to the config file."""
        if d is None:
            d = {}
        protected_rewrites = self._protect_rewrites
        self._protect_rewrites = True

        self._config.update(d, **kwargs)

        self._protect_rewrites = protected_rewrites

        self._rewrite()

    def change_file(self, filename):
        """Change the current configuration to use the given YAML file.

        Parameters
        ----------
        filename: str or path-like
            The YAML file to switch to using for config.

        Example
        -------
        .. code-block:: python

            eniric.config.change_file('new_config.yaml')

        """
        self._path = filename
        with open(self._path, "r") as fd:
            self._config = yaml.safe_load(fd)

    def copy_file(self, directory=None, switch=True):
        """Copies the master config file to the given directory.

        Parameters
        ----------
        directory: str or path-like
            The directory to copy the ``config.yaml`` file to. Default is os.getcwd().
        switch: bool
            If True, will switch the current config to use the copied file. Default is True

        Example
        -------
        .. code-block:: python

            eniric.config.copy_file()

        """
        if directory is None:
            directory = os.getcwd()
        outname = os.path.join(directory, "config.yaml")
        shutil.copy(DEFAULT_CONFIG_FILE, outname)
        if switch:
            self.change_file(outname)

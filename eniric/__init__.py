__version__ = "1.0rc1"

import os
import warnings

eniric_dir = os.path.dirname(__file__)
DEFAULT_CONFIG_FILE = os.path.join(eniric_dir, "config.yaml")

from ._config import Config

if (
    os.path.exists("config.yaml")
    and os.path.abspath("config.yaml") != DEFAULT_CONFIG_FILE
):
    config = Config("config.yaml")
else:
    warnings.warn(
        "Using the default config.yaml file located at {0}. This is likely NOT what you want and "
        "you will not be able to change any of the config values. Please use config.copy_file(<path>) to copy a "
        "version of the default config for your own project.".format(
            DEFAULT_CONFIG_FILE
        ),
        UserWarning,
    )
    config = Config(DEFAULT_CONFIG_FILE)

__all__ = [
    "atmosphere",
    "broaden",
    "config",
    "corrections",
    "legacy",
    "io_module",
    "precision",
    "resample",
    "snr_normalization",
    "utilities",
]

__version__ = "1.0rc1"

import os
import warnings

from ._config import Config

if os.path.exists("config.yaml"):
    config = Config("config.yaml")
else:
    base_dir = os.path.dirname(__file__)
    default = os.path.join(base_dir, "config.yaml")
    warnings.warn(
        "Using the default config.yaml file located at {0}. "
        "This is likely NOT what you want. Please create a similar "
        "'config.yaml' file in your current working directory.".format(default),
        UserWarning,
    )
    config = Config(default)

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

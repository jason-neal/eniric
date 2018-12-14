__version__ = "0.6-beta"
__all__ = [
    "atmosphere",
    "io_module",
    "Qcalculator",
    "resample",
    "snr_normalization",
    "utilities",
]


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


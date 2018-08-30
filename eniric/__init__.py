__version__ = "0.6-beta"
__all__ = [
    "atmosphere",
    "io_module",
    "Qcalculator",
    "resample",
    "snr_normalization",
    "utilities",
]

# Read the users config.yaml file.
# If it doesn't exist, print a useful help message

import os

import yaml

try:
    with open("config.yaml") as f:
        config = yaml.safe_load(f)
except FileNotFoundError as e:
    default = __file__[:-11] + "config.yaml"
    import warnings

    warnings.warn(
        "Using the default config.yaml file located at {0}. "
        "This is likely NOT what you want. Please create a similar "
        "'config.yaml' file in your current working directory.".format(default),
        UserWarning,
    )
    with open(default) as f:
        config = yaml.safe_load(f)

# Read the YAML variables into package-level dictionaries to be used by the other programs.
name = config["name"]
paths = config["paths"]
bands = config["bands"]
custom_bands = config["custom_bands"]
atmmodel = config["atmmodel"]
cache = config["cache"]

# Turn list into path
for key, value in paths.items():
    if isinstance(value, list):
        paths[key] = os.path.join(*value)

if (cache["location"] is None) or (cache["location"] == "None"):
    cache["location"] = None  # Disables caching
else:
    cache["location"] = os.path.join(*cache["location"])

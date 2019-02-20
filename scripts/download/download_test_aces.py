#!/usr/bin/env python
import os
from itertools import product

import numpy as np

import eniric
try:
    from Starfish.grid_tools import download_PHOENIX_models
except ImportError:
    raise ImportError("You need to install Starfish.")

models = np.array(list(product(
    (2600, 2800, 3500, 3900),
    (4.5,),
    (0.0,)
)))

outdir = eniric.config.paths["phoenix_raw"]

download_PHOENIX_models(models, outdir)

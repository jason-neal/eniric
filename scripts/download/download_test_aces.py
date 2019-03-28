#!/usr/bin/env python
from itertools import product
from os.path import join

import numpy as np

from eniric import config

try:
    from Starfish.grid_tools import download_PHOENIX_models
except ImportError:
    raise ImportError("You need to install Starfish.")


if __name__ == "__main__":
    models = np.array(list(product((2600, 2800, 3500, 3900), (4.5,), (0.0,))))

    outdir = join(config.pathdir, config.paths["phoenix_raw"])

    download_PHOENIX_models(models, outdir)

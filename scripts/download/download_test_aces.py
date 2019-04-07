#!/usr/bin/env python
import argparse
import sys
from itertools import product
from os.path import join

import numpy as np

from eniric import config

try:
    from Starfish.grid_tools import download_PHOENIX_models
except ImportError:
    raise ImportError("You need to install Starfish.")


def _parser():
    """Take care of all the argparse stuff."""
    parser = argparse.ArgumentParser(
        description="Download Phoenix ACES test spectra. A single spectra each for an M0, M3, M6, and M9 M-dwarf."
    )

    parser.add_argument(
        "-o",
        "--outdir",
        help="Output directory. Default is configured by the config.yaml.",
        type=str,
        default=join(config.pathdir, config.paths["phoenix_raw"]),
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = _parser()
    outdir = args.outdir

    models = np.array(list(product((2600, 2800, 3500, 3900), (4.5,), (0.0,))))
    download_PHOENIX_models(models, outdir)
    sys.exit(0)

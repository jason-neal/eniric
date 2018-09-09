#!/usr/bin/env python
"""To be run on travis to generate the results data.

Compare to published results and run other tests on it.
Don't do to many.
"""
from typing import List, Tuple

import eniric
from eniric.obsolete.nIR_run import main as nir_run
from eniric.obsolete.prepare_data import main as prepare_data

if __name__ == "__main__":
    print("Eniric paths: {}".format(eniric.paths))

    # subprocess.call("python eniric_scripts/prepare_data.py -s M0 M3 M6 M9 -l 4.50 -m 0.0", shell=True)
    prepare_data(
        startype=["M0", "M3", "M6", "M9"],
        temp=[],
        logg=[4.50],
        metallicity=[0],
        alpha=[0],
    )

    parameters = [
        (["M0"], ["Z", "J"], [1], ["60k", "100K"]),
        (["M0"], ["K", "J"], [5], ["60k"]),
        (["M3"], ["Y", "J"], [5], ["80k"]),
        (["M3"], ["K", "J"], [10], ["100k"]),
        (["M6"], ["J", "H"], [10], ["100k"]),
        (["M6"], ["J", "H"], [1], ["80k"]),
        (["M9"], ["J", "K"], [1], ["80k"]),
    ]  # type: List[Tuple[List[str], List[str], List[int], List[str]]]

    for (sptype, band, vel, res) in parameters:
        # subprocess.call(["python eniric_scripts/nIR_run.py
        #       -s {0} -b {1} -R {2} -v {3}".format(sptype, band, res, vel)], shell=True)

        nir_run(startype=sptype, vsini=vel, resolution=res, band=band)

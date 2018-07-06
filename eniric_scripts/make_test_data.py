#!/usr/bin/env python
"""To be run on travis to generate the results data.

Compare to published results and run other tests on it.
Don't do to many.
"""
# import subprocess
from datetime import datetime
from typing import List, Tuple
import eniric
from eniric_scripts.nIR_run import main as nir_run
from eniric_scripts.prepare_data import main as prepare_data

if __name__ == "__main__":
    print("Eniric paths: {}".format(eniric.paths))

    # subprocess.call("python eniric_scripts/prepare_data.py -s M0 M3 M6 M9 -l 4.50 -m 0.0", shell=True)
    prepare_data(
        startype=["M0", "M3", "M6", "M9"], temp=[], logg=[4.50], metallicity=[0], alpha=[0]
    )

    parameters = [
        (["M0"], ["Z", "H"], [1, 10], ["60k", "100k"]),
        (["M3"], ["Y", "K"], [10], ["100k"]),
        (["M0"], ["K"], [1, 5], ["60k", "100k"]),
        (["M6", "M9"], ["H", "K"], [1], ["80k"]),
        (["M0", "M3", "M6", "M9"], ["J"], [1, 5, 10], ["60k", "80k", "100k"]),
    ]  # type: List[Tuple[List[str], List[str], List[int], List[str]]]

    counter = 0
    start_time = datetime.now()
    for (sptype, band, vel, res) in parameters:
        # subprocess.call(["python eniric_scripts/nIR_run.py
        #       -s {0} -b {1} -R {2} -v {3}".format(sptype, band, res, vel)], shell=True)

        nir_run(startype=sptype, vsini=vel, resolution=res, band=band)
        counter += 1
    end_time = datetime.now()

    # List data files
    # print([d for d in os.walk("data")])

    print("Preformed {} convolutions in {}".format(counter, end_time - start_time))

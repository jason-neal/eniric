#!/usr/bin/env python
"""To be run on travis to generate the results data.

Compare to published results and run other tests on it.
Don't do to many.
"""
import os
import subprocess
from datetime import datetime
from eniric_scripts.nIR_run import main as nir_run
from eniric_scripts.prepare_data import main as prepare_data

# subprocess.call("python eniric_scripts/prepare_data.py -s M0 M3 M6 M9 -l 4.50 -m 0.0", shell=True)
prepare_data(startype=["M0", "M3", "M6", "M9"], temp=[], logg=[4.50], metallicity=[0], alpha=[0])

parameters = [("M0", "Z", 1, "60k"),
              ("M0", "H", 1, "60k"),
              ("M0", "Y", 10, "100k"),
              ("M0", "K", 5, ["60k", "100k"]),
              ("M6", "H", 1, "80k"),
              ("M9", "K", 5, "60k"),
              ("M9", "H", 1, "100k"),
              ("M6", "J", 10, "100k"),
              ("M3", "Y", 5, "80k")]

counter = 0
start_time = datetime.now()
for sptype, band, vel, res in parameters:
    # subprocess.call(["python eniric_scripts/nIR_run.py -s {0} -b {1} -R {2} -v {3}".format(sptype, band, res, vel)], shell=True)
    if not isinstance(res, list):
        res = [res]
    nir_run(startype=[sptype], vsini=[vel], resolution=res, band=[band])
    counter += 1
end_time = datetime.now()

# List data files
print([d for d in os.walk("data")])

print("Preformed {} convolutions in {}".format(counter, end_time-start_time))
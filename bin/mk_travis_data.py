import subprocess

"""To be run on travis to generate the results data.

Copare to published results and run other tests on it.

Don't do to many."""


subprocess.call("python bin/prepare_data.py -s M0 M3 M6 M9 -l 4.50 -m 0.0", shell=True)

parameters = [("M0", "Z", 1, "60k"),
              ("M0", "Y", 10, "100k"),
              ("M0", "K", 5, "60k 100k"),
              ("M6", "H", 1, "80k"),
              ("M9", "K", 5, "60k"),
              ("M9", "H", 1, "100k"),
              ("M6", "J", 10, "100k"),
              ("M3", "Y", 5, "80k")]

for sptype, band, vel, res in parameters:

    subprocess.call(["python bin/nIR_run.py -s {0} -b {1} -R {2} -v {3}".format(sptype, band, res, vel)], shell=True)

# List data files
subprocess.call(["ls -R data"], shell=True)

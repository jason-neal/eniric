"""This extracts the test data compressed fits.tar.gz files.

May not be needed if Starfish is able to handle .tar.gz files.
"""
import eniric
import tarfile
from os.path import exists, join


def unzip(dir: str, filename: str, overwrite: bool = False):
    joint_name = join(dir, filename)

    if (not exists(joint_name)) or (overwrite is True):
        with tarfile.open(join(joint_name + ".tar.gz"), "r") as tar:
            tar.extractall(dir)
    else:
        print("Did not overwrite file: {0}".format(joint_name))

if __name__ == "__main__":
    path = eniric.paths["phoenix_raw"]

    unzip(path, "WAVE_PHOENIX-ACES-AGSS-COND-2011.fits")

    aces_spec = [
        "lte02600-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits",
        "lte02800-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits",
        "lte03500-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits",
        "lte03900-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits",
    ]
    for aces in aces_spec:
        unzip(join(path, "Z-0.0"), aces)
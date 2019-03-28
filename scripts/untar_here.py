#!/usr/bin/env python

"""
untar_here.py
-------------
Bundled script to un-tar the eniric data downloaded.

Uses the tarfile module to extract the data.

"""
import sys
import tarfile

if __name__ == "__main__":
    filename = sys.argv[1]

    with tarfile.open(filename, "r") as tar:
        tar.extractall()

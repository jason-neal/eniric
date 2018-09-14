#!/usr/bin/env python

"""Bundled script to un-tar the eniric data downloaded."""
import sys
import tarfile

if __name__ == "__main__":
    filename = sys.argv[1]

    with tarfile.open(filename, "r") as tar:
        tar.extractall()

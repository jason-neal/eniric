#!/usr/bin/env python

"""
untar_here.py
-------------
Bundled script to un-tar the eniric data downloaded.

Uses the tarfile module to extract the data.

"""
import argparse
import sys
import tarfile


def _parser():
    """Take care of all the argparse stuff."""
    parser = argparse.ArgumentParser(description="Extract from a tar file.")

    parser.add_argument("filename", help="File to untar.", type=str, default="")
    return parser.parse_args()


if __name__ == "__main__":
    filename = _parser().filename

    with tarfile.open(filename, "r") as tar:
        tar.extractall()
    sys.exit(0)

#!/usr/bin/env python
"""
csv2tsv.py
----------
Simple script using pandas to transform from \*.csv to \*.tsv.

"""

import argparse
import sys

import pandas


def _parser():
    """Take care of all the argparse stuff."""
    parser = argparse.ArgumentParser(description="Transform a *.csv file into a *.tsv.")

    parser.add_argument("name", help="Filename to transform.", type=str, default="")
    return parser.parse_args()


if __name__ == "__main__":
    args = _parser()
    name = args.name

    df = pandas.read_csv(name, header=0)
    df.to_csv(name.replace(".csv", ".tsv"), sep="\t", header=True, index=False)
    sys.exit(0)

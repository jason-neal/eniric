#!/usr/bin/env python
"""
csv2tsv.py
----------
Simple script using pandas to transform from \*.csv to \*.tsv.

"""

import sys

import pandas

if __name__ == "__main__":
    name = sys.argv[1]

    df = pandas.read_csv(name, header=0)

    df.to_csv(name.replace(".csv", ".tsv"), sep="\t", header=True, index=False)

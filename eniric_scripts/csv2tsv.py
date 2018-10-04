#!/usr/bin/env python
import sys
import pandas

"""Simple script using pandas to transform from *.csv to *.tsv"""

if __name__ == "__main__":
    name = sys.argv[1]

    df = pandas.read_csv(name, header=0)

    df.to_csv(name.replace(".csv", ".tsv"), sep="\t", header=True, index=False)

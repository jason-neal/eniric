#!/usr/bin/env python
"""Script to process weekly TAPAS spectra to form an average atmospheric spectrum.

Created on Fri Feb  6 00:36:14 2015

@author: pfigueira
"""
import os
import string

import numpy as np

import eniric.IOmodule as io

# dirmodels = "/home/pfigueira/data/tapas/nIRanalysispaper/"
dirmodels = "/home/pfigueira/data/tapas/nIRanalysis_visible/"
list_files = "list_tapas_models.txt"

outdir = "atmmodel/"


def read_tapas(filename):
    """Read in tapas file."""
    file_conv = io.read_fullcol(filename)
    properties_dict = {}
    lambdas = []
    flux = []
    for line in file_conv:
        if line[0] == "|":
            # this is to ignore the header of the plot
            continue
        elif line[0] == "\\":
            # add everything starting with a \\ to the dictionary
            # remove the "\" form the beginning and the "\n" from the end
            key_dic, value_dic = string.split(line[1:-1], "=")
            properties_dict[key_dic] = value_dic
        else:
            # these are the values of the file itself
            lambda_act, flux_act = string.split(line)
            lambdas.append(float(lambda_act))
            flux.append(float(flux_act))

    lambdas = np.array(lambdas)
    flux = np.array(flux)
    # note that the lambdas are provided in reversed order
    return [properties_dict, lambdas[::-1], flux[::-1]]


def read_allfiles(mask_limit=0.02):
    """Read all the files in list_files."""
    files = io.read_fullcol(os.path.join(dirmodels, list_files))
    print("Reading the files...")
    atm_models = [
        read_tapas(os.path.join(dirmodels, dirmodels + file_act[:-1]))
        for file_act in files
    ]
    print("done.")

    wav = atm_models[0][1]
    mean_flux = []
    std_flux = []
    mask = []
    for i, __ in enumerate(wav):
        flux_at_wav = [model[2][i] for model in atm_models]
        mean_flux.append(np.average(flux_at_wav))
        std_flux.append(np.std(flux_at_wav))
        if mean_flux[-1] > (1.0 - mask_limit):
            # if transmission above threshold do not mask, otherwise mask
            mask.append(1.0)
        else:
            mask.append(0.0)

    write_4col_eeed(
        os.path.join(outdir, "Average_TAPAS_2014_visible.txt"),
        wav,
        mean_flux,
        std_flux,
        mask,
    )


def write_4col_eeed(filename, data1, data2, data3, data4):
    """Write data in 3 columns separated by tabs in a "filename" file."""
    f = open(filename, "w")

    for i, __ in enumerate(data1):
        # f.write("\t"+str(data1[i])+"\t\t"+str(data2[i])+"\t\t"+str(data3[i])+"\n")
        f.write(
            "\t{0:e}\t\t{1:e}\t\t{2:e}\t\t{3:d}\n".format(
                data1[i], data2[i], data3[i], data4[i]
            )
        )

    f.close()


################################################################################
#
#    Functions to read column-separated files
#
################################################################################
import pandas as pd
import numpy as np


def pdread_2col(filename, noheader=False):
    """ Read in a 2 column file with pandas.

    Faster then read_2col

    Parameters
    ----------
    filename: str
        Name of file to read.
    noheader: bool
        Flag indicating if there is no column names given in file.
        Default = False.

    Returns
    -------
    col1: ndarray
        First column as float64.
    col2: ndarray
        Second column as float64.
    """
    if noheader:
        data = pd.read_table(filename, comment='#', names=["col1", "col2"],
                         header=None, dtype=np.float64, delim_whitespace=True)
    else:
        data = pd.read_table(filename, comment='#', names=["col1", "col2"],
                         dtype=np.float64, delim_whitespace=True)


    return data["col1"].values, data["col2"].values


def pdread_3col(filename, noheader=False):
    """ Read in a 3 column file with pandas.

    Faster then read_3col

    Parameters
    ----------
    filename: str
        Name of file to read.
    noheader: bool
        Flag indicating if there is no column names given in file

    Returns
    -------
    col1: ndarray
        First column as float64.
    col2: ndarray
        Second column as float64.
    col3: ndarray
        Third column as float64.
    """
    if noheader:
        data = pd.read_table(filename, comment='#', names=["col1", "col2", "col3"],
                             header=None, dtype=np.float64, delim_whitespace=True)
    else:
        data = pd.read_table(filename, comment='#', names=["col1", "col2", "col3"],
                             dtype=np.float64, delim_whitespace=True)

    return data["col1"].values, data["col2"].values, data["col3"].values


def read_2col(filename):
    # The same as the previous, but returns 2 vectors, corresponding each
    # one to a column.THE RESULTS ARE FLOAT PYTHON VECTORS.
    # Note that in python all "float" are in fact "double-precision".

    list_data = read_col(filename)

    col1 = []
    col2 = []

    for i in range(len(list_data)):
        # checking if the line is valid
        if(list_data[i][0][0] != '#'):
            col1.append(float(list_data[i][0]))
            col2.append(float(list_data[i][1]))

    return [col1, col2]


def read_3col(filename):
    # The same as the previous, but returns 3 columns

    list_data = read_col(filename)

    col1 = []
    col2 = []
    col3 = []

    for i in range(len(list_data)):
        # checking if the line is valid
        if(list_data[i][0][0] != '#'):
            col1.append(float(list_data[i][0]))
            col2.append(float(list_data[i][1]))
            col3.append(float(list_data[i][2]))

    return [col1, col2, col3]


################################################################################
#
#    Functions to write files in column-separated formats
#
################################################################################















def write_2col(filename, data1, data2):
    # Writes data in 2 columns separated by tabs in a "filename" file.

    f = open(filename, "w")

    for i in range(len(data1)):
        f.write("\t"+str(data1[i])+"\t\t"+str(data2[i])+"\n")

    f.close()


def write_3col(filename, data1, data2, data3):
    # Writes data in 2 columns separated by tabs in a "filename" file.

    f = open(filename, "w")

    for i in range(len(data1)):
        f.write("\t"+str(data1[i])+"\t\t"+str(data2[i])+"\t\t"+str(data3[i])+"\n")

    f.close()



def write_e_2col(filename, data1, data2):
    # Writes data in 2 columns separated by tabs in a "filename" file.

    f = open(filename, "w")

    for i in range(len(data1)):
        # f.write("\t"+str(data1[i])+"\t\t"+str(data2[i])+"\t\t"+str(data3[i])+"\n")
        f.write("\t{0:e}\t\t{1:e}\n".format(data1[i], data2[i]))

    f.close()


def write_e_3col(filename, data1, data2, data3):
    # Writes data in 3 columns separated by tabs in a "filename" file.

    f = open(filename, "w")

    for i in range(len(data1)):
        # f.write("\t"+str(data1[i])+"\t\t"+str(data2[i])+"\t\t"+str(data3[i])+"\n")
        f.write("\t{0:e}\t\t{1:e}\t\t{2:e}\n".format(data1[i], data2[i], data3[i]))

    f.close()

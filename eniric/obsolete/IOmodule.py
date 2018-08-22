from typing import List

from numpy.core.multiarray import ndarray


def read_col(filename: str) -> List[List[str]]:
    """This program reads column formatted data from a file and
    returns a list in which each sublist correspond to the line's elements.
    THE RESULT IS A LIST OF STRINGS!"""

    f = open(filename, "r")

    list_data = []

    while 1:
        line = f.readline()

        if line == "":
            break
        if line[0] == "#":
            continue

        list_data.append(line.strip().split())

    f.close()

    return list_data


def read_2col(filename: str) -> List[List[float]]:
    """The same as the previous, but returns 2 vectors, corresponding each
    one to a column.THE RESULTS ARE FLOAT PYTHON VECTORS.
    Note that in python all "float" are in fact "double-precision"."""

    list_data = read_col(filename)

    col1 = []
    col2 = []

    for i, __ in enumerate(list_data):
        # checking if the line is valid
        if list_data[i][0][0] != "#":
            col1.append(float(list_data[i][0]))
            col2.append(float(list_data[i][1]))

    return [col1, col2]


def read_3col(filename: str) -> List[List[float]]:
    """The same as the previous, but returns 3 columns."""

    list_data = read_col(filename)

    col1 = []
    col2 = []
    col3 = []

    for i, __ in enumerate(list_data):
        # checking if the line is valid
        if list_data[i][0][0] != "#":
            col1.append(float(list_data[i][0]))
            col2.append(float(list_data[i][1]))
            col3.append(float(list_data[i][2]))

    return [col1, col2, col3]


def read_4col(filename: str) -> List[List[float]]:
    """The same as the previous, but returns 4 columns."""

    list_data = read_col(filename)

    col1 = []
    col2 = []
    col3 = []
    col4 = []

    for i, __ in enumerate(list_data):
        # checking if the line is valid
        if list_data[i][0][0] != "#":
            col1.append(float(list_data[i][0]))
            col2.append(float(list_data[i][1]))
            col3.append(float(list_data[i][2]))
            col4.append(float(list_data[i][3]))

    return [col1, col2, col3, col4]


################################################################################
#
#    Functions to write files in column-separated formats
#
################################################################################


def write_2col(filename, data1, data2):
    """Writes data in 2 columns separated by tabs in a "filename" file."""

    f = open(filename, "w")

    for i, __ in enumerate(data1):
        f.write("\t" + str(data1[i]) + "\t\t" + str(data2[i]) + "\n")

    f.close()


def write_3col(filename, data1, data2, data3):
    """Writes data in 2 columns separated by tabs in a "filename" file."""

    f = open(filename, "w")

    for i, __ in enumerate(data1):
        f.write(
            "\t"
            + str(data1[i])
            + "\t\t"
            + str(data2[i])
            + "\t\t"
            + str(data3[i])
            + "\n"
        )

    f.close()


def write_e_2col(filename: str, data1: ndarray, data2: ndarray) -> None:
    """Writes data in 2 columns separated by tabs in a "filename" file."""

    f = open(filename, "w")

    for i, __ in enumerate(data1):
        # f.write("\t"+str(data1[i])+"\t\t"+str(data2[i])+"\t\t"+str(data3[i])+"\n")
        f.write("\t{0:e}\t\t{1:e}\n".format(data1[i], data2[i]))

    f.close()


def write_e_3col(filename: str, data1: ndarray, data2: ndarray, data3: ndarray) -> None:
    """Writes data in 3 columns separated by tabs in a "filename" file."""

    f = open(filename, "w")

    for i, __ in enumerate(data1):
        # f.write("\t"+str(data1[i])+"\t\t"+str(data2[i])+"\t\t"+str(data3[i])+"\n")
        f.write("\t{0:e}\t\t{1:e}\t\t{2:e}\n".format(data1[i], data2[i], data3[i]))

    f.close()


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


def read_fullcol(filename):
    # This program reads column formatted data from a file and
    # returns a list in which each sublist correspond to the line's elements.
    # THE RESULT IS A LIST OF STRINGS!

    f = open(filename, "r")

    list_data = []

    while 1:
        line = f.readline()

        if line == "":
                    break
        if line[0] == '#':
                    continue

        list_data.append(line)

    f.close()

    return list_data


def read_col(filename):
    # This program reads column formatted data from a file and
    # returns a list in which each sublist correspond to the line's elements.
    # THE RESULT IS A LIST OF STRINGS!

    f = open(filename, "r")

    list_data = []

    while 1:
        line = f.readline()

        if line == "":
                    break
        if line[0] == '#':
                    continue

        list_data.append(line.strip().split())

    f.close()

    return list_data


def read_col_charsplit(filename, sepchar):
    # This program reads column formatted data from a file and
    # returns a list in which each sublist correspond to the line's elements separated by sepchar.
    # THE RESULT IS A LIST OF STRINGS!

    f = open(filename, "r")

    list_data = []

    while 1:
        line = f.readline()

        if line == "":
                    break
        if line[0] == '#':
                    continue

        list_data.append(line.strip().split(sepchar))

    f.close()

    return list_data


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


def read_2col1str(filename):
    # The same as the previous, but returns 2 columns and the first is a string

    list_data = read_col(filename)

    col1 = []
    col2 = []

    for i in range(len(list_data)):
        # checking if the line is valid
        if(list_data[i][0][0] != '#'):
            col1.append(list_data[i][0])
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


def read_3col1str(filename):
    # The same as the previous, but returns 3 columns and the first is a string

    list_data = read_col(filename)

    col1 = []
    col2 = []
    col3 = []

    for i in range(len(list_data)):
        # checking if the line is valid
        if(list_data[i][0][0] != '#'):
            col1.append(list_data[i][0])
            col2.append(float(list_data[i][1]))
            col3.append(float(list_data[i][2]))

    return [col1, col2, col3]


def read_3col2str(filename):
    # The same as the previous, but returns 3 columns and the first two are strings

    list_data = read_col(filename)

    col1 = []
    col2 = []
    col3 = []

    for i in range(len(list_data)):
        # checking if the line is valid
        if(list_data[i][0][0] != '#'):
            col1.append(list_data[i][0])
            col2.append(list_data[i][1])
            col3.append(float(list_data[i][2]))

    return [col1, col2, col3]


def read_4col(filename):
    # The same as the previous, but returns 4 columns

    list_data = read_col(filename)

    col1 = []
    col2 = []
    col3 = []
    col4 = []

    for i in range(len(list_data)):
        # checking if the line is valid
        if(list_data[i][0][0] != '#'):
            col1.append(float(list_data[i][0]))
            col2.append(float(list_data[i][1]))
            col3.append(float(list_data[i][2]))
            col4.append(float(list_data[i][3]))

    return [col1, col2, col3, col4]


def read_4col1str(filename):
    # The same as the previous, but returns 4 columns and the first is a string

    list_data = read_col(filename)

    col1 = []
    col2 = []
    col3 = []
    col4 = []

    for i in range(len(list_data)):
        # checking if the line is valid
        if(list_data[i][0][0] != '#'):
            col1.append(list_data[i][0])
            col2.append(float(list_data[i][1]))
            col3.append(float(list_data[i][2]))
            col4.append(float(list_data[i][3]))

    return [col1, col2, col3, col4]


def read_5col(filename):
    # The same as the previous, but returns 5 columns

    list_data = read_col(filename)

    col1 = []
    col2 = []
    col3 = []
    col4 = []
    col5 = []

    for i in range(len(list_data)):
                # checking if the line is valid
        if(list_data[i][0][0] != '#'):
            col1.append(float(list_data[i][0]))
            col2.append(float(list_data[i][1]))
            col3.append(float(list_data[i][2]))
            col4.append(float(list_data[i][3]))
            col5.append(float(list_data[i][4]))

    return [col1, col2, col3, col4, col5]


def read_5col1str(filename):
    # The same as the previous, but returns 5 columns of which one is a string

    list_data = read_col(filename)

    col1 = []
    col2 = []
    col3 = []
    col4 = []
    col5 = []

    for i in range(len(list_data)):
                # checking if the line is valid
        if(list_data[i][0][0] != '#'):
            col1.append(list_data[i][0])
            col2.append(float(list_data[i][1]))
            col3.append(float(list_data[i][2]))
            col4.append(float(list_data[i][3]))
            col5.append(float(list_data[i][4]))

    return [col1, col2, col3, col4, col5]


def read_6col1str(filename):
    # The same as the previous, but returns 6 columns and the first is a string

    list_data = read_col(filename)

    col1 = []
    col2 = []
    col3 = []
    col4 = []
    col5 = []
    col6 = []

    for i in range(len(list_data)):
                # checking if the line is valid
        if(list_data[i][0][0] != '#'):
            col1.append(list_data[i][0])
            col2.append(float(list_data[i][1]))
            col3.append(float(list_data[i][2]))
            col4.append(float(list_data[i][3]))
            col5.append(float(list_data[i][4]))
            col6.append(float(list_data[i][5]))

    return [col1, col2, col3, col4, col5, col6]


def read_6col(filename):
    # The same as the previous, but returns 6 columns

    list_data = read_col(filename)

    col1 = []
    col2 = []
    col3 = []
    col4 = []
    col5 = []
    col6 = []

    for i in range(len(list_data)):
                # checking if the line is valid
        if(list_data[i][0][0] != '#'):
            col1.append(float(list_data[i][0]))
            col2.append(float(list_data[i][1]))
            col3.append(float(list_data[i][2]))
            col4.append(float(list_data[i][3]))
            col5.append(float(list_data[i][4]))
            col6.append(float(list_data[i][5]))

    return [col1, col2, col3, col4, col5, col6]


def read_7col(filename):
    # The same as the previous, but returns 5 columns

    list_data = read_col(filename)

    col1 = []
    col2 = []
    col3 = []
    col4 = []
    col5 = []
    col6 = []
    col7 = []

    for i in range(len(list_data)):
                # checking if the line is valid
        if(list_data[i][0][0] != '#'):
            col1.append(float(list_data[i][0]))
            col2.append(float(list_data[i][1]))
            col3.append(float(list_data[i][2]))
            col4.append(float(list_data[i][3]))
            col5.append(float(list_data[i][4]))
            col6.append(float(list_data[i][5]))
            col7.append(float(list_data[i][6]))

    return [col1, col2, col3, col4, col5, col6, col7]


def read_7col1str(filename):
    # The same as the previous, but returns 7 columns being the first a string

    list_data = read_col(filename)

    col1 = []
    col2 = []
    col3 = []
    col4 = []
    col5 = []
    col6 = []
    col7 = []

    for i in range(len(list_data)):
                # checking if the line is valid
        if(list_data[i][0][0] != '#'):
            col1.append(list_data[i][0])
            col2.append(float(list_data[i][1]))
            col3.append(float(list_data[i][2]))
            col4.append(float(list_data[i][3]))
            col5.append(float(list_data[i][4]))
            col6.append(float(list_data[i][5]))
            col7.append(float(list_data[i][6]))

    return [col1, col2, col3, col4, col5, col6, col7]


def read_8col(filename):
    # The same as the previous, but returns 8 columns

    list_data = read_col(filename)

    col1 = []
    col2 = []
    col3 = []
    col4 = []
    col5 = []
    col6 = []
    col7 = []
    col8 = []

    for i in range(len(list_data)):
        # checking if the line is valid
        if(list_data[i][0][0] != '#'):
            col1.append(float(list_data[i][0]))
            col2.append(float(list_data[i][1]))
            col3.append(float(list_data[i][2]))
            col4.append(float(list_data[i][3]))
            col5.append(float(list_data[i][4]))
            col6.append(float(list_data[i][5]))
            col7.append(float(list_data[i][6]))
            col8.append(float(list_data[i][7]))

    return [col1, col2, col3, col4, col5, col6, col7, col8]


def read_8col1str(filename):
    # The same as the previous, but returns 8 columns the first being a string

    list_data = read_col(filename)

    col1 = []
    col2 = []
    col3 = []
    col4 = []
    col5 = []
    col6 = []
    col7 = []
    col8 = []

    for i in range(len(list_data)):
        # checking if the line is valid
        if(list_data[i][0][0] != '#'):
            col1.append(list_data[i][0])
            col2.append(float(list_data[i][1]))
            col3.append(float(list_data[i][2]))
            col4.append(float(list_data[i][3]))
            col5.append(float(list_data[i][4]))
            col6.append(float(list_data[i][5]))
            col7.append(float(list_data[i][6]))
            col8.append(float(list_data[i][7]))

    return [col1, col2, col3, col4, col5, col6, col7, col8]


def read_9col(filename):
    # The same as the previous, but returns 9 columns the first being a string

    list_data = read_col(filename)

    col1 = []
    col2 = []
    col3 = []
    col4 = []
    col5 = []
    col6 = []
    col7 = []
    col8 = []
    col9 = []

    for i in range(len(list_data)):
        # checking if the line is valid
        if(list_data[i][0][0] != '#'):
            col1.append(float(list_data[i][0]))
            col2.append(float(list_data[i][1]))
            col3.append(float(list_data[i][2]))
            col4.append(float(list_data[i][3]))
            col5.append(float(list_data[i][4]))
            col6.append(float(list_data[i][5]))
            col7.append(float(list_data[i][6]))
            col8.append(float(list_data[i][7]))
            col9.append(float(list_data[i][8]))

    return [col1, col2, col3, col4, col5, col6, col7, col8, col9]


def read_9col1str(filename):
    # The same as the previous, but returns 9 columns the first being a string

    list_data = read_col(filename)

    col1 = []
    col2 = []
    col3 = []
    col4 = []
    col5 = []
    col6 = []
    col7 = []
    col8 = []
    col9 = []

    for i in range(len(list_data)):
        # checking if the line is valid
        if(list_data[i][0][0] != '#'):
            col1.append(list_data[i][0])
            col2.append(float(list_data[i][1]))
            col3.append(float(list_data[i][2]))
            col4.append(float(list_data[i][3]))
            col5.append(float(list_data[i][4]))
            col6.append(float(list_data[i][5]))
            col7.append(float(list_data[i][6]))
            col8.append(float(list_data[i][7]))
            col9.append(float(list_data[i][8]))

    return [col1, col2, col3, col4, col5, col6, col7, col8, col9]


def read_10col(filename):
    # The same as the previous, but returns 11 columns

    list_data = read_col(filename)

    col1 = []
    col2 = []
    col3 = []
    col4 = []
    col5 = []
    col6 = []
    col7 = []
    col8 = []
    col9 = []
    col10 = []

    for i in range(len(list_data)):
        # checking if the line is valid
        if(list_data[i][0][0] != '#'):
            col1.append(float(list_data[i][0]))
            col2.append(float(list_data[i][1]))
            col3.append(float(list_data[i][2]))
            col4.append(float(list_data[i][3]))
            col5.append(float(list_data[i][4]))
            col6.append(float(list_data[i][5]))
            col7.append(float(list_data[i][6]))
            col8.append(float(list_data[i][7]))
            col9.append(float(list_data[i][8]))
            col10.append(float(list_data[i][9]))

    return [col1, col2, col3, col4, col5, col6, col7, col8, col9, col10]


def read_10col1strl(filename):
    # The same as the previous, but returns 10 columns with the last column being a string

    list_data = read_col(filename)

    col1 = []
    col2 = []
    col3 = []
    col4 = []
    col5 = []
    col6 = []
    col7 = []
    col8 = []
    col9 = []
    col10 = []

    for i in range(len(list_data)):
        # checking if the line is valid
        if(list_data[i][0][0] != '#'):
            col1.append(float(list_data[i][0]))
            col2.append(float(list_data[i][1]))
            col3.append(float(list_data[i][2]))
            col4.append(float(list_data[i][3]))
            col5.append(float(list_data[i][4]))
            col6.append(float(list_data[i][5]))
            col7.append(float(list_data[i][6]))
            col8.append(float(list_data[i][7]))
            col9.append(float(list_data[i][8]))
            col10.append(list_data[i][9])

    return [col1, col2, col3, col4, col5, col6, col7, col8, col9, col10]


def read_11col(filename):
    # The same as the previous, but returns 11 columns

    list_data = read_col(filename)

    col1 = []
    col2 = []
    col3 = []
    col4 = []
    col5 = []
    col6 = []
    col7 = []
    col8 = []
    col9 = []
    col10 = []
    col11 = []

    for i in range(len(list_data)):
                # checking if the line is valid
        if(list_data[i][0][0] != '#'):
            col1.append(float(list_data[i][0]))
            col2.append(float(list_data[i][1]))
            col3.append(float(list_data[i][2]))
            col4.append(float(list_data[i][3]))
            col5.append(float(list_data[i][4]))
            col6.append(float(list_data[i][5]))
            col7.append(float(list_data[i][6]))
            col8.append(float(list_data[i][7]))
            col9.append(float(list_data[i][8]))
            col10.append(float(list_data[i][9]))
            col11.append(float(list_data[i][10]))

    return [col1, col2, col3, col4, col5, col6, col7, col8, col9, col10, col11]


def read_11col1str(filename):
    # The same as the previous, but returns 11 columns the first being a string

    list_data = read_col(filename)

    col1 = []
    col2 = []
    col3 = []
    col4 = []
    col5 = []
    col6 = []
    col7 = []
    col8 = []
    col9 = []
    col10 = []
    col11 = []

    for i in range(len(list_data)):
                # checking if the line is valid
        if(list_data[i][0][0] != '#'):
            col1.append(list_data[i][0])
            col2.append(float(list_data[i][1]))
            col3.append(float(list_data[i][2]))
            col4.append(float(list_data[i][3]))
            col5.append(float(list_data[i][4]))
            col6.append(float(list_data[i][5]))
            col7.append(float(list_data[i][6]))
            col8.append(float(list_data[i][7]))
            col9.append(float(list_data[i][8]))
            col10.append(float(list_data[i][9]))
            col11.append(float(list_data[i][10]))

    return [col1, col2, col3, col4, col5, col6, col7, col8, col9, col10, col11]


def read_12col(filename):
    # The same as the previous, but returns 12 columns the first being a string

    list_data = read_col(filename)

    col1 = []
    col2 = []
    col3 = []
    col4 = []
    col5 = []
    col6 = []
    col7 = []
    col8 = []
    col9 = []
    col10 = []
    col11 = []
    col12 = []

    for i in range(len(list_data)):
                # checking if the line is valid
        if(list_data[i][0][0] != '#'):
            col1.append(float(list_data[i][0]))
            col2.append(float(list_data[i][1]))
            col3.append(float(list_data[i][2]))
            col4.append(float(list_data[i][3]))
            col5.append(float(list_data[i][4]))
            col6.append(float(list_data[i][5]))
            col7.append(float(list_data[i][6]))
            col8.append(float(list_data[i][7]))
            col9.append(float(list_data[i][8]))
            col10.append(float(list_data[i][9]))
            col11.append(float(list_data[i][10]))
            col12.append(float(list_data[i][11]))

    return [col1, col2, col3, col4, col5, col6, col7, col8, col9, col10, col11, col12]


def read_12col1str(filename):
    # The same as the previous, but returns 12 columns the first being a string

    list_data = read_col(filename)

    col1 = []
    col2 = []
    col3 = []
    col4 = []
    col5 = []
    col6 = []
    col7 = []
    col8 = []
    col9 = []
    col10 = []
    col11 = []
    col12 = []

    for i in range(len(list_data)):
                # checking if the line is valid
        if(list_data[i][0][0] != '#'):
            col1.append(list_data[i][0])
            col2.append(float(list_data[i][1]))
            col3.append(float(list_data[i][2]))
            col4.append(float(list_data[i][3]))
            col5.append(float(list_data[i][4]))
            col6.append(float(list_data[i][5]))
            col7.append(float(list_data[i][6]))
            col8.append(float(list_data[i][7]))
            col9.append(float(list_data[i][8]))
            col10.append(float(list_data[i][9]))
            col11.append(float(list_data[i][10]))
            col12.append(float(list_data[i][11]))

    return [col1, col2, col3, col4, col5, col6, col7, col8, col9, col10, col11, col12]


def read_2col_rdb(filename):
    # The same as the previous, but returns 2 columns
    # This is particularly usefull to read the "*_coralie.rdb" files

    list_data = read_col(filename)

    col1 = []
    col2 = []

    for i in range(2, len(list_data)):
                # checking if the line is valid
        if(list_data[i][0][0] != '#'):
            col1.append(float(list_data[i][0]))
            col2.append(float(list_data[i][1]))

    return [col1, col2]


def read_3col_rdb(filename):
    # The same as the previous, but returns 3 columns
    # This is particularly usefull to read the "*_coralie.rdb" files

    list_data = read_col(filename)

    col1 = []
    col2 = []
    col3 = []

    for i in range(2, len(list_data)):
                # checking if the line is valid
        if(list_data[i][0][0] != '#'):
            col1.append(float(list_data[i][0]))
            col2.append(float(list_data[i][1]))
            col3.append(float(list_data[i][2]))

    return [col1, col2, col3]


def read_4col_rdb(filename):
    # The same as the previous, but returns 6 columns
    # This is particularly usefull to read the "*_coralie.rdb" files

    list_data = read_col(filename)

    col1 = []
    col2 = []
    col3 = []
    col4 = []

    for i in range(2, len(list_data)):
                # checking if the line is valid
        if(list_data[i][0][0] != '#'):
            col1.append(float(list_data[i][0]))
            col2.append(float(list_data[i][1]))
            col3.append(float(list_data[i][2]))
            col4.append(float(list_data[i][3]))

    return [col1, col2, col3, col4]


def read_5col_rdb(filename):
    # The same as the previous, but returns 6 columns
    # This is particularly usefull to read the "*_coralie.rdb" files

    list_data = read_col(filename)

    col1 = []
    col2 = []
    col3 = []
    col4 = []
    col5 = []

    for i in range(2, len(list_data)):
                # checking if the line is valid
        if(list_data[i][0][0] != '#'):
            col1.append(float(list_data[i][0]))
            col2.append(float(list_data[i][1]))
            col3.append(float(list_data[i][2]))
            col4.append(float(list_data[i][3]))
            col5.append(float(list_data[i][4]))

    return [col1, col2, col3, col4, col5]


def read_6col_rdb(filename):
    # The same as the previous, but returns 6 columns
    # This is particularly usefull to read the "*_coralie.rdb" files

    list_data = read_col(filename)

    col1 = []
    col2 = []
    col3 = []
    col4 = []
    col5 = []
    col6 = []

    for i in range(2, len(list_data)):
                # checking if the line is valid
        if(list_data[i][0][0] != '#'):
            col1.append(float(list_data[i][0]))
            col2.append(float(list_data[i][1]))
            col3.append(float(list_data[i][2]))
            col4.append(float(list_data[i][3]))
            col5.append(float(list_data[i][4]))
            col6.append(float(list_data[i][5]))

    return [col1, col2, col3, col4, col5, col6]


def read_7col_rdb(filename):
    # The same as the previous, but returns 7 columns
    # This is particularly usefull to read the "*_coralie.rdb" files

    list_data = read_col(filename)

    col1 = []
    col2 = []
    col3 = []
    col4 = []
    col5 = []
    col6 = []
    col7 = []

    for i in range(2, len(list_data)):
                # checking if the line is valid
        if(list_data[i][0][0] != '#'):
            col1.append(float(list_data[i][0]))
            col2.append(float(list_data[i][1]))
            col3.append(float(list_data[i][2]))
            col4.append(float(list_data[i][3]))
            col5.append(float(list_data[i][4]))
            col6.append(float(list_data[i][5]))
            col7.append(float(list_data[i][6]))

    return [col1, col2, col3, col4, col5, col6, col7]


def read_10col_rdb(filename):
    # The same as the previous, but returns 10 columns
    # This is particularly usefull to read the "*_coralie.rdb" files

    list_data = read_col(filename)

    col1 = []
    col2 = []
    col3 = []
    col4 = []
    col5 = []
    col6 = []
    col7 = []
    col8 = []
    col9 = []
    col10 = []

    for i in range(2, len(list_data)):
                # checking if the line is valid
        if(list_data[i][0][0] != '#'):
            col1.append(float(list_data[i][0]))
            col2.append(float(list_data[i][1]))
            col3.append(float(list_data[i][2]))
            col4.append(float(list_data[i][3]))
            col5.append(float(list_data[i][4]))
            col6.append(float(list_data[i][5]))
            col7.append(float(list_data[i][6]))
            col8.append(float(list_data[i][7]))
            col9.append(float(list_data[i][8]))
            col10.append(float(list_data[i][9]))

    return [col1, col2, col3, col4, col5, col6, col7, col8, col9, col10]


def read_15col_rdb(filename):
    # The same as the previous, but returns 10 columns
    # This is particularly usefull to read the "*_coralie.rdb" files

    list_data = read_col(filename)

    col1 = []
    col2 = []
    col3 = []
    col4 = []
    col5 = []
    col6 = []
    col7 = []
    col8 = []
    col9 = []
    col10 = []
    col11 = []
    col12 = []
    col13 = []
    col14 = []
    col15 = []

    for i in range(2, len(list_data)):
                # checking if the line is valid
        if(list_data[i][0][0] != '#'):
            col1.append(float(list_data[i][0]))
            col2.append(float(list_data[i][1]))
            col3.append(float(list_data[i][2]))
            col4.append(float(list_data[i][3]))
            col5.append(float(list_data[i][4]))
            col6.append(float(list_data[i][5]))
            col7.append(float(list_data[i][6]))
            col8.append(float(list_data[i][7]))
            col9.append(float(list_data[i][8]))
            col10.append(float(list_data[i][9]))
            col11.append(float(list_data[i][10]))
            col12.append(float(list_data[i][11]))
            col13.append(float(list_data[i][12]))
            col14.append(float(list_data[i][13]))
            col15.append(float(list_data[i][14]))

    return [col1, col2, col3, col4, col5, col6, col7, col8, col9,
            col10, col11, col12, col13, col14, col15]


def read_Gauss():
    # The same as the previous, but returns 10 columns

    list_data = read_col(("/scisoft/i386/Packages/Python-2.4.3/myownpackages/GaussDist4py.txt"))
    Gaussdata = []

    for i in range(len(list_data)):
        for j in range(len(list_data[0])):
            Gaussdata.append(float(list_data[i][j]))

    return Gaussdata

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


def write_4col(filename, data1, data2, data3, data4):
    # Writes data in 2 columns separated by tabs in a "filename" file.

    f = open(filename, "w")

    for i in range(len(data1)):
        f.write("\t"+str(data1[i])+"\t\t"+str(data2[i])+"\t\t"+str(data3[i])+"\t\t"+str(data4[i])+"\n")

    f.close()


def write_5col(filename, data1, data2, data3, data4, data5):
    # Writes data in 5 columns separated by tabs in a "filename" file.

    f = open(filename, "w")

    for i in range(len(data1)):
        f.write("\t"+str(data1[i])+"\t\t"+str(data2[i])+"\t\t"+str(data3[i])+"\t\t"+str(data4[i])+"\t\t"+str(data5[i])+"\n")

    f.close()

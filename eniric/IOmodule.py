
################################################################################
#
#    Functions to read column-separated files
#
################################################################################
import pandas as pd
import numpy as np


def pdread_2col(filename, noheader=False):
    """Read in a 2 column file with pandas.

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
    try:
        if noheader:
            data = pd.read_table(filename, comment='#', names=["col1", "col2"],
                                 header=None, dtype=np.float64,
                                 delim_whitespace=True)
        else:
            data = pd.read_table(filename, comment='#', names=["col1", "col2"],
                                 dtype=np.float64, delim_whitespace=True)
    except:
        print("There was an error trying to read in the file \n{}".format(filename))
        raise

    return data["col1"].values, data["col2"].values


def pdread_3col(filename, noheader=False):
    """Read in a 3 column file with pandas.

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
    try:
        if noheader:
            data = pd.read_table(filename, comment='#', names=["col1", "col2", "col3"],
                                 header=None, dtype=np.float64, delim_whitespace=True)
        else:
            data = pd.read_table(filename, comment='#', names=["col1", "col2", "col3"],
                                 dtype=np.float64, delim_whitespace=True)
    except:
        print("There was an error trying to read in the file \n{}".format(filename))
        raise

    return data["col1"].values, data["col2"].values, data["col3"].values


def pdread_4col(filename, noheader=False):
    """Read in a 4 column file with pandas.

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
    col4: ndarray
        Fourth column as float64.
    """
    try:
        if noheader:
            data = pd.read_table(filename, comment='#', names=["col1", "col2", "col3", "col4"],
                                 header=None, dtype=np.float64, delim_whitespace=True)
        else:
            data = pd.read_table(filename, comment='#', names=["col1", "col2", "col3", "col4"],
                                 dtype=np.float64, delim_whitespace=True)
    except:
        print("There was an error trying to read in the file \n{}".format(filename))
        raise

    return data["col1"].values, data["col2"].values, data["col3"].values, data["col4"].values


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


def read_2col(filename):
    # The same as the previous, but returns 2 vectors, corresponding each
    # one to a column.THE RESULTS ARE FLOAT PYTHON VECTORS.
    # Note that in python all "float" are in fact "double-precision".

    list_data = read_col(filename)

    col1 = []
    col2 = []

    for i, __ in enumerate(list_data):
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

    for i, __ in enumerate(list_data):
        # checking if the line is valid
        if(list_data[i][0][0] != '#'):
            col1.append(float(list_data[i][0]))
            col2.append(float(list_data[i][1]))
            col3.append(float(list_data[i][2]))

    return [col1, col2, col3]


def read_4col(filename):
    # The same as the previous, but returns 4 columns

    list_data = read_col(filename)

    col1 = []
    col2 = []
    col3 = []
    col4 = []

    for i, __ in enumerate(list_data):
        # checking if the line is valid
        if(list_data[i][0][0] != '#'):
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


def pdwrite_2col(filename, data1, data2, sep="\t", header=False, float_format=None):
    """Write out a 2 column file with pandas.

    Faster then write_2col, uses pandas.DataFrame.to_csv()

    Parameters
    ----------
    filename: str
        Name of file to write.
    data1: ndarray or list, array-like
        The data for the first column
    data2: ndarray or list, array-like
        The data for the second column
    sep: str
        Character separation between values.
    header: list of strings or bool, default False
        Header strings to apply to columns.
    float_format: str default None
        Specify floating point string format.

    Returns
    -------
    flag: bool
        Returns 0 if successful.
    """
    if header:
        df = pd.DataFrame({"# " + header[0]: data1, header[1]: data2})
    else:
        df = pd.DataFrame({"# x": data1, "y": data2})

    # Write dataframe to file
    df.to_csv(filename, sep=sep, header=header, index=False, float_format=float_format)  # header=False

    return 0


def pdwrite_3col(filename, data1, data2, data3, sep="\t", header=False, float_format=None):
    """Write out a 3 column file with pandas.

    Faster then write_3col, uses pandas.DataFrame.to_csv()

    Parameters
    ----------
    filename: str
        Name of file to write.
    data1: ndarray or list, array-like
        The data for the first column
    data2: ndarray or list, array-like
        The data for the second column
    data3: ndarray or list, array-like
        The data for the third column
    sep: str
        Character separation between values.
    header: list of strings or bool, default False
        Header strings to apply to columns.
    float_format: str default None
        Specify floating point string format.

    Returns
    -------
    flag: bool
        Returns 0 if successful.
    """
    if header:
        df = pd.DataFrame({"# " + header[0]: data1, header[1]: data2, header[2]: data3})
    else:
        df = pd.DataFrame({"# x": data1, "y": data2, "z": data3})

    # Write dataframe to file
    df.to_csv(filename, sep=sep, header=header, index=False, float_format=float_format)  # header=False

    return 0


def write_2col(filename, data1, data2):
    # Writes data in 2 columns separated by tabs in a "filename" file.

    f = open(filename, "w")

    for i, __ in enumerate(data1):
        f.write("\t"+str(data1[i])+"\t\t"+str(data2[i])+"\n")

    f.close()


def write_3col(filename, data1, data2, data3):
    # Writes data in 2 columns separated by tabs in a "filename" file.

    f = open(filename, "w")

    for i, __ in enumerate(data1):
        f.write("\t"+str(data1[i])+"\t\t"+str(data2[i])+"\t\t"+str(data3[i])+"\n")

    f.close()


def write_e_2col(filename, data1, data2):
    # Writes data in 2 columns separated by tabs in a "filename" file.

    f = open(filename, "w")

    for i, __ in enumerate(data1):
        # f.write("\t"+str(data1[i])+"\t\t"+str(data2[i])+"\t\t"+str(data3[i])+"\n")
        f.write("\t{0:e}\t\t{1:e}\n".format(data1[i], data2[i]))

    f.close()


def write_e_3col(filename, data1, data2, data3):
    # Writes data in 3 columns separated by tabs in a "filename" file.

    f = open(filename, "w")

    for i, __ in enumerate(data1):
        # f.write("\t"+str(data1[i])+"\t\t"+str(data2[i])+"\t\t"+str(data3[i])+"\n")
        f.write("\t{0:e}\t\t{1:e}\t\t{2:e}\n".format(data1[i], data2[i], data3[i]))

    f.close()


def pdwrite_cols(filename, *data, **kwargs):
    """Write out a csv file with pandas, variable columns possible.

    Uses pandas.DataFrame.to_csv()

    Parameters
    ----------
    filename: str
        Name of file to write.
    *data: ndarray or list, array-like
        Variable number of data columns to be writen in the given order.
    **kwargs: dict
        Keyword args for pandas
    sep: str, default="\t"
        Character separation between values.
    header: list of strings or bool, default False
        Header strings to apply to columns. Must be equal to number
        of data columns provided.

    Returns
    -------
    flag: bool
        Returns 0 if successful.
    """

    # unpack keyword args, second argument is the defualt if not found.
    header = kwargs.pop('header', False)
    sep = kwargs.pop('sep', "\t")
    index = kwargs.pop('index', False)
    float_format = kwargs.pop('float_format', '%.6f')
    # TODO: See about passing any extra keywords into pandas call
    if kwargs:   # check for unwanted kewords
        raise TypeError('Unexpected **kwargs: {!r}'.format(kwargs))

    if header:
        if len(header) != len(data):
            raise ValueError("Size of data and header does not match.")

    data_dict = {}
    for i, data_i in enumerate(data):
        data_dict[i] = data[i]    # keys are assigned the index value from enumerate

        if len(data[i]) != len(data[0]):
            raise ValueError("The length of the data columns are not equal")

    df = pd.DataFrame(data_dict)

    write_sequence = range(len(data))  # key values to write data in order

    # Write dataframe to file
    df.to_csv(filename, columns=write_sequence, sep=sep, header=header, index=index, float_format=float_format)

    return 0

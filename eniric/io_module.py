"""Functions to read column-separated files.

These are a wrapper around pandas.
"""
from typing import List, Optional, Tuple

import pandas as pd
from numpy import ndarray


# noinspection SpellCheckingInspection,SpellCheckingInspection
def pdread_2col(filename: str, noheader: bool = False) -> Tuple[ndarray, ndarray]:
    """Read in a 2 column file with pandas.

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
        First column as float.
    col2: ndarray
        Second column as float.
    """
    try:
        if noheader:
            data = pd.read_table(
                filename,
                comment="#",
                names=["col1", "col2"],
                header=None,
                dtype=float,
                delim_whitespace=True,
            )
        else:
            data = pd.read_table(
                filename,
                comment="#",
                names=["col1", "col2"],
                dtype=float,
                delim_whitespace=True,
            )
    except Exception as e:
        print("There was an error trying to read in the file \n{}".format(filename))
        raise e

    return data["col1"].values, data["col2"].values


def pdread_3col(
    filename: str, noheader: bool = False
) -> Tuple[ndarray, ndarray, ndarray]:
    """Read in a 3 column file with pandas.

    Parameters
    ----------
    filename: str
        Name of file to read.
    noheader: bool
        Flag indicating if there is no column names given in file

    Returns
    -------
    col1: ndarray
        First column as float.
    col2: ndarray
        Second column as float.
    col3: ndarray
        Third column as float.
    """
    try:
        if noheader:
            data = pd.read_table(
                filename,
                comment="#",
                names=["col1", "col2", "col3"],
                header=None,
                dtype=float,
                delim_whitespace=True,
            )
        else:
            data = pd.read_table(
                filename,
                comment="#",
                names=["col1", "col2", "col3"],
                dtype=float,
                delim_whitespace=True,
            )
    except Exception as e:
        print("There was an error trying to read in the file \n{}".format(filename))
        raise e

    return data["col1"].values, data["col2"].values, data["col3"].values


def pdread_4col(
    filename: str, noheader: bool = False
) -> Tuple[ndarray, ndarray, ndarray, ndarray]:
    """Read in a 4 column file with pandas.

    Parameters
    ----------
    filename: str
        Name of file to read.
    noheader: bool
        Flag indicating if there is no column names given in file

    Returns
    -------
    col1: ndarray
        First column as float.
    col2: ndarray
        Second column as float.
    col3: ndarray
        Third column as float.
    col4: ndarray
        Fourth column as float.
    """
    try:
        if noheader:
            data = pd.read_table(
                filename,
                comment="#",
                names=["col1", "col2", "col3", "col4"],
                header=None,
                dtype=float,
                delim_whitespace=True,
            )
        else:
            data = pd.read_table(
                filename,
                comment="#",
                names=["col1", "col2", "col3", "col4"],
                dtype=float,
                delim_whitespace=True,
            )
    except Exception as e:
        print("There was an error trying to read in the file \n{}".format(filename))
        raise e

    return (
        data["col1"].values,
        data["col2"].values,
        data["col3"].values,
        data["col4"].values,
    )


def pdwrite_2col(
    filename: str,
    data1: ndarray,
    data2: ndarray,
    sep: str = "\t",
    header: Optional[List[str]] = None,
    float_format: Optional[str] = None,
) -> int:
    """Write out a 2 column file with pandas.

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
    header: Optional list of strings
        Header strings to apply to columns.
    float_format: str default None
        Specify floating point string format.

    Returns
    -------
    flag: bool
        Returns 0 if successful.
    """
    if header is not None:
        df = pd.DataFrame({"# {}".format(header[0]): data1, header[1]: data2})
    else:
        df = pd.DataFrame({"# x": data1, "y": data2})

    # Write DataFrame to file
    df.to_csv(
        filename, sep=sep, header=header, index=False, float_format=float_format
    )  # header=False

    return 0


def pdwrite_3col(
    filename: str,
    data1: ndarray,
    data2: ndarray,
    data3: ndarray,
    sep: str = "\t",
    header: Optional[List[str]] = None,
    float_format: Optional[str] = None,
) -> int:
    """Write out a 3 column file with pandas.

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
    header: optional list of strings
        Header strings to apply to columns.
    float_format: str default None
        Specify floating point string format.

    Returns
    -------
    flag: bool
        Returns 0 if successful.
    """
    if header is not None:
        df = pd.DataFrame(
            {"# {}".format(header[0]): data1, header[1]: data2, header[2]: data3}
        )
    else:
        df = pd.DataFrame({"# x": data1, "y": data2, "z": data3})

    # Write DataFrame to file
    df.to_csv(
        filename, sep=sep, header=header, index=False, float_format=float_format
    )  # header=False

    return 0


def pdwrite_cols(filename: str, *data, **kwargs) -> int:
    """Write out a csv file with pandas, variable columns possible.

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
    header: optional list of strings or bool
        Header strings to apply to columns. Must be equal to number
        of data columns provided.

    Returns
    -------
    flag: bool
        Returns 0 if successful.
    """

    # unpack keyword args, second argument is the default if not found.
    header = kwargs.pop("header", None)
    sep = kwargs.pop("sep", "\t")
    index = kwargs.pop("index", False)
    float_format = kwargs.pop("float_format", "%.6f")
    # TODO: See about passing any extra keywords into pandas call
    if kwargs:  # check for unwanted key words
        raise TypeError("Unexpected **kwargs: {!r}".format(kwargs))

    if header is not None:
        if len(header) != len(data):
            raise ValueError("Size of data and header does not match.")

    data_dict = {}
    for i, data_i in enumerate(data):
        data_dict[i] = data_i  # keys are assigned the index value from enumerate

        if len(data[i]) != len(data[0]):
            raise ValueError("The length of the data columns are not equal")

    df = pd.DataFrame(data_dict)

    write_sequence = range(len(data))  # key values to write data in order

    # Write DataFrame to file
    df.to_csv(
        filename,
        columns=write_sequence,
        sep=sep,
        header=header,
        index=index,
        float_format=float_format,
    )

    return 0

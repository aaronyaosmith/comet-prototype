import pandas as pd


def read_marker_csv(path):
    """Parses marker expression CSV data into a DataFrame.

    Path must contain a CSV formatted file with ',' delimiters, '.' decimal
    markings. Rows correspond to cell name, columns to gene name. Value is gene
    expression expressed as a floating point number.

    :param path: Location of marker expression CSV file.
    :returns: A DataFrame containing marker expression data; cells are rows,
    genes columns.
    :rtype: DataFrame

    """

    return pd.DataFrame()

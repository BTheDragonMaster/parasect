# -*- coding: utf-8 -*-

"""Module for reading tabular data files."""

import os
from collections import OrderedDict
from typing import List, Union, Any


def write_tabular(dictionaries: list[dict[str, Any]], header: list[str], out_file: str) -> None:
    """
    Write tabular file from list of dictionaries, sorted by the first column

    :param dictionaries: list of dictionaries. Keys of all dictionaries should be the same
    :type dictionaries: list[dict[str, Any]]
    :param header: list of categories for labelling the header. Length must be one more than the list of dictionaries.
    First element represents category describing dictionary keys
    :type header: list[str]
    :param out_file: path to output file
    :type out_file: str
    """

    assert len(dictionaries) + 1 == len(header)
    assert len(header) >= 2
    assert len(dictionaries) >= 1

    with open(out_file, 'w') as out:

        header_string = '\t'.join(header)
        out.write(f"{header_string}\n")

        keys = list(dictionaries[0].keys())
        keys.sort()

        for key in keys:
            row = [key]
            for dictionary in dictionaries:
                row.append(dictionary[key])

            row_string = '\t'.join(row)
            out.write(f"{row_string}\n")


class Tabular:
    """Class for reading and storing tabular data files."""

    def __init__(self, path_in: str, separator: str = '\t') -> None:
        """Initialize the Tabular class.

        :param path_in: Path to tabular file.
        :type path_in: str
        :param separator: Separator used in the tabular file. Default: tab
        :type separator: str
        :raises FileNotFoundError: If the file at the specified path does not exist.
        :raises ValueError: If the number of columns in a row does not match the
            length of the number of columns in the header.
        :raises ValueError: If there is a duplicate row ID in the data.
        """
        self.index_of_column_with_id = 0  # default index of column with ID
        self.column_names = []
        self.rows: OrderedDict = OrderedDict()

        # check if the file exists
        if not os.path.exists(path_in):
            raise FileNotFoundError(f"file not found: {path_in}")

        # read the tabular file
        with open(path_in, "r") as fo:
            self.column_names = [n.strip() for n in fo.readline().strip().split(separator)]

            for line_idx, line in enumerate(fo):
                # split the line into a list of values
                row = [v.strip() for v in line.strip().split(separator)]

                # check if the row has the same number of columns as the header
                if len(row) != len(self.column_names):
                    msg = f"row {line_idx + 2} has a different number of columns than the header"  # noqa: E501
                    raise ValueError(msg)

                # parse out the row ID
                row_id = row[self.index_of_column_with_id]

                # check for duplicate row IDs
                if row_id in self.rows:
                    msg = f"duplicate row ID when reading {path_in}: {row_id}"
                    raise ValueError(msg)

                # if the row ID is unique, add the row to the data dictionary
                self.rows[row_id] = OrderedDict()

                for value_idx, value in enumerate(row):
                    column_name = self.column_names[value_idx]
                    self.rows[row_id][column_name] = value

    def get_column_values(self, column_name: str) -> List[Union[int, float, str]]:
        """Return a list of values from a specified column.

        :param column_name: Name of the column.
        :type column_name: str
        :return: List of values from the specified column.
        :rtype: List[Union[int, float, str]]
        :raises KeyError: If the specified column name is not found in the data.
        """
        # check if the column name is in the column names
        if column_name not in self.column_names:
            raise KeyError(f"cannot find category {column_name} in data")

        # get the values from the specified column
        column_values = []
        for row_id in self.rows:
            column_values.append(self.get_row_value(row_id, column_name))

        return column_values

    def get_row_values(self, row_id: str) -> List[Union[int, float, str]]:
        """Return a row of values from the data.

        :param row_id: ID of the row.
        :type row_id: str
        :return: Row of values from the data.
        :rtype: List[Union[int, float, str]]
        :raises KeyError: If the specified row ID is not found in the data.
        """
        # check if the row ID is in the data
        if row_id not in self.rows:
            raise KeyError(f"cannot find data ID {row_id} in data")

        # get the values from the specified row
        row_values = []
        for category in self.column_names:
            row_values.append(self.get_row_value(row_id, category))

        return row_values

    def get_row_value(self, row_id: str, column_name: str) -> Union[int, float, str]:
        """Return a value from a specified row and column.

        :param row_id: ID of the row.
        :type row_id: str
        :param column_name: Name of the column.
        :type column_name: str
        :return: Value from the specified row and column.
        :rtype: Union[int, float, str]
        """
        return self.rows[row_id][column_name]

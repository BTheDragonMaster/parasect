# -*- coding: utf-8 -*-

"""Parsers module for PARASECT."""

import itertools
import logging
import os
from collections import OrderedDict
from typing import Dict, List, Tuple, Union

from Bio import SeqIO

from parasect.core.chem import smiles_to_fingerprint


class Tabular:
    """Class for reading and storing tabular data files."""

    def __init__(self, path_in: str, separator: str = "\t") -> None:
        """Initialize the Tabular class.

        :param path_in: Path to tabular file.
        :type path_in: str
        :param separator: Separator used in the tabular file (default is tab).
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


def parse_substrate_list(path_in: str) -> List[str]:
    """Parse a list of substrate names from a file.

    :param path_in: Path to input file.
    :type path_in: str
    :return: List of substrates.
    :rtype: List[str]
    :raises FileNotFoundError: If the file at the specified path does not exist.
    """
    # check if the file exists
    if not os.path.exists(path_in):
        raise FileNotFoundError(f"substrate list file not found: {path_in}")

    # initialize the list of substrates
    substrate_names = []

    # read the file
    with open(path_in, "r") as fo:
        for line in fo:
            substrate_name = line.strip()

            # check if the line is not empty
            if substrate_name:
                substrate_names.append(substrate_name)

    return substrate_names


def parse_morgan_fingerprint_file(path_in: str) -> Dict[str, List[int]]:
    """Parse a dictionary of molecule name to Morgan fingerprint from a file.

    :param path_in: Path to input file.
    :type path_in: str
    :return: Dictionary of molecule name to Morgan fingerprint.
    :rtype: Dict[str, List[int]]
    :raises FileNotFoundError: If the file at the specified path does not exist.
    """
    # check if the file exists
    if not os.path.exists(path_in):
        raise FileNotFoundError(f"morgan fingerprint file not found: {path_in}")

    # read the tabular file
    data = Tabular(path_in)

    # parse the Morgan fingerprints from the file
    amino_acid_fingerprints = OrderedDict()
    for row_id in data.rows:
        # skip the molecule name and SMILESMILES
        fingerprint = [int(v) for v in data.get_row_values(row_id)[2:]]

        # the row_id is the molecule name
        amino_acid_fingerprints[row_id] = fingerprint

    # print(amino_acid_fingerprints)
    return amino_acid_fingerprints


def parse_fasta_file(path_in: str) -> Dict[str, str]:
    """Read fasta file and parses it into a dictionary format.

    :param path_in: Path to fasta file.
    :type path_in: str
    :return: Dictionary of fasta sequences, where the key is the sequence header
        and the value is the sequence.
    :rtype: Dict[str, str]
    :raises FileNotFoundError: If the file at the specified path does not exist.
    """
    # check if the file exists
    if not os.path.exists(path_in):
        raise FileNotFoundError(f"File not found: {path_in}")

    # initialize the dictionary
    fasta_dict: Dict[str, str] = {}

    with open(path_in, "r") as fo:

        # initialize the header and sequence list
        header = ""
        sequence: List[str] = []

        for line in fo:
            line = line.strip()

            if line.startswith(">"):
                if sequence:
                    # if the sequence list is not empty, join the list into a string
                    # and add it to the dictionary with the header as the key
                    fasta_dict[header] = "".join(sequence)

                    header = line[1:]  # remove the ">" from the header
                    sequence = []  # reset the sequence list
                else:
                    header = line[1:]

            else:
                sequence.append(line)

        # if the header is not None, this means there is a sequence that has not
        # yet been added to the dictionary
        if len(header) > 0:
            fasta_dict[header] = "".join(sequence)

    return fasta_dict


def write_fasta_file(fasta_dict: Dict[str, str], path_out: str) -> None:
    """Write a dictionary of fasta sequences to a file.

    :param fasta_dict: Dictionary of fasta sequences, where the key is the sequence
        header and the value is the sequence.
    :type fasta_dict: Dict[str, str]
    :param path_out: Path to output fasta file.
    :type path_out: str
    """
    with open(path_out, "w") as fo:

        # iterate over the dictionary items
        for header, sequence in fasta_dict.items():
            fo.write(f">{header}\n{sequence}\n")


def parse_genbank_file(path_in: str, path_out: str) -> None:
    """Parse protein sequences from a GenBank file and writes them to a fasta file.

    :param path_in: Path to input GenBank file.
    :type path_in: str
    :param path_out: Path to output fasta file.
    :type path_out: str
    :raises FileNotFoundError: If the file at the specified path does not exist.
    """
    # check if the file exists
    if not os.path.exists(path_in):
        raise FileNotFoundError(f"File not found: {path_in}")

    # initialize a counter for generating gene IDs
    counter = itertools.count()

    # initialize a dictionary to store the fasta sequences
    fasta_dict = {}

    # parse the GenBank file
    for record in SeqIO.parse(path_in, "genbank"):
        for feature in record.features:

            # check if the feature is a coding sequence (CDS)
            if feature.type == "CDS":

                # check if the feature has a translation
                if "translation" in feature.qualifiers:
                    sequence = feature.qualifiers["translation"]

                    # check if the feature has a protein ID, gene ID, or locus tag
                    # to use as the sequence ID. If not, generate a gene ID
                    if "protein_id" in feature.qualifiers:
                        seq_id = feature.qualifiers["protein_id"][0]
                    elif "gene_id" in feature.qualifiers:
                        seq_id = feature.qualifiers["gene_id"][0]
                    elif "locus_tag" in feature.qualifiers:
                        seq_id = feature.qualifiers["locus_tag"][0]
                    else:
                        seq_id = f"gene_{next(counter)}"

                    # add the sequence to the dictionary
                    fasta_dict[seq_id] = sequence

    # write the fasta sequences to a file
    write_fasta_file(fasta_dict=fasta_dict, path_out=path_out)


def bitvectors_from_substrate_names(
    substrate_names: List[str],
    path_in_fingerprint_file: str,
) -> Tuple[List[str], List[List[int]]]:
    """Return substrate names and fingerprints for all substrate names for which a fingerprint could be found.

    :param substrate_names: Substrate names. If the
        substrate name is not in the fingerprint file, a warning will be logged.
    :type substrate_names: List[str]
    :param path_in_fingerprint_file: Path to file containing precomputed fingerprints,
        with one substrate per row and one substructure per column.
    :type path_in_fingerprint_file: str
    :return: Substrate names and fingerprints.
    :rtype: Tuple[List[str], List[List[int]]]

    .. note:: The index of each fingerprint matches the index of the corresponding
        substrate name in substrates.
    """
    logger = logging.getLogger(__name__)

    substrate_to_fingerprint = parse_morgan_fingerprint_file(path_in_fingerprint_file)

    ordered_substrate_names = []
    ordered_fingerprints = []

    for substrate_name in substrate_names:

        if substrate_name in substrate_to_fingerprint:
            ordered_substrate_names.append(substrate_name)
            ordered_fingerprints.append(substrate_to_fingerprint[substrate_name])

        else:
            msg = (
                f"Could not find a fingerprint for substrate name "
                f"{substrate_name}. Excluded from analysis."
            )
            logger.info(msg)

    return ordered_substrate_names, ordered_fingerprints


def bitvector_from_smiles(smiles: str, path_in_bitvector_file: str) -> List[int]:
    """Return a bitvector from a SMILES string.

    :param smiles: SMILES string.
    :type smiles: str
    :param path_in_bitvector_file: Path to file containing substrate bitvectors on which
        paras and parasect have been trained.
    :type path_in_bitvector_file: str
    :return: Bitvector, denoting the presence or absence of substructures in the SMILES string.
    :rtype: List[int]
    """
    fingerprint = smiles_to_fingerprint(smiles)

    with open(path_in_bitvector_file, "r") as bitvectors:
        substructure_hashes = list(map(int, bitvectors.readline().split("\t")[2:]))

    vector = []

    for substructure_hash in substructure_hashes:
        if substructure_hash in fingerprint:
            vector.append(1)
        else:
            vector.append(0)

    return vector

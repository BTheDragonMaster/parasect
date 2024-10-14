# -*- coding: utf-8 -*-

"""Constants used throughout the PARASECT package."""

import os
from typing import Dict, List, Union

import parasect.data
from parasect.core.tabular import Tabular


def _parse_amino_acid_properties_file(path_in: str) -> Dict[Union[int, float, str], List[float]]:
    """Parse a dictionary of amino acid to property vector from a file.

    :param path_in: Path to input file.
    :type path_in: str
    :return: Dictionary of amino acid to property vector.
    :rtype: Dict[Union[int, float, str], List[float]]
    :raises FileNotFoundError: If the file at the specified path does not exist.
    :raises ValueError: If the property vector for an amino acid is not of length 15.
    """
    # check if the file exists
    if not os.path.exists(path_in):
        raise FileNotFoundError(f"amino acid properties file not found: {path_in}")

    # read the tabular file
    data = Tabular(path_in=path_in, separator="\t")

    # parse the amino acid properties from the file
    amino_acid_properties = {}
    for row_id in data.rows:

        # get the amino acid ID and property vector
        amino_acid_id = data.get_row_value(row_id, column_name="AA")
        amino_acid_props = [float(v) for v in data.get_row_values(row_id)[1:]]

        # check if the property vector is of length 15
        if len(amino_acid_props) != 15:
            raise ValueError(
                f"amino acid property vector for {amino_acid_id} is not of length 15"
            )  # noqa: E501

        # add the amino acid ID and property vector to the dictionary
        amino_acid_properties[amino_acid_id] = amino_acid_props

    return amino_acid_properties


def _read_positions(path_in: str, start_position: int) -> List[int]:
    """Parse positions from a tab-separated file.

    :param path_in: Path to tab-separated input positions file.
    :type path_in: str
    :param start_position: Relative start position to adjust all positions by.
    :type start_position: int
    :return: List of positions, each relative to the start position.
    :rtype: List[int]
    """
    with open(path_in, "r") as fo:
        text = fo.read().strip()

        # split the text by tabs and adjust the positions by the start position
        positions = []
        for position in text.split("\t"):
            positions.append(int(position) - start_position)

    return positions


DATA_DIR = os.path.dirname(parasect.data.__file__)


def get_path(file_name: str) -> str:
    """Return the path to a file in the data directory.

    :param file_name: Name of the file.
    :type file_name: str
    :return: Path to the file.
    :rtype: str
    """
    return os.path.join(DATA_DIR, file_name)


ALIGNMENT_FILE = get_path("structure_alignment.fasta")
A_POSITION_FILE = get_path("stachelhaus.txt")
A_POSITION_FILE_34 = get_path("active_site.txt")
A_POSITION_FILE_HMM2 = get_path("stachelhaus_hmm2.txt")
A_POSITION_FILE_34_HMM2 = get_path("active_site_hmm2.txt")
SMILES_FILE = get_path("smiles.tsv")

FINGERPRINTS_FILE = os.path.join(DATA_DIR, "fingerprints.txt")
INCLUDED_SUBSTRATES_FILE = os.path.join(DATA_DIR, "included_substrates.txt")

HMM2_POSITIONS_SIGNATURE = _read_positions(A_POSITION_FILE_HMM2, 0)
HMM2_POSITIONS_EXTENDED_SIGNATURE = _read_positions(A_POSITION_FILE_34_HMM2, 0)
POSITIONS_SIGNATURE = _read_positions(A_POSITION_FILE, 66)
POSITIONS_EXTENDED_SIGNATURE = _read_positions(A_POSITION_FILE_34, 66)

HMM2_FILE = get_path("AMP-binding_hmmer2.hmm")
PROPERTIES = _parse_amino_acid_properties_file(get_path("physicochemical_properties.txt"))

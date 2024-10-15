# -*- coding: utf-8 -*-

"""Parsers module for PARASECT."""

import itertools
import logging
import os
from typing import Dict, List, Tuple

from Bio import SeqIO

from parasect.core.chem import smiles_to_fingerprint
from parasect.core.constants import FINGERPRINTS_FILE
from parasect.core.tabular import Tabular


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


def data_from_substrate_names(
    substrate_names: List[str],
) -> Tuple[List[str], List[str], List[List[int]]]:
    """Return substrate names and fingerprints for all substrate names for which a fingerprint could be found.

    :param substrate_names: Substrate names. If the
        substrate name is not in the fingerprint file, a warning will be logged.
    :type substrate_names: List[str]
    :return: Substrate names, SMILES strings, and fingerprints.
    :rtype: Tuple[List[str], List[str], List[List[int]]]

    .. note:: The index of each fingerprint matches the index of the corresponding
        substrate name in substrates.
    """
    logger = logging.getLogger(__name__)

    data = Tabular(path_in=FINGERPRINTS_FILE, separator="\t")

    ordered_substrate_names = []
    ordered_substrate_smiles = []
    ordered_fingerprints = []

    for name in substrate_names:

        if name in data.rows:
            # get the SMILES string and fingerprint for the substrate
            row_values = data.get_row_values(name)
            smiles = str(row_values[1])
            fingerprint = [int(v) for v in row_values[2:]]

            # add the substrate name, SMILES string, and fingerprint to the lists
            ordered_substrate_names.append(name)
            ordered_substrate_smiles.append(smiles)
            ordered_fingerprints.append(fingerprint)

        else:
            msg = f"could not find structure information for substrate name {name}"
            logger.info(msg)

    return (
        ordered_substrate_names,
        ordered_substrate_smiles,
        ordered_fingerprints,
    )


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

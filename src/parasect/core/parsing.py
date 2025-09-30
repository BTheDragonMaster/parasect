# -*- coding: utf-8 -*-

"""Parsers module for PARASECT."""

import logging
import os

from typing import Dict, List, Tuple, Optional, Generator
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from parasect.core.chem import smiles_to_fingerprint
from parasect.core.constants import FINGERPRINTS_FILE
from parasect.core.tabular import Tabular


@dataclass
class SubstrateData:
    name: str
    smiles: str


@dataclass
class TaxonomyData:
    domain: str
    kingdom: str
    phylum: str
    cls: str
    order: str
    family: str
    genus: str
    species: str
    strain: Optional[str]

    def __hash__(self):

        return hash((self.domain, self.kingdom, self.phylum, self.cls, self.order, self.family, self.genus, self.species, self.strain))

def parse_taxonomy_file(taxonomy_file: str) -> dict[str, TaxonomyData]:
    """Return dictionary of protein id to taxonomy

    :param taxonomy_file: taxonomy file
    :type taxonomy_file: str
    :return: dict of protein to taxonomy data object
    :rtype: dict[str, TaxonomyData
    """

    protein_to_taxonomy: dict[str, TaxonomyData] = {}
    tax_info = Tabular(taxonomy_file, separator='\t')
    for protein in tax_info.rows:
        domain = tax_info.get_row_value(protein, "domain")
        kingdom = tax_info.get_row_value(protein, "kingdom")
        phylum = tax_info.get_row_value(protein, "phylum")
        cls = tax_info.get_row_value(protein, "class")
        order = tax_info.get_row_value(protein, "order")
        family = tax_info.get_row_value(protein, "family")
        genus = tax_info.get_row_value(protein, "genus")
        species = tax_info.get_row_value(protein, "species")
        strain = tax_info.get_row_value(protein, "strain")
        if strain == "Unknown":
            strain = None
        protein_to_taxonomy[protein] = TaxonomyData(domain, kingdom, phylum, cls, order, family, genus, species, strain)

    return protein_to_taxonomy


def parse_raw_taxonomy(taxonomy_file: str) -> dict[str, list[str]]:
    """Return dictionary of protein id to taxonomy

    :param taxonomy_file: taxonomy file
    :type taxonomy_file: str
    :return: dictionary of protein ID to taxonomy
    :rtype: dict[str, list[str]]
    """
    protein_to_taxonomy: dict[str, list[str]] = {}
    with open(taxonomy_file, 'r') as tax_info:
        for line in tax_info:
            tax_data = line.strip().split('\t')
            protein = tax_data[0]
            taxonomy = tax_data[1:]
            protein_to_taxonomy[protein] = taxonomy
    return protein_to_taxonomy


def parse_list(list_file: str, sort: bool = True, unique: bool = True) -> list[str]:
    """Return list of things from file containing one item per line

    :param list_file: path to file containing one string per line
    :type list_file: str
    :param sort: if True, sort output list
    :type sort: bool
    :param unique: if True, return only unique items
    :type unique: bool
    :return: list of things
    :rtype: list[str]
    """
    string_list: list[str] = []
    with open(list_file, 'r') as list_of_things:
        for line in list_of_things:
            line = line.strip()
            if line:
                string_list.append(line)

    if unique:
        string_list = list(set(string_list))

    if sort:
        string_list.sort()

    return string_list


def parse_pcs(pca_file: str) -> tuple[list[str], NDArray[np.float64]]:
    """Return list of domain names and numpy array of precomputed principal components for each domain

    :param pca_file: file containing precomputed principal components
    :type pca_file: str
    :return: list of domain names and array of precomputed principal components for each domain
    :rtype: tuple[list[str], NDArray[np.float64]]
    """

    with open(pca_file, 'r') as pca_data:
        header = pca_data.readline()
        pc_nr = int(header.split('\t')[-1].split('_')[-1])
        domain_nr = 0
        for line in pca_data:
            line = line.strip()
            if line:
                domain_nr += 1

    array = np.zeros(shape=(domain_nr, pc_nr))
    domains = []

    with open(pca_file, 'r') as pca_data:
        pca_data.readline()
        counter = 0
        for line in pca_data:
            line = line.strip()
            if line:
                line_data = line.split('\t')
                pcs = list(map(float, line_data[1:]))
                array[counter] = pcs
                domains.append(line_data[0])
                counter += 1

    return domains, array


def parse_esm_embedding(embedding_path: str) -> NDArray[np.float64]:
    """Return ESM embedding from file

    :param embedding_path: path to file containing precomputed ESM embeddings
    :type embedding_path: str
    :return: Numpy array of full, ~34,000 feature embedding
    :rtype: numpy.ndarray
    """
    embedding: list[float] = []
    with open(embedding_path, 'r') as embedding_file:

        embedding_file.readline()
        for line in embedding_file:
            line = line.strip()
            if line:
                features = line.split('\t')[1:]
                embedding.extend(list(map(float, features)))

    return np.array(embedding, dtype=np.float64)


def parse_parasect_data(parasect_path: str) -> tuple[dict[str, str], dict[str, list[str]]]:
    """Return two dictionaries, with ids as keys and sequences/specificities as values from parasect data file

    :param parasect_path: path to parasect data file
    :type parasect_path: str
    :return: Two dictionaries, one mapping id to sequence and the other mapping id to specificity
    :rtype: tuple[dict[str, str], dict[str, list[str]]]
    """
    id_to_seq: dict[str, str] = {}
    id_to_spec: dict[str, list[str]] = {}

    parasect_data = Tabular(parasect_path, separator='\t')
    for domain_id in parasect_data.rows:
        id_to_seq[domain_id] = parasect_data.get_row_value(domain_id, "sequence")
        id_to_spec[domain_id] = parasect_data.get_row_value(domain_id, "specificity").split('|')

    return id_to_seq, id_to_spec


def parse_smiles_mapping(path_in: str) -> List[SubstrateData]:
    """Return substrates from SMILES mapping

    :param path_in: input path containing SMILES mapping
    :type path_in: str
    :return: list of substrates
    :rtype: list[SubstrateData]
    :raises FileNotFoundError: If the file at the specified path does not exist
    """
    if not os.path.exists(path_in):
        raise FileNotFoundError(f"SMILES mapping file not found: {path_in}")

    substrates = []
    with open(path_in, 'r') as input_file:
        input_file.readline()
        for line in input_file:
            line = line.strip()
            if line:
                name, smiles = line.split('\t')
                substrate = SubstrateData(name, smiles)
                substrates.append(substrate)

    return substrates


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


def iterate_over_dir(directory: str, extension: Optional[str] = None, get_dirs: bool = False) -> Generator[tuple[str, str], None, None]:
    """Yield the file/dir name and file/dir path of each file/dir in the specified path

    :param directory: path to directory to iterate over
    :type directory: str
    :param extension: extension of files to retrieve
    :type extension: Optional[str]
    :param get_dirs: if True, retrieve subdirectories of directory. If False, retrieve files in directory.
    :type get_dirs: bool
    :return: generator of tuples of file/dir name and file/dir path
    :rtype: Generator[str, str]
    """
    for file_name in os.listdir(directory):
        if not get_dirs:
            if file_name.endswith(extension):
                file_label = file_name.split(extension)[0]
                file_path = os.path.join(directory, file_name)
                yield file_label, file_path
        else:
            file_path = os.path.join(directory, file_name)
            if os.path.isdir(file_path):
                yield file_name, file_path

# -*- coding: utf-8 -*-

"""Module for interacting with pikachu library."""

import os
from typing import List, Dict
from math import isclose
from collections import Counter

from pikachu.fingerprinting.ecfp_4 import ECFP
from pikachu.general import read_smiles
from pikachu.fingerprinting.similarity import get_jaccard_index
from pikachu.errors import StructureError

from parasect.database.build_database import Substrate


def parse_smiles_mapping(path_in: str) -> Dict[str, str]:
    """Parse a name to SMILES mapping from a file.

    :param path_in: Path to input file
    :type path_in: str
    :return: Dict of {substrate_name: smiles_string, ->}
    :rtype: Dict[str, str]

    """
    if not os.path.exists(path_in):
        raise FileNotFoundError(f"SMILES mapping file not found: {path_in}")

    name_to_smiles = {}
    with open(path_in, 'r') as input_file:
        input_file.readline()
        for line in input_file:
            line = line.strip()
            name, smiles = line.split('\t')
            name_to_smiles[name] = smiles

    return name_to_smiles


def is_same_molecule(smiles_1, smiles_2):
    structure_1 = read_smiles(smiles_1)
    structure_2 = read_smiles(smiles_2)

    atom_number_1 = len(structure_1.get_atoms())
    atom_number_2 = len(structure_2.get_atoms())

    if atom_number_1 != atom_number_2:
        return False

    if get_jaccard_index(structure_1, structure_2) == 1.0:
        return True

    return False


def is_same_molecule_fingerprint(fingerprint_1: set[int], fingerprint_2: set[int]) -> bool:
    """Return True if the Jaccard index between the molecules is 1.0, False otherwise

    :param fingerprint_1: ECFP fingerprint of molecule 1
    :type fingerprint_1: set[int]
    :param fingerprint_2: ECFP fingerprint of molecule 2
    :type fingerprint_2: set[int]
    :return: True if molecules are the same, False otherwise
    :rtype: bool
    """

    jaccard_index = len(fingerprint_1.intersection(fingerprint_2)) / len(fingerprint_1.union(fingerprint_2))

    if isclose(jaccard_index, 1.0):
        return True

    return False


def smiles_to_fingerprint(smiles: str) -> set[int]:
    """Convert SMILES string to ECFP fingerprint.

    :param smiles: SMILES string.
    :type smiles: str
    :return: ECFP fingerprint.
    :rtype: List[int]
    """
    try:
        structure = read_smiles(smiles)
    except Exception:
        raise StructureError(f"Invalid SMILES string: {smiles}")

    ecfp = ECFP(structure)
    fingerprint = ecfp.fingerprint

    return fingerprint


def get_fingerprint_hashes(substrates: set[Substrate], fingerprint_size: int) -> list[int]:
    """Return list of fingerprint hashes sorted by count

    :param substrates: set of substrates from the PARASECT database
    :type substrates: set[Substrate]
    :param fingerprint_size: size of fingerprint
    :type fingerprint_size: int
    :return: list of hashes
    :rtype: list[int]
    """
    hashes = []

    for substrate in substrates:
        for chem_hash in substrate.fingerprint:
            hashes.append(chem_hash)

    counts = Counter(hashes)

    hashes = []

    for fingerprint_hash, count in counts.most_common(fingerprint_size):
        hashes.append(fingerprint_hash)

    return hashes


def fingerprint_to_bitvector(hashes: list[int], fingerprint: set[int]) -> list[int]:
    """Return bitvector from fingerprint and list of included hashes

    :param hashes: list of included hashes representing substrate training set
    :type hashes: list[int]
    :param fingerprint: fingerprint for substrate of interest
    :type fingerprint: set[int]
    :return: bitvector representing substrate of interest
    :rtype: list[int]
    """
    bitvector = []
    for fingerprint_hash in hashes:
        if fingerprint_hash in fingerprint:
            bitvector.append(1)
        else:
            bitvector.append(0)

    return bitvector


def smiles_to_bitvector(hashes: list[int], smiles: str) -> list[int]:
    """Return bitvector from SMILES string

    :param hashes: list[int]
    :type hashes: list
    :param smiles: SMILES string
    :type smiles: str
    :return: bitvector
    :rtype: list[int]
    """

    fingerprint = smiles_to_fingerprint(smiles)
    bitvector = fingerprint_to_bitvector(hashes, fingerprint)

    return bitvector

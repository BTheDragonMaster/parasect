# -*- coding: utf-8 -*-

"""Module for interacting with pikachu library."""

import os
from typing import List, Optional, Dict
from dataclasses import dataclass

from pikachu.fingerprinting.ecfp_4 import ECFP
from pikachu.general import read_smiles
from pikachu.fingerprinting.similarity import get_jaccard_index
from pikachu.errors import StructureError


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


@dataclass
class Substrate:
    name: str
    smiles: str


def is_same_molecule(smiles_1, smiles_2):
    structure_1 = read_smiles(smiles_1)
    structure_2 = read_smiles(smiles_2)

    if get_jaccard_index(structure_1, structure_2) == 1.0:
        atom_number_1 = len(structure_1.get_atoms())
        atom_number_2 = len(structure_2.get_atoms())
        if atom_number_1 == atom_number_2:
            return True

    return False


def smiles_to_fingerprint(smiles: str) -> List[int]:
    """Convert SMILES string to ECFP fingerprint.

    :param smiles: SMILES string.
    :type smiles: str
    :return: ECFP fingerprint.
    :rtype: List[int]
    """
    structure = read_smiles(smiles)

    ecfp = ECFP(structure)
    fingerprint = ecfp.fingerprint

    return fingerprint


def molecule_in_dataset(smiles_string: str, substrate_dataset: str) -> Optional[Substrate]:
    """Return True if molecule already exists within dataset, False otherwise

    :param smiles_string: SMILES string
    :type smiles_string: str
    :param substrate_dataset: path to file containing substrate dataset
    :return: Substrate if match was found, None otherwise
    :rtype: Optional[Substrate]
    """
    try:
        read_smiles(smiles_string)
    except Exception:
        raise StructureError("Invalid SMILES string")

    substrate_to_smiles = parse_smiles_mapping(substrate_dataset)
    for substrate, smiles in substrate_to_smiles.items():
        if is_same_molecule(smiles_string, smiles):
            return Substrate(substrate, smiles)

    return None


def substrate_name_in_dataset(substrate_name: str, substrate_dataset: str) -> Optional[Substrate]:
    """Return True if substrate name is already in dataset, False otherwise

    :param substrate_name: name of molecule
    :type substrate_name: str
    :param substrate_dataset: path to file containing substrate dataset
    :type substrate_dataset: str
    :return: Substrate if match was found, None otherwise
    :rtype: Optional[Substrate]
    """
    substrate_to_smiles = parse_smiles_mapping(substrate_dataset)
    for substrate, smiles in substrate_to_smiles.items():
        if substrate_name == substrate:
            return Substrate(substrate, smiles)

    return None

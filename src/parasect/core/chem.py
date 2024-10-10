# -*- coding: utf-8 -*-

"""Module for interacting with pikachu library."""

from typing import List

from pikachu.fingerprinting.ecfp_4 import ECFP
from pikachu.general import read_smiles


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

import os
from typing import TYPE_CHECKING
from sys import argv
import warnings

from pikachu.fingerprinting.ecfp_4 import build_ecfp_bitvector, ECFP
from pikachu.general import read_smiles

if TYPE_CHECKING:
    from pikachu.chem.structure import Structure

from paras.scripts.parsers.tabular import Tabular
from paras.scripts.parsers.parsers import parse_morgan_fingerprint


def bitvector_from_smiles(smiles_string: str, bitvector_file: str) -> list[int]:
    """
    Return a bitvector from a SMILES string

    Parameters
    ----------
    smiles_string: str, must be a valid SMILES string
    bitvector_file: path to file containing bitvectors on which paras and parasect have been trained

    Returns
    -------
    vector: list of int, with each int either 0 or 1, 0 denoting absence and 1 denoting prescence of bitvector. Vector
        is given in the same order as in bitvector_file

    """
    structure: Structure = read_smiles(smiles_string)

    ecfp: ECFP = ECFP(structure)
    fingerprint: set[int] = ecfp.fingerprint

    with open(bitvector_file, 'r') as bitvectors:
        substructure_hashes = list(map(int, bitvectors.readline().split('\t')[2:]))

    vector: list[int] = []

    for substructure_hash in substructure_hashes:
        if substructure_hash in fingerprint:
            vector.append(1)
        else:
            vector.append(0)

    return vector


def bitvectors_from_substrate_names(substrate_names, fingerprint_file):
    """
    Return substrate names and fingerprints for all substrate names for which a fingerprint could be found

    Parameters
    ----------
    substrate_names: list of [str, ->], with each str a substrate name. If the substrate name is not in the fingerprint
        file, a warning will be printed.
    fingerprint_file: str, path to file containing precomputed fingerprints, with one substrate per row and one
        substructure per column

    Returns
    -------
    substrates: list of [str, ->], substrate names which were present in the fingerprint file
    fingerprints: list of [[int, ->], ->], with each list of integers a fingerprint. The index of each fingerprint
        matches the index of the corresponding substrate name in substrates.

    """
    substrate_to_fingerprint: dict[str, list[int]] = parse_morgan_fingerprint(fingerprint_file)
    substrates: list[str] = []
    fingerprints: list[list[int]] = []
    for substrate in substrate_names:
        if substrate in substrate_to_fingerprint:
            substrates.append(substrate)
            fingerprints.append(substrate_to_fingerprint[substrate])
        else:
            warnings.warn(f"Could not find a fingerprint for substrate name {substrate}. Excluded from analysis.")
    return substrates, fingerprints


def bitvector_from_substrates(substrate_file, out_dir):
    """
    Compute bitvectors from SMILES strings and write them to a file

    Parameters
    ----------
    substrate_file: path to file containing one column labelled 'substrate' and another labelled 'smiles'
    out_dir: path to output file

    """
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    substrate_data: Tabular = Tabular(substrate_file, [0])
    smiles_strings: list[str] = substrate_data.get_column('smiles')
    names: list[str] = substrate_data.get_column('substrate')
    structures: list[Structure] = []

    for smiles_string in smiles_strings:
        structures.append(read_smiles(smiles_string))

    bitvector, mapping, fingerprints = build_ecfp_bitvector(structures)

    vectors: list[list[int]] = []

    for fingerprint in fingerprints:
        vector: list[int] = []
        for substructure_hash in bitvector:
            if substructure_hash in fingerprint:
                vector.append(1)
            else:
                vector.append(0)
        vectors.append(vector)

    image_dir: str = os.path.join(out_dir, 'feature_images')
    if not os.path.exists(image_dir):
        os.mkdir(image_dir)

    for bit, substructure in mapping.items():
        file_name: str = os.path.join(image_dir, f"{bit}.svg")
        substructure.draw(file_name)

    out_file: str = os.path.join(out_dir, 'fingerprints.txt')

    with open(out_file, 'w') as out:
        out.write("substrate_name\tsmiles\t")
        out.write("\t".join(map(str, bitvector)))
        out.write('\n')

        for i, name in enumerate(names):
            smiles = smiles_strings[i]
            vector = vectors[i]
            out.write(f"{name}\t{smiles}\t")
            out.write("\t".join(map(str, vector)))
            out.write('\n')


if __name__ == "__main__":
    bitvector_from_substrates(substrate_file=argv[1], out_dir=argv[2])

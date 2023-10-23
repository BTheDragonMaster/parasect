from pikachu.fingerprinting.ecfp_4 import build_ecfp_bitvector, ECFP
from pikachu.general import read_smiles
from paras.scripts.parsers.tabular import Tabular
from paras.scripts.parsers.parsers import parse_morgan_fingerprint
from sys import argv
import warnings


def bitvector_from_smiles(smiles_string, bitvector_file):
    structure = read_smiles(smiles_string)

    ecfp = ECFP(structure)
    fingerprint = ecfp.fingerprint

    with open(bitvector_file, 'r') as bitvectors:
        substructure_hashes = bitvectors.readline().split('\t')[2:]

    vector = []

    for substructure_hash in substructure_hashes:
        if substructure_hash in fingerprint:
            vector.append(1)
        else:
            vector.append(0)

    return vector


def bitvectors_from_substrate_names(substrate_names, fingerprint_file):
    substrate_to_fingerprint = parse_morgan_fingerprint(fingerprint_file)
    substrates = []
    fingerprints = []
    for substrate in substrate_names:
        if substrate in substrate_to_fingerprint:
            substrates.append(substrate)
            fingerprints.append(substrate_to_fingerprint[substrate])
        else:
            warnings.warn(f"Could not find a fingerprint for substrate name {substrate}. Excluded from analysis.")
    return substrates, fingerprints


def bitvector_from_substrates(substrate_file, out_file):
    substrate_data = Tabular(substrate_file, [0])
    smiles_strings = substrate_data.get_column('smiles')
    names = substrate_data.get_column('substrate')
    structures = []

    for smiles_string in smiles_strings:
        structures.append(read_smiles(smiles_string))

    bitvector, mapping, fingerprints = build_ecfp_bitvector(structures)

    vectors = []

    for fingerprint in fingerprints:
        vector = []
        for substructure_hash in bitvector:
            if substructure_hash in fingerprint:
                vector.append(1)
            else:
                vector.append(0)
        vectors.append(vector)

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
    bitvector_from_substrates(substrate_file=argv[1], out_file=argv[2])

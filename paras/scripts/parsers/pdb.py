from collections import OrderedDict
from sys import argv

from paras.scripts.math.shapes import Vector3D, Sphere
from paras.data.sequence_data.amino_acid_properties.amino_acid_dictionaries import AA_TABLE
from paras.scripts.parsers.atom import Atom, Residue


VAN_DER_WAALS_RADII = {"H": 1.2,
                       "C": 1.7,
                       "N": 1.55,
                       "O": 1.52,
                       "P": 1.8,
                       "S": 1.8,
                       "Mg": 1.73}


def atoms_from_pdb(pdb_file):
    atoms = []
    with open(pdb_file, 'r') as pdb:
        for line in pdb:
            if line.startswith("HETATM") or line.startswith("ATOM"):
                atom_name = line[12:16].strip()
                atom_type = line[76:78].strip().title()
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                location = Vector3D(x, y, z)
                residue_name = line[17:20].strip()
                residue_nr = int(line[22:26].strip())
                chain_id = line[21]
                if residue_name in AA_TABLE:
                    residue_type = AA_TABLE[residue_name]
                else:
                    residue_type = residue_name

                if atom_type in VAN_DER_WAALS_RADII:
                    atom_radius = VAN_DER_WAALS_RADII[atom_type]
                else:
                    raise ValueError(f"Unknown atom type: {atom_type}")

                sphere = Sphere(location, atom_radius)

                residue = Residue(residue_type, residue_nr)

                atom = Atom(atom_type, atom_name, residue, chain_id, sphere)
                atoms.append(atom)

    return atoms


def sequence_from_pdb(pdb_file):
    chain_to_residues = OrderedDict()

    with open(pdb_file, 'r') as pdb:
        for line in pdb:
            if line.startswith("HETATM") or line.startswith("ATOM"):

                residue_name = line[17:20].strip()
                residue_nr = int(line[22:26].strip())
                chain_id = line[21]
                if chain_id not in chain_to_residues:
                    chain_to_residues[chain_id] = []

                if residue_name in AA_TABLE:
                    residue_type = AA_TABLE[residue_name]
                    residue = Residue(residue_type, residue_nr)
                    if residue not in chain_to_residues[chain_id]:
                        chain_to_residues[chain_id].append(residue)

    sequence = []
    residue_nrs = []
    for chain_id, residues in chain_to_residues.items():
        if chain_id == 'A':
            for residue in residues:
                sequence.append(residue.type)
                residue_nrs.append(residue.nr)
    return ''.join(sequence), residue_nrs


if __name__ == "__main__":
    print(sequence_from_pdb(argv[1]))


import os
from sys import argv

from pikachu.general import svg_from_smiles

from paras.scripts.parsers.parsers import parse_substrate_list
from paras.scripts.parsers.parsers import parse_substrate_smiles


def draw_compounds(included_substrates, smiles_file, out_dir):
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    compound_to_smiles = parse_substrate_smiles(smiles_file)
    for compound in parse_substrate_list(included_substrates):
        svg_from_smiles(compound_to_smiles[compound], os.path.join(out_dir, f"{compound}.svg"))


if __name__ == "__main__":
    draw_compounds(argv[1], argv[2], argv[3])
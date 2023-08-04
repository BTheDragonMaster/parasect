from paras.scripts.parsers.parsers import parse_feature_importances
from paras.scripts.data_processing.structure_processing.add_dummy_substrate import add_dummy_substrate
import argparse
from pymol import cmd, stored
import os
from dataclasses import asdict

number_to_alphabet = {'0': "A",
                      '1': "B",
                      '2': "C",
                      '3': "D",
                      '4': "E",
                      '5': "F",
                      '6': "G",
                      '7': "H",
                      '8': "I",
                      '9': "J"}
alphabet_to_number = {"A": "0",
                      "B": "1",
                      "C": "2",
                      "D": "3",
                      "E": "4",
                      "F": "5",
                      "G": "6",
                      "H": "7",
                      "I": "8",
                      "J": "9"}


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type=str, required=True, help="Path to importances file")
    parser.add_argument('-o', type=str, required=True, help="Path to output directory")
    parser.add_argument('-w', type=float, default=1.0, help="Voxel width")
    parser.add_argument('-pdb', type=str, required=True, help="PDB file to superpose voxels on.")
    parser.add_argument('-d', type=str, default=None, help="Link to file with PDB dummy substrate")
    args = parser.parse_args()

    return args


def run():
    args = parse_arguments()
    if not os.path.exists(args.o):
        os.mkdir(args.o)

    voxels = parse_feature_importances(args.i, voxel_size=args.w)
    normalize_importances(voxels)
    if args.d:
        pdb_dummy = os.path.join(args.o, "domain_with_substrate.pdb")
        add_dummy_substrate(args.pdb, args.d, pdb_dummy)
        pdb_in = pdb_dummy
    else:
        pdb_in = args.pdb
    pdb_out = os.path.join(args.o, "pdbs")
    pymol_out = os.path.join(args.o, "visualisations")
    if not os.path.exists(pdb_out):
        os.mkdir(pdb_out)
    if not os.path.exists(pymol_out):
        os.mkdir(pymol_out)
    add_voxels_to_pdb(pdb_in, voxels, pdb_out)

    visualise_in_pymol(pdb_out, voxels, pymol_out)


def normalize_importances(voxels):
    max_importance = 0.0

    for voxel in voxels:
        for feature, importance in asdict(voxel.importances).items():
            if abs(importance) > max_importance:
                max_importance = abs(importance)

    scaling_factor = 1.0 / max_importance

    for voxel in voxels:
        for field in voxel.importances.__dataclass_fields__:
            value = getattr(voxel.importances, field)
            voxel.update_importances(field, value * scaling_factor)


def importance_to_alphabetic(importance):
    importance_int = int(abs(importance) * 100) - 1
    if importance_int < 0:
        importance_int = 0

    assert 0 <= importance_int < 100
    string = f'{importance_int:02d}'
    assert len(string) == 2
    alphabetic = ''
    for character in string:
        alphabetic += number_to_alphabet[character]

    if importance < 0:
        alphabetic += 'R'
    else:
        alphabetic += 'B'

    return alphabetic


def alphabetic_to_colour(alphabetic):
    importance = alphabetic_to_importance(alphabetic)

    if importance < 0:
        colour = (1.0, 1 - abs(importance), 1 - abs(importance))
    else:
        colour = (1 - abs(importance), 1 - abs(importance), 1.0)

    return colour


def alphabetic_to_importance(alphabetic):
    number = ''
    for character in alphabetic[:-1]:
        number += alphabet_to_number[character]

    number = float(number) / 100

    if alphabetic[-1] == 'R':
        number *= -1

    return number


def add_voxels_to_pdb(pdb_in, voxels, pdb_out):
    categories = voxels[0].importances.__dataclass_fields__
    for i, category in enumerate(categories):
        pdb_file = os.path.join(pdb_out, f"{category}.pdb")

        with open(pdb_in, 'r') as pdb:
            with open(pdb_file, 'w') as out:
                for line in pdb:
                    if line.strip() != "END":
                        out.write(line)

                for j, voxel in enumerate(voxels):
                    importance = getattr(voxel.importances, category)
                    alphabetic = importance_to_alphabetic(importance)

                    pdb_line = f"HETATM{j + 1:5d} {alphabetic}  STP F{1:4d}    {voxel.cube.midpoint.x:8.3f}{voxel.cube.midpoint.y:8.3f}{voxel.cube.midpoint.z:8.3f}  0.00  0.00          Ve\n"
                    out.write(pdb_line)
                if voxels:
                    out.write("TER\n")

                out.write("END\n")


def visualise_in_pymol(pdb_in, voxels, pse_out):
    categories = voxels[0].importances.__dataclass_fields__

    for i, category in enumerate(categories):
        out_file = os.path.join(pse_out, f"{category}.pse")
        print(f"Colouring category {category}")
        pdb_file = os.path.join(pdb_in, f"{category}.pdb")
        cmd.load(pdb_file)
        cmd.hide('lines', 'resn STP')
        cmd.hide('sticks', 'resn STP')

        cmd.select("voxels", f"resn STP")

        stored.letters = []
        cmd.iterate("resn STP", "stored.letters.append(name)")
        letter_pairs = list(set(stored.letters))
        colours = []
        for letter_pair in letter_pairs:
            colours.append(alphabetic_to_colour(letter_pair))

        for j, letter_pair in enumerate(letter_pairs):
            selection_name = f"{letter_pair}_{category}"
            cmd.select(selection_name, f"resn STP and name {letter_pair}")
            if not letter_pair.startswith('A') and not letter_pair.startswith('B'):
                colour = colours[j]
                cmd.set_color(letter_pair, list(colour))
                cmd.color(letter_pair, letter_pair)

                # cmd.set("sphere_transparency", f"{1 - abs(importance)}", letter_pair)
            else:
                cmd.remove(selection_name)
            cmd.show("spheres", selection_name)
            cmd.delete(selection_name)

            cmd.hide('nb_spheres', 'resn STP')
            cmd.set("sphere_scale", "0.3", 'resn STP')
        cmd.save(out_file)
        cmd.reinitialize()


if __name__ == "__main__":
    run()

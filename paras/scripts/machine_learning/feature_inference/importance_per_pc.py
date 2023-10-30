from paras.scripts.parsers.parsers import parse_feature_importances
from paras.scripts.data_processing.structure_processing.add_dummy_substrate import add_dummy_substrate, get_highest_serial_number
from paras.scripts.machine_learning.feature_inference.plot_importances_pymol import normalize_importances, \
    importance_to_alphabetic, alphabetic_to_colour, alphabetic_to_importance
import argparse
import os
from pymol import cmd, stored


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type=str, required=True, help="Path to importances directory")
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

    for file_name in os.listdir(args.i):
        file_path = os.path.join(args.i, file_name)
        if file_name.endswith('.txt') and os.path.isfile(file_path):
            pc_nr = int(file_name.split('_')[1]) + 1
            voxels = parse_feature_importances(file_path, voxel_size=args.w)
            normalize_importances(voxels)
            if args.d:
                pdb_dummy = os.path.join(args.o, "domain_with_substrate.pdb")
                add_dummy_substrate(args.pdb, args.d, pdb_dummy)
                pdb_in = pdb_dummy
            else:
                pdb_in = args.pdb
            pdb_out = os.path.join(args.o, f"pc_{pc_nr}.pdb")
            pymol_out = os.path.join(args.o, f"pc_{pc_nr}.pse")

            add_voxels_to_pdb(pdb_in, voxels, pdb_out)

            visualise_in_pymol(pdb_out, voxels, pymol_out)


def add_voxels_to_pdb(pdb_in, voxels, pdb_out):
    categories = voxels[0].importances.__dataclass_fields__
    highest_serial_number = get_highest_serial_number(pdb_in)

    with open(pdb_in, 'r') as pdb:
        with open(pdb_out, 'w') as out:
            for i, category in enumerate(categories):

                for line in pdb:
                    if line.strip() != "END":
                        out.write(line)

                for j, voxel in enumerate(voxels):
                    highest_serial_number += 1
                    importance = getattr(voxel.importances, category)
                    alphabetic = importance_to_alphabetic(importance)

                    pdb_line = f"HETATM{highest_serial_number:5d}  {alphabetic} STP F{i + 1:4d}    {voxel.cube.midpoint.x:8.3f}{voxel.cube.midpoint.y:8.3f}{voxel.cube.midpoint.z:8.3f}  0.00  0.00          Ve\n"
                    out.write(pdb_line)

                if voxels:
                    out.write("TER\n")

            out.write("END\n")


def visualise_in_pymol(pdb_in, voxels, pse_out):
    categories = voxels[0].importances.__dataclass_fields__

    pdb_file = os.path.join(pdb_in)
    cmd.load(pdb_file, object="1amu")
    cmd.hide('lines', 'resn STP')
    cmd.hide('sticks', 'resn STP')

    stored.letters = []
    cmd.iterate("resn STP", "stored.letters.append(name)")
    letter_pairs = list(set(stored.letters))

    for letter_pair in letter_pairs:
        colour = alphabetic_to_colour(letter_pair)
        cmd.set_color(letter_pair, list(colour))

    for i, category in enumerate(categories):
        cmd.select(category, f"resn STP and resi {i + 1}")
        print(f"Colouring category {category}")

        stored.pairs = []

        for j, letter_pair in enumerate(letter_pairs):
            selection_name = f"{letter_pair}_{category}"
            cmd.select(selection_name, f"resn STP and name {letter_pair} and resi {i + 1}")
            if not letter_pair.startswith('A') and not letter_pair.startswith('B') and not letter_pair.startswith('C')\
                    and not letter_pair.startswith('D') and not letter_pair.startswith('E'):
                cmd.color(letter_pair, letter_pair)
                cmd.iterate(letter_pair, "stored.pairs.append(name)")
                importance = alphabetic_to_importance(letter_pair)
                cmd.set("sphere_transparency", f"{1 - abs(importance)}", letter_pair)
            else:
                cmd.remove(selection_name)
            cmd.show("spheres", selection_name)
            cmd.delete(selection_name)

            cmd.hide('nb_spheres', 'resn STP')
            cmd.set("sphere_scale", "0.3", 'resn STP')

        if not stored.pairs:
            print(category)
            cmd.delete(category)
    cmd.save(pse_out)
    cmd.reinitialize()


if __name__ == "__main__":
    run()

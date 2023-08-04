from paras.scripts.feature_extraction.structure_feature_extraction.voxel_pockets import *
import argparse
import os
from paras.scripts.math.shapes import Vector3D


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type=str, required=True, help="Input directory containing pdb files of A-domains.")
    parser.add_argument('-o', type=str, required=True, help="Output directory.")

    args = parser.parse_args()
    return args


def run():

    args = parse_arguments()

    if not os.path.exists(args.o):
        os.mkdir(args.o)

    offset = Vector3D(0.0, 0.0, 0.0)

    counter = 0

    for pdb_file in os.listdir(args.i):
        pdb_path = os.path.join(args.i, pdb_file)

        if os.path.isfile(pdb_path) and pdb_path.endswith('.pdb'):
            counter += 1
            pdb_name = pdb_file.split('.pdb')[0]
            voxel_grid = pdb_to_voxel_representation(pdb_path, 1.2, 10.0, SUBSTRATE_LOCATION, offset)
            pocket = get_empty_voxels(voxel_grid)
            adjacent_voxels = label_adjacent_voxels(pocket, voxel_grid, 1.2)
            pdb_with_voxels = os.path.join(args.o, pdb_file)
            add_voxels_to_pdb(pdb_path, pocket, adjacent_voxels, pdb_with_voxels)
            visualise_in_pymol(pdb_with_voxels, pdb_name, args.o)

            if counter % 10 == 0:
                print(f"Processed {counter} pockets.")


if __name__ == "__main__":
    run()
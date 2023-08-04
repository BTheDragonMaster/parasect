from typing import List

import argparse
import os

from paras.scripts.math.shapes import Cube, Sphere, Vector3D
from paras.scripts.math.distance_calculations import intersection_cube_sphere
from paras.scripts.parsers.pdb import atoms_from_pdb
from paras.scripts.parsers.voxel import Voxel


BETA_CARBON_LOCATION = Vector3D(32.833,  97.884,  32.711)
ALPHA_CARBON_LOCATION = Vector3D(32.341,  99.049,  31.853)

SUBSTRATE_LOCATION = BETA_CARBON_LOCATION

FILTER_ORB = Sphere(Vector3D(40.243, 107.507, 26.930), 10.0)


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type=str, required=True, help="Input directory containing pdb files of A-domains.")
    parser.add_argument('-o', type=str, required=True, help="Output file.")
    parser.add_argument('-w', type=float, default=1.0, help="Voxel width")
    parser.add_argument('-r', type=float, default=10.0, help="Radius around substrate")
    parser.add_argument('-x', type=float, default=0.0, help="Offset of voxel sphere in x direction")
    parser.add_argument('-y', type=float, default=0.0, help="Offset of voxel sphere in y direction")
    parser.add_argument('-z', type=float, default=0.0, help="Offset of voxel sphere in z direction")
    args = parser.parse_args()
    return args


def run():

    args = parse_arguments()

    offset = Vector3D(args.x, args.y, args.z)

    counter = 0

    with open(args.o, 'w') as features:

        wrote_header = False

        for pdb_file in os.listdir(args.i):

            pdb_path = os.path.join(args.i, pdb_file)

            if os.path.isfile(pdb_path) and pdb_path.endswith('.pdb'):
                counter += 1
                pdb_name = pdb_file.split('.pdb')[0]
                voxel_grid = pdb_to_voxel_representation(pdb_path, args.w, args.r, SUBSTRATE_LOCATION, offset)
                if not wrote_header:
                    features.write('domain_name\t')
                    for i, voxel in enumerate(voxel_grid):
                        features.write(
                            f"{voxel.cube.midpoint.__repr__()}|nr_aromatic\t")
                        features.write(
                            f"{voxel.cube.midpoint.__repr__()}|nr_pos_charge\t")
                        features.write(
                            f"{voxel.cube.midpoint.__repr__()}|nr_neg_charge\t")
                        features.write(
                            f"{voxel.cube.midpoint.__repr__()}|nr_pos_polar\t")
                        features.write(
                            f"{voxel.cube.midpoint.__repr__()}|nr_neg_polar\t")
                        features.write(
                            f"{voxel.cube.midpoint.__repr__()}|nr_sulfur\t")
                        features.write(
                            f"{voxel.cube.midpoint.__repr__()}|nr_hydrophobic")

                        if i != len(voxel_grid) - 1:
                            features.write('\t')

                    features.write('\n')
                    wrote_header = True

                feature_vector = []
                for voxel in voxel_grid:
                    voxel.to_feature_vector()
                    voxel_vector = voxel.vector.as_list(["nr_aromatic",
                                                         "nr_pos_charge",
                                                         "nr_neg_charge",
                                                         "nr_pos_polar",
                                                         "nr_neg_polar",
                                                         "nr_sulfur",
                                                         "nr_hydrophobic"])
                    feature_vector += voxel_vector

                features.write(f'{pdb_name}\t')
                features.write('\t'.join(map(str, feature_vector)))
                features.write('\n')

                if counter % 10 == 0:
                    print(f"Processed {counter} structures.")


def get_empty_voxels(voxel_grid):
    empty_voxels = []
    for voxel in voxel_grid:
        if not voxel.atoms:
            empty_voxels.append(voxel)
    return empty_voxels


def make_voxel_grid(voxel_width: float, radius: float, center: Vector3D, offset: Vector3D) -> List[Voxel]:
    voxels = []

    nr_voxels = int(radius / voxel_width)
    min_x = center.x - nr_voxels * voxel_width + offset.x
    min_y = center.y - nr_voxels * voxel_width + offset.y
    min_z = center.z - nr_voxels * voxel_width + offset.z

    for i in range(nr_voxels * 2):
        for j in range(nr_voxels * 2):
            for k in range(nr_voxels * 2):
                midpoint = Vector3D(min_x + i * voxel_width + 0.5 * voxel_width,
                                    min_y + j * voxel_width + 0.5 * voxel_width,
                                    min_z + k * voxel_width + 0.5 * voxel_width)
                cube = Cube(midpoint, voxel_width / 2)
                voxel = Voxel(cube, [])
                voxels.append(voxel)

    return voxels


def pdb_to_voxel_representation(pdb_file: str, voxel_width: float,
                                distance_to_substrate: float,
                                substrate_location: Vector3D,
                                offset: Vector3D) -> List[Voxel]:
    voxel_grid = make_voxel_grid(voxel_width, distance_to_substrate, substrate_location, offset)
    atoms = atoms_from_pdb(pdb_file)

    candidate_atoms = []
    voxel_cube = Cube(substrate_location, distance_to_substrate)

    for atom in atoms:
        if intersection_cube_sphere(voxel_cube, atom.sphere):
            candidate_atoms.append(atom)

    for voxel in voxel_grid:
        for atom in candidate_atoms:
            if intersection_cube_sphere(voxel.cube, atom.sphere):
                voxel.atoms.append(atom)

    return voxel_grid


if __name__ == "__main__":
    run()

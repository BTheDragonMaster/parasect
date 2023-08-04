from typing import List

import argparse
import os

from pymol import cmd

from paras.scripts.math.shapes import Cube, Sphere, Vector3D
from paras.scripts.math.distance_calculations import intersection_cube_sphere
from paras.scripts.parsers.pdb import atoms_from_pdb
from paras.scripts.parsers.voxel import Voxel, VOXEL_TYPE_TO_ABBREVIATION
from paras.scripts.data_processing.structure_processing.voxel_adjacency_matrix import make_adjacency_matrix


BETA_CARBON_LOCATION = Vector3D(32.833,  97.884,  32.711)
ALPHA_CARBON_LOCATION = Vector3D(32.341,  99.049,  31.853)

SUBSTRATE_LOCATION = BETA_CARBON_LOCATION

FILTER_ORB = Sphere(Vector3D(40.243, 107.507, 26.930), 10.0)


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type=str, required=True, help="Input directory containing pdb files of A-domains.")
    parser.add_argument('-o', type=str, required=True, help="Output directory.")
    parser.add_argument('-w', type=float, default=1.0, help="Voxel width")
    parser.add_argument('-r', type=float, default=10.0, help="Radius around substrate")
    parser.add_argument('-x', type=float, default=0.0, help="Offset of voxel sphere in x direction")
    parser.add_argument('-y', type=float, default=0.0, help="Offset of voxel sphere in y direction")
    parser.add_argument('-z', type=float, default=0.0, help="Offset of voxel sphere in z direction")
    args = parser.parse_args()
    return args


def run():

    args = parse_arguments()

    if not os.path.exists(args.o):
        os.mkdir(args.o)

    offset = Vector3D(args.x, args.y, args.z)

    features_out = os.path.join(args.o, 'pocket_features.txt')
    features_simplified_out = os.path.join(args.o, 'pocket_features_simplified.txt')
    counter = 0

    with open(features_out, 'w') as features:
        with open(features_simplified_out, 'w') as features_simplified:
            wrote_header = False

            for pdb_file in os.listdir(args.i):
                counter += 1
                pdb_path = os.path.join(args.i, pdb_file)

                if os.path.isfile(pdb_path) and pdb_path.endswith('.pdb'):
                    pdb_name = pdb_file.split('.pdb')[0]
                    voxel_grid = pdb_to_voxel_representation(pdb_path, args.w, args.r, SUBSTRATE_LOCATION, offset)
                    if not wrote_header:
                        features.write('domain_name\t')
                        features_simplified.write('domain_name\t')
                        for i, voxel in enumerate(voxel_grid):
                            features.write(
                                f"{voxel.__repr__()}|is_empty\t")
                            features.write(
                                f"{voxel.cube.midpoint.__repr__()}|is_pocket\t")
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
                                f"{voxel.cube.midpoint.__repr__()}|nr_hydrophobic\t")
                            features.write(
                                f"{voxel.cube.midpoint.__repr__()}|nr_metal\t")
                            features.write(
                                f"{voxel.cube.midpoint.__repr__()}|nr_phosphorus")
                            features_simplified.write(f"{voxel.__repr__()}")
                            if i != len(voxel_grid) - 1:
                                features.write('\t')
                                features_simplified.write('\t')

                        features.write('\n')
                        features_simplified.write('\n')
                        wrote_header = True
                    pocket = get_empty_voxels(voxel_grid)
                    adjacent_voxels = label_adjacent_voxels(pocket, voxel_grid, args.w)
                    pdb_with_voxels = os.path.join(args.o, pdb_file)
                    add_voxels_to_pdb(pdb_path, pocket, adjacent_voxels, pdb_with_voxels)
                    visualise_in_pymol(pdb_with_voxels, pdb_name, args.o)
                    simple_vector = []
                    feature_vector = []
                    for voxel in voxel_grid:
                        feature_vector += voxel.to_feature_vector()
                        if voxel.atoms:

                            simple_vector.append(0)
                        else:
                            simple_vector.append(1)
                    features.write(f'{pdb_name}\t')
                    features_simplified.write(f'{pdb_name}\t')

                    features.write('\t'.join(map(str, feature_vector)))
                    features_simplified.write('\t'.join(map(str, simple_vector)))

                    features.write('\n')
                    features_simplified.write('\n')

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
    sphere = Sphere(center, radius)

    for i in range(nr_voxels * 2 + 1):
        for j in range(nr_voxels * 2 + 1):
            for k in range(nr_voxels * 2 + 1):
                midpoint = Vector3D(min_x + i * voxel_width,
                                    min_y + j * voxel_width,
                                    min_z + k * voxel_width)
                cube = Cube(midpoint, voxel_width / 2)
                if intersection_cube_sphere(cube, sphere):
                    voxel = Voxel(cube, [])
                    voxels.append(voxel)

    return voxels


def add_voxels_to_pdb(pdb_in, pocket_voxels, adjacent_voxels, pdb_out):

    with open(pdb_in, 'r') as pdb:
        with open(pdb_out, 'w') as out:
            for line in pdb:
                if line.strip() != "END":
                    out.write(line)

            for i, voxel in enumerate(pocket_voxels):

                pdb_line = f"HETATM{i + 1:5d} EMPT STP E{1:4d}    {voxel.cube.midpoint.x:8.3f}{voxel.cube.midpoint.y:8.3f}{voxel.cube.midpoint.z:8.3f}  0.00  0.00          Ve\n"
                out.write(pdb_line)
            if pocket_voxels:
                out.write("TER\n")

            for i, voxel in enumerate(adjacent_voxels):

                pdb_line = f"HETATM{i + 1:5d} {VOXEL_TYPE_TO_ABBREVIATION[voxel.type]} STP E{2:4d}    {voxel.cube.midpoint.x:8.3f}{voxel.cube.midpoint.y:8.3f}{voxel.cube.midpoint.z:8.3f}  0.00  0.00          Ve\n"
                out.write(pdb_line)
            if adjacent_voxels:
                out.write("TER\n")

            out.write("END\n")


def visualise_in_pymol(pdb_file, pdb_name, out_dir):
    cmd.load(pdb_file)

    cmd.hide('lines', 'resn STP')

    # show spheres, resn STP
    cmd.select("pocket", "resn STP and resi 1 and name EMPT")

    cmd.select("pocket_lining", "resn STP and resi 2")
    cmd.select("cysteines", "resn STP and resi 2 and name CYST")
    cmd.select("nitrogens", "resn STP and resi 2 and name POSI")
    cmd.select("oxygens", "resn STP and resi 2 and name NEGA")
    cmd.select("hydrophobic", "resn STP and resi 2 and name APOL")
    cmd.select("polar", "resn STP and resi 2 and name POL")
    cmd.select("metal", "resn STP and resi 2 and name MET")

    cmd.show("spheres", "pocket")
    cmd.set("sphere_scale", "0.3", "pocket")

    cmd.show("spheres", "pocket")
    cmd.set("sphere_scale", "0.3", "hydrophobic")
    cmd.show("spheres", "hydrophobic")
    cmd.set("sphere_scale", "0.3", "cysteines")
    cmd.show("spheres", "cysteines")
    cmd.set("sphere_scale", "0.3", "oxygens")
    cmd.show("spheres", "oxygens")
    cmd.set("sphere_scale", "0.3", "nitrogens")
    cmd.show("spheres", "nitrogens")
    cmd.set("sphere_scale", "0.3", "polar")
    cmd.show("spheres", "metal")
    cmd.set("sphere_scale", "0.3", "metal")

    cmd.color("white", "pocket")

    cmd.color("grey", "hydrophobic")
    cmd.color("red", "oxygens")
    cmd.color("blue", "nitrogens")
    cmd.color("yellow", "cysteines")
    cmd.color("purple", "polar")
    cmd.color("green", "metal")

    cmd.set("sphere_transparency", "0.1", "pocket")

    cmd.hide("sticks", "pocket")

    cmd.hide("sticks", "cysteines")
    cmd.hide("sticks", "nitrogens")
    cmd.hide("sticks", "oxygens")
    cmd.hide("sticks", "hydrophobic")
    cmd.hide("sticks", "polar")
    cmd.hide("sticks", "metal")

    cmd.save(os.path.join(out_dir, f"{pdb_name}.pse"))
    cmd.reinitialize()


def label_adjacent_voxels(pocket, all_voxels, voxel_width):
    involved_residues_sidechain = set()
    involved_residues_mainchain = set()
    full_voxels = []
    for voxel in all_voxels:
        if voxel.atoms:
            full_voxels.append(voxel)

    adjacency_matrix_pocket = make_adjacency_matrix(pocket, full_voxels, voxel_width)
    adjacent_voxels = set()

    for voxel in pocket:
        voxel.type = 'empty'
        voxel.pocket = True

        for adjacent_voxel in adjacency_matrix_pocket[voxel]:
            voxel.pocket = True
            counts = {'nitrogen': 0,
                      'oxygen': 0,
                      'sulphur': 0,
                      'hydrophobic': 0,
                      'metal': 0}
            total_count = 0
            for atom in adjacent_voxel.atoms:
                if atom.type not in {'C', 'CA', 'O', 'OXT'}:
                    involved_residues_sidechain.add(atom.residue)
                else:
                    involved_residues_mainchain.add(atom.residue)
                if atom.type in {'N', 'O', 'C', 'S', 'Mg'}:
                    total_count += 1
                if atom.type == 'N':
                    counts['nitrogen'] += 1
                elif atom.type == 'C':
                    counts['hydrophobic'] += 1
                elif atom.type == 'O':
                    counts['oxygen'] += 1
                elif atom.type == 'S':
                    counts['sulphur'] += 1
                elif atom.type == 'Mg':
                    counts['metal'] += 1

            if total_count == 0:
                adjacent_voxel.type = 'hydrophobic'
            elif float(counts['hydrophobic']) / float(total_count) > 0.50:
                adjacent_voxel.type = 'hydrophobic'
            else:
                best_count = 0
                best_candidates = []
                for atom_type, count in counts.items():
                    if atom_type != 'hydrophobic':
                        if count > best_count:
                            best_candidates = [atom_type]
                            best_count = count
                        elif count == best_count:
                            best_candidates.append(atom_type)
                if len(best_candidates) == 0:
                    adjacent_voxel.type = 'hydrophobic'
                elif len(best_candidates) == 1:
                    adjacent_voxel.type = best_candidates[0]
                else:
                    adjacent_voxel.type = 'polar'
            adjacent_voxels.add(adjacent_voxel)

    return adjacent_voxels


def pdb_to_voxel_representation(pdb_file: str, voxel_width: float,
                                distance_to_substrate: float,
                                substrate_location: Vector3D,
                                offset: Vector3D) -> List[Voxel]:
    voxel_grid = make_voxel_grid(voxel_width, distance_to_substrate, substrate_location, offset)
    atoms = atoms_from_pdb(pdb_file)
    apolar_atoms = {'C', 'H'}
    for voxel in voxel_grid:
        for atom in atoms:
            if intersection_cube_sphere(voxel.cube, atom.sphere):
                voxel.atoms.append(atom)

    for voxel in voxel_grid:
        for atom in voxel.atoms:
            if atom.type not in apolar_atoms:
                voxel.polar = True
                break

    return voxel_grid


def find_closest_voxel(location, voxels):
    min_distance = 10000.0
    closest_voxel = None
    for voxel in voxels:
        distance = location.euclidean_distance(voxel.cube.midpoint)
        if distance < min_distance:
            closest_voxel = voxel
            min_distance = distance

    return closest_voxel


def trim_pocket(pocket_voxels, voxel_width):
    origin_voxel = None

    for voxel in pocket_voxels:
        if voxel.cube.midpoint == SUBSTRATE_LOCATION:
            origin_voxel = voxel

    adjacency_matrix = make_adjacency_matrix(pocket_voxels, pocket_voxels, voxel_width)

    if origin_voxel not in adjacency_matrix:
        origin_voxel = find_closest_voxel(SUBSTRATE_LOCATION, list(adjacency_matrix.keys()))
        if origin_voxel is None:
            raise RuntimeError("No pocket around substrate location")

    pocket = set()
    pocket.add(origin_voxel)
    seen = set()
    pocket_length = 1

    while True:
        for voxel in list(pocket):
            if voxel not in seen:
                for neighbour in adjacency_matrix[voxel]:
                    pocket.add(neighbour)
                    seen.add(voxel)
        new_pocket_length = len(pocket)
        if new_pocket_length == pocket_length:
            break
        pocket_length = new_pocket_length

    return pocket


if __name__ == "__main__":
    run()

import os
from argparse import ArgumentParser

import torch
from torch.utils.data import Dataset

from paras.scripts.parsers.voxel import voxels_from_file, Voxel
from paras.scripts.math.shapes import Vector3D, Cube
from paras.scripts.parsers.parsers import parse_domain_list, parse_specificities, parse_substrate_list
from paras.scripts.data_processing.relabel_data import relabel_from_substrate_list
from paras.scripts.data_analysis.count_substrates import count_substrates_from_file

BETA_CARBON_LOCATION = Vector3D(32.833,  97.884,  32.711)
SUBSTRATE_LOCATION = BETA_CARBON_LOCATION

OFFSET = Vector3D(0.0, 0.0, 0.0)


def parse_arguments():
    parser = ArgumentParser()
    parser.add_argument('-r', type=float, default=10.0, help="Radius.")
    parser.add_argument('-w', type=float, default=1.2, help="Voxel width.")
    parser.add_argument('-i', type=str, required=True, help="Path to voxel file.")
    parser.add_argument('-o', type=str, required=True, help="Path to output_directory.")
    args = parser.parse_args()

    return args


def run():
    args = parse_arguments()
    if not os.path.exists(args.o):
        os.mkdir(args.o)

    for domain_name, grid in voxels_as_grid(args.i, radius=args.r, voxel_width=args.w):
        out_path = os.path.join(args.o, f"{domain_name}.txt")
        write_grid(grid, out_path)


def encode_specificity(specificity, all_specificities):
    target = torch.zeros(len(all_specificities))
    for spec in specificity:
        if spec in all_specificities:
            index = all_specificities.index(spec)
            target[index] = 1
    return torch.FloatTensor(target)


class VoxelDataset(Dataset):
    def __init__(self, data_dir, domain_list, specificities_file, substrate_file, grid_size=17):
        self.grid_size = grid_size
        self.data_dir = data_dir
        self.domain_list = parse_domain_list(domain_list)
        self.specificities = []
        self.weights = []

        domain_to_specificity = parse_specificities(specificities_file)
        domain_to_filtered = relabel_from_substrate_list(self.domain_list, domain_to_specificity, substrate_file)

        self.raw_specificities = []

        for domain, specificities in domain_to_specificity.items():
            self.raw_specificities += specificities

        self.specificity_counts = count_substrates_from_file(domain_list, specificities_file, substrate_file)

        self.all_specificities = parse_substrate_list(substrate_file)

        self.all_specificities.sort()

        for domain in self.domain_list:
            self.specificities.append(domain_to_filtered[domain])

    def __len__(self):
        return len(self.domain_list)

    def __getitem__(self, idx):
        file_name = os.path.join(self.data_dir, f"{self.domain_list[idx]}.txt")
        tensor = read_grid(file_name, x=self.grid_size, y=self.grid_size, z=self.grid_size)
        specificity = encode_specificity(self.specificities[idx], self.all_specificities)

        return tensor, specificity

    def get_weights(self):
        true_weights = []
        for specificity in self.all_specificities:
            nr_true = self.specificity_counts[specificity]
            nr_false = len(self) - nr_true
            true_weights.append(nr_false / nr_true)

        return true_weights


class VoxelSequenceDataset(VoxelDataset):
    def __init__(self, data_dir, domain_list, specificities_file, substrate_file):
        super().__init__(data_dir, domain_list, specificities_file, substrate_file)


class SiameseVoxelDataset(Dataset):
    def __init__(self, data_dir, domain_list, specificities_file, smiles_file):
        self.data_dir = data_dir
        self.domain_list = parse_domain_list(domain_list)
        self.specificities = []
        domain_to_specificity = parse_specificities(specificities_file)
        for domain in self.domain_list:
            self.specificities.append(domain_to_specificity[domain])

        self.similarity_matrix = self.get_similarity_matrix(smiles_file)

    def __getitem__(self, idx):
        pass

    def get_similarity_matrix(self, smiles_file):
        similarity_matrix = []

        for i, specificity_1 in enumerate(self.specificities):
            similarity_row = []
            for j, specificity_2 in enumerate(self.specificities):
                pass


def voxel_to_grid_position(voxel, min_x, min_y, min_z, voxel_width=1.2):
    i = round((voxel.cube.midpoint.x - (min_x + 0.5 * voxel_width)) / voxel_width)
    j = round((voxel.cube.midpoint.y - (min_y + 0.5 * voxel_width)) / voxel_width)
    k = round((voxel.cube.midpoint.z - (min_z + 0.5 * voxel_width)) / voxel_width)

    return i, j, k


def voxels_as_grid(voxel_file, radius=10.0, voxel_width=1.2):

    center = SUBSTRATE_LOCATION
    offset = OFFSET

    nr_voxels = int(radius / voxel_width)
    min_x = center.x - nr_voxels * voxel_width + offset.x
    min_y = center.y - nr_voxels * voxel_width + offset.y
    min_z = center.z - nr_voxels * voxel_width + offset.z

    for domain_name, voxels in voxels_from_file(voxel_file, pocket_only=False):

        grid = []

        for i in range(nr_voxels * 2):
            plane = []
            for j in range(nr_voxels * 2):
                row = []
                for k in range(nr_voxels * 2):
                    midpoint = Vector3D(min_x + i * voxel_width + 0.5 * voxel_width,
                                        min_y + j * voxel_width + 0.5 * voxel_width,
                                        min_z + k * voxel_width + 0.5 * voxel_width)
                    cube = Cube(midpoint, voxel_width / 2)
                    voxel = Voxel(cube, [])
                    row.append(voxel)
                plane.append(row)
            grid.append(plane)

        for voxel in voxels:
            i, j, k = voxel_to_grid_position(voxel, min_x, min_y, min_z, voxel_width=voxel_width)
            grid[i][j][k] = voxel

        yield domain_name, grid


def write_grid(grid, out_file):
    with open(out_file, 'w') as out:

        for i, plane in enumerate(grid):
            for j, row in enumerate(plane):
                for k, voxel in enumerate(row):
                    out.write(f"{i}\t{j}\t{k}\t")
                    attributes = voxel.vector.as_list(categories=["nr_aromatic",
                                                                  "nr_hydrophobic",
                                                                  "nr_pos_charge",
                                                                  "nr_neg_charge",
                                                                  "nr_pos_polar",
                                                                  "nr_neg_polar",
                                                                  "nr_sulfur"],
                                                      data_type='int')
                    out.write('\t'.join(map(str, attributes)))
                    out.write('\n')


def read_grid(grid_file, x=17, y=17, z=17, in_channels=7):
    empty_grid = []
    for a in range(in_channels):
        empty_grid.append([])
        for i in range(x):
            empty_grid[a].append([])
            for j in range(y):
                empty_grid[a][i].append([])
                for k in range(z):
                    empty_grid[a][i][j].append(0)

    with open(grid_file, 'r') as grid:
        for line in grid:
            voxel_data = line.split('\t')
            voxel_positions = list(map(int, voxel_data[:3]))
            attribute_data = list(map(int, voxel_data[3:]))
            for i, attribute in enumerate(attribute_data):
                empty_grid[i][voxel_positions[0]][voxel_positions[1]][voxel_positions[2]] = attribute

    tensor = torch.Tensor(empty_grid)
    return tensor


if __name__ == "__main__":
    run()

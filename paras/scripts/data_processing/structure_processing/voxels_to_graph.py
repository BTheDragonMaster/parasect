from torch_geometric.data import Data, InMemoryDataset
import torch
from typing import List, Tuple
from paras.scripts.parsers.voxel import Voxel, voxels_from_file
from paras.scripts.data_processing.structure_processing.voxel_adjacency_matrix import make_adjacency_matrix
import networkx as nx
from sys import argv
from dataclasses import astuple
from scipy.spatial.distance import euclidean


def get_distance(cloud_1, cloud_2):
    min_distance = 1000000000.0
    voxels = []
    for voxel_1 in cloud_1:
        for voxel_2 in cloud_2:
            distance = euclidean([voxel_1.cube.midpoint.x,
                                  voxel_1.cube.midpoint.y,
                                  voxel_1.cube.midpoint.z],
                                 [voxel_2.cube.midpoint.x,
                                  voxel_2.cube.midpoint.y,
                                  voxel_2.cube.midpoint.z])
            if distance < min_distance:
                min_distance = distance
                voxels = (voxel_1, voxel_2)

    assert voxels

    return min_distance, voxels


def connect_voxel_clouds(voxel_clouds):

    connections = []
    distances = []

    while len(voxel_clouds) > 1:
        min_distance = 1000000000.0
        cloud_indices = []
        connection = []

        for i, cloud_1 in enumerate(voxel_clouds):
            for j, cloud_2 in enumerate(voxel_clouds):
                if i > j:
                    distance, voxels = get_distance(cloud_1, cloud_2)
                    if distance < min_distance:
                        min_distance = distance
                        cloud_indices = [i, j]
                        connection = voxels
        connections.append(connection)
        distances.append(min_distance)
        # as i is always larger than j, remove them in this order
        cloud_1 = voxel_clouds.pop(cloud_indices[0])
        cloud_2 = voxel_clouds.pop(cloud_indices[1])
        new_cloud = cloud_1 + cloud_2
        voxel_clouds.append(new_cloud)

    return connections, distances


def voxels_to_edges(voxels):
    source_nodes = []
    target_nodes = []

    idx_to_voxel = {}
    voxel_to_idx = {}

    adjacency_matrix = make_adjacency_matrix(voxels, voxels, 1.2)
    nx_graph = nx.from_dict_of_lists(adjacency_matrix)
    voxel_clouds = [list(c) for c in nx.connected_components(nx_graph)]
    connections, distances = connect_voxel_clouds(voxel_clouds)

    for voxel_1, voxel_2 in connections:
        adjacency_matrix[voxel_1].append(voxel_2)
        adjacency_matrix[voxel_2].append(voxel_1)

    sorted_voxels = sorted(adjacency_matrix.keys(), key=lambda vox: (vox.cube.midpoint.x,
                                                                     vox.cube.midpoint.y,
                                                                     vox.cube.midpoint.z))
    for i, voxel in enumerate(sorted_voxels):
        idx_to_voxel[i] = voxel
        voxel_to_idx[voxel] = i

    edge_features = []

    for voxel, neighbours in adjacency_matrix.items():
        for neighbour in neighbours:
            source_nodes.append(voxel_to_idx[voxel])
            target_nodes.append(voxel_to_idx[neighbour])
            if (voxel, neighbour) in connections:
                edge_weight = distances[connections.index((voxel, neighbour))]
            elif (neighbour, voxel) in connections:
                edge_weight = distances[connections.index((neighbour, voxel))]
            else:
                edge_weight = 1.0

            edge_features.append(edge_weight)

    return source_nodes, target_nodes, edge_features, sorted_voxels


def voxels_to_graph(voxels) -> Tuple[torch.tensor, torch.tensor, torch.tensor]:
    """
    """
    source_nodes, target_nodes, edge_features, sorted_voxels = voxels_to_edges(voxels)
    edge_index = torch.tensor([source_nodes, target_nodes], dtype=torch.long)

    node_features = []
    for voxel in sorted_voxels:
        node_features.append([voxel.cube.midpoint.x,
                              voxel.cube.midpoint.y,
                              voxel.cube.midpoint.z] + list(astuple(voxel.vector)))

    node_features = torch.tensor(node_features, dtype=torch.float)
    edge_features = torch.tensor(edge_features, dtype=torch.float)

    return node_features, edge_index, edge_features


class EnzymeData(Data):

    def __init__(self, voxels: List[Voxel]) -> None:
        node_features, edge_index, edge_features = voxels_to_graph(voxels)
        super().__init__(x=node_features, edge_index=edge_index, edge_attr=edge_features)


class EnzymeDataset(InMemoryDataset):
    """
    Create a dataset of enzyme active sites in graph format.
    """

    def __init__(self, voxel_file: str) -> None:
        """
        Initialize an enzyme dataset.

        Arguments
        ---------
        voxel_file: str, path to file containing voxel data
        """
        super().__init__(".", None, None, None)
        data_list = []
        idx_to_domain = {}
        for i, domain_name_and_voxels in enumerate(voxels_from_file(voxel_file)):
            domain_name, voxels = domain_name_and_voxels
            enzyme_data = EnzymeData(voxels)
            data_list.append(enzyme_data)
            idx_to_domain[i] = domain_name

        self.data, self.slices = self.collate(data_list)


def test(voxel_file):
    for i, domain_name_and_voxels in enumerate(voxels_from_file(voxel_file)):
        domain_name, voxels = domain_name_and_voxels
        enzyme_data = EnzymeData(voxels)
        break
    print("Done")


if __name__ == "__main__":
    test(argv[1])


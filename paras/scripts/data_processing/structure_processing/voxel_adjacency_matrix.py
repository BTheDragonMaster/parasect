from math import isclose


def make_adjacency_matrix(voxels_1, voxels_2, voxel_width):
    voxel_to_neighbours = {}
    for i, voxel_1 in enumerate(voxels_1):
        if voxel_1 not in voxel_to_neighbours:
            voxel_to_neighbours[voxel_1] = []
        for j, voxel_2 in enumerate(voxels_2):
            if i != j:
                if isclose(abs(voxel_1.cube.midpoint.x - voxel_2.cube.midpoint.x), voxel_width, rel_tol=0.00005) \
                        and isclose(voxel_1.cube.midpoint.y - voxel_2.cube.midpoint.y, 0.0, rel_tol=0.00005) \
                        and isclose(voxel_1.cube.midpoint.z - voxel_2.cube.midpoint.z, 0.0, rel_tol=0.00005):
                    voxel_to_neighbours[voxel_1].append(voxel_2)
                elif isclose(voxel_1.cube.midpoint.x - voxel_2.cube.midpoint.x, 0.0, rel_tol=0.00005) \
                        and isclose(abs(voxel_1.cube.midpoint.y - voxel_2.cube.midpoint.y), voxel_width, rel_tol=0.00005) \
                        and isclose(voxel_1.cube.midpoint.z - voxel_2.cube.midpoint.z, 0.0, rel_tol=0.00005):
                    voxel_to_neighbours[voxel_1].append(voxel_2)
                elif isclose(voxel_1.cube.midpoint.x - voxel_2.cube.midpoint.x, 0.0, rel_tol=0.00005) \
                        and isclose(voxel_1.cube.midpoint.y - voxel_2.cube.midpoint.y, 0.0, rel_tol=0.00005) \
                        and isclose(abs(voxel_1.cube.midpoint.z - voxel_2.cube.midpoint.z), voxel_width, rel_tol=0.00005):
                    voxel_to_neighbours[voxel_1].append(voxel_2)
    return voxel_to_neighbours

from paras.scripts.math.shapes import Cube, Sphere, Vector3D
import math


def intersection_cube_sphere(voxel: Cube, sphere: Sphere):
    sphere_x = sphere.center.x - voxel.midpoint.x
    sphere_y = sphere.center.y - voxel.midpoint.y
    sphere_z = sphere.center.z - voxel.midpoint.z
    scaling_factor = 1.0 / voxel.radius

    sphere_x *= scaling_factor
    sphere_y *= scaling_factor
    sphere_z *= scaling_factor
    sphere_radius = sphere.radius * scaling_factor

    distance = math.sqrt(max([0, abs(sphere_x) - 1])**2
                         + max([0, abs(sphere_y) - 1])**2
                         + max([0, abs(sphere_z) - 1])**2)

    if distance < sphere_radius:
        return True
    else:
        return False

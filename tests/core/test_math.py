# -*- coding: utf-8 -*-

"""Unit tests for the math module."""

import math

import parasect.core.math as pcm


def test_point3d_construction():
    """Test the construction of a Point3D object."""
    point = pcm.Point3D(1.0, 2.0, 3.0)
    if (point.x != 1.0) or (point.y != 2.0) or (point.z != 3.0):
        raise


def test_point3d_construction_invalid_x_str_input():
    """Test the construction of a Point3D object with an invalid str x-coordinate."""
    try:
        pcm.Point3D("1.0", 2.0, 3.0)
    except TypeError:
        pass
    else:
        raise


def test_point3d_construction_invalid_x_int_input():
    """Test the construction of a Point3D object with an invalid int x-coordinate."""
    try:
        pcm.Point3D(1, 2.0, 3.0)
    except TypeError:
        pass
    else:
        raise


def test_point3d_construction_invalid_y__str_input():
    """Test the construction of a Point3D object with an invalid str y-coordinate."""
    try:
        pcm.Point3D(1.0, "2.0", 3.0)
    except TypeError:
        pass
    else:
        raise


def test_point3d_construction_invalid_y_int_input():
    """Test the construction of a Point3D object with an invalid int y-coordinate."""
    try:
        pcm.Point3D(1.0, 2, 3.0)
    except TypeError:
        pass
    else:
        raise


def test_point3d_construction_invalid_z_str_input():
    """Test the construction of a Point3D object with an invalid str z-coordinate."""
    try:
        pcm.Point3D(1.0, 2.0, "3.0")
    except TypeError:
        pass
    else:
        raise


def test_point3d_construction_invalid_z_int_input():
    """Test the construction of a Point3D object with an invalid int z-coordinate."""
    try:
        pcm.Point3D(1.0, 2.0, 3)
    except TypeError:
        pass
    else:
        raise


def test_point3d_equality():
    """Test the equality of two Point3D objects."""
    point1 = pcm.Point3D(1.0, 2.0, 3.0)
    point2 = pcm.Point3D(1.0, 2.0, 3.0)
    point3 = pcm.Point3D(1.0, 2.0, 3.00001)
    point4 = pcm.Point3D(2.0, 3.0, 4.0)

    if point1 != point2 or point1 != point3 or point1 == point4:
        raise


def test_point3d_hash():
    """Test the hash value of a Point3D object."""
    point1 = pcm.Point3D(1.0, 2.0, 3.0)
    point2 = pcm.Point3D(1.0, 2.0, 3.0)
    point3 = pcm.Point3D(1.0, 2.0, 3.00005)
    point4 = pcm.Point3D(1.0, 2.0, 3.00006)
    point5 = pcm.Point3D(2.0, 3.0, 4.0)

    if (
        hash(point1) != hash(point2)
        or hash(point1) != hash(point3)
        or hash(point1) == hash(point4)
        or hash(point1) == hash(point5)
    ):
        raise


def test_point3d_representation():
    """Test the string representation of a Point3D object."""
    point = pcm.Point3D(1.0, 2.0, 3.0)
    if repr(point) != "Point3D(x=1.0, y=2.0, z=3.0)":
        raise


def test_pont3d_calc_euclidean_distance():
    """Test the calculation of the Euclidean distance between two Point3D objects."""
    point1 = pcm.Point3D(1.0, 2.0, 3.0)
    point2 = pcm.Point3D(4.0, 6.0, 8.0)
    calc_dist = point1.calc_euclidean_distance(point2)
    if not math.isclose(calc_dist, 7.07, rel_tol=0.01):
        raise


def test_pont3d_calc_euclidean_distance_same_point():
    """Test the calculation of the Euclidean distance between two identical Point3D objects."""
    point1 = pcm.Point3D(1.0, 2.0, 3.0)
    calc_dist = point1.calc_euclidean_distance(point1)
    if not math.isclose(calc_dist, 0.0, rel_tol=0.01):
        raise


def test_pont3d_calc_euclidean_distance_invalid_other_type():
    """Test the calculation of the Euclidean distance between a Point3D object and an invalid type."""
    point1 = pcm.Point3D(1.0, 2.0, 3.0)
    try:
        point1.calc_euclidean_distance([1.0, 2.0, 3.0])
    except TypeError:
        pass
    else:
        raise


def test_sphere_construction():
    """Test the construction of a Sphere object."""
    center = pcm.Point3D(1.0, 2.0, 3.0)
    sphere = pcm.Sphere(center, 4.0)
    if (sphere.center != center) or (sphere.radius != 4.0):
        raise


def test_cube_construction():
    """Test the construction of a Cube object."""
    center = pcm.Point3D(1.0, 2.0, 3.0)
    cube = pcm.Cube(center, 4.0)
    if (cube.center != center) or (cube.radius != 4.0):
        raise


def test_cube_and_sphere_intersect_intersection():
    """Test the intersection of a Cube and a Sphere."""
    cube = pcm.Cube(pcm.Point3D(0.0, 0.0, 0.0), 1.0)
    sphere = pcm.Sphere(pcm.Point3D(0.0, 0.0, 0.0), 1.0)
    if not pcm.calc_if_cube_and_sphere_intersect(cube, sphere):
        raise


def test_cube_and_sphere_intersect_no_intersection():
    """Test the intersection of a Cube and a Sphere."""
    cube = pcm.Cube(pcm.Point3D(0.0, 0.0, 0.0), 1.0)
    sphere = pcm.Sphere(pcm.Point3D(2.0, 0.0, 0.0), 1.0)
    if pcm.calc_if_cube_and_sphere_intersect(cube, sphere):
        raise


def test_cube_and_sphere_interect_invalid_input():
    """Test the intersection of a Cube and a Sphere with invalid input."""
    cube = pcm.Cube(pcm.Point3D(0.0, 0.0, 0.0), 1.0)
    sphere = pcm.Sphere(pcm.Point3D(0.0, 0.0, 0.0), 1.0)
    try:
        pcm.calc_if_cube_and_sphere_intersect(cube, sphere.center)
    except TypeError:
        pass
    else:
        raise
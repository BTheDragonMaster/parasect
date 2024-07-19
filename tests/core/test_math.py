# -*- coding: utf-8 -*-

"""Unit tests for the math module."""

import math

import parasect.core.math as pcm


def test_point3d_construction():
    """Test the construction of a Point3D object."""
    point = pcm.Point3D(1.0, 2.0, 3.0)
    assert point.x == 1.0
    assert point.y == 2.0
    assert point.z == 3.0


def test_point3d_equality():
    """Test the equality of two Point3D objects."""
    point1 = pcm.Point3D(1.0, 2.0, 3.0)
    point2 = pcm.Point3D(1.0, 2.0, 3.0)
    point3 = pcm.Point3D(1.0, 2.0, 3.00001)
    point4 = pcm.Point3D(2.0, 3.0, 4.0)

    assert point1 == point2
    assert point1 == point3
    assert point1 != point4


def test_point3d_hash():
    """Test the hash value of a Point3D object."""
    point1 = pcm.Point3D(1.0, 2.0, 3.0)
    point2 = pcm.Point3D(1.0, 2.0, 3.0)
    point3 = pcm.Point3D(1.0, 2.0, 3.00005)
    point4 = pcm.Point3D(1.0, 2.0, 3.00006)
    point5 = pcm.Point3D(2.0, 3.0, 4.0)

    assert hash(point1) == hash(point2)
    assert hash(point1) == hash(point3)
    assert hash(point1) != hash(point4)
    assert hash(point1) != hash(point5)


def test_point3d_representation():
    """Test the string representation of a Point3D object."""
    point = pcm.Point3D(1.0, 2.0, 3.0)
    assert repr(point) == "Point3D(x=1.0, y=2.0, z=3.0)"


def test_pont3d_calc_euclidean_distance():
    """Test the calculation of the Euclidean distance between two Point3D objects."""
    point1 = pcm.Point3D(1.0, 2.0, 3.0)
    point2 = pcm.Point3D(4.0, 6.0, 8.0)
    calc_dist = point1.calc_euclidean_distance(point2)
    assert math.isclose(calc_dist, 7.07, rel_tol=0.01)


def test_pont3d_calc_euclidean_distance_same_point():
    """Test the calculation of the Euclidean distance between two identical Point3D objects."""
    point1 = pcm.Point3D(1.0, 2.0, 3.0)
    calc_dist = point1.calc_euclidean_distance(point1)
    assert math.isclose(calc_dist, 0.0, rel_tol=0.01)


def test_sphere_construction():
    """Test the construction of a Sphere object."""
    center = pcm.Point3D(1.0, 2.0, 3.0)
    sphere = pcm.Sphere(center, 4.0)
    assert sphere.center == center
    assert sphere.radius == 4.0


def test_cube_construction():
    """Test the construction of a Cube object."""
    center = pcm.Point3D(1.0, 2.0, 3.0)
    cube = pcm.Cube(center, 4.0)
    assert cube.center == center
    assert cube.radius == 4.0


def test_cube_and_sphere_intersect_intersection():
    """Test the intersection of a Cube and a Sphere."""
    cube = pcm.Cube(pcm.Point3D(0.0, 0.0, 0.0), 1.0)
    sphere = pcm.Sphere(pcm.Point3D(0.0, 0.0, 0.0), 1.0)
    assert pcm.calc_if_cube_and_sphere_intersect(cube, sphere)


def test_cube_and_sphere_intersect_no_intersection():
    """Test the intersection of a Cube and a Sphere."""
    cube = pcm.Cube(pcm.Point3D(0.0, 0.0, 0.0), 1.0)
    sphere = pcm.Sphere(pcm.Point3D(2.0, 0.0, 0.0), 1.0)
    assert not pcm.calc_if_cube_and_sphere_intersect(cube, sphere)

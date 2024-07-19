# -*- coding: utf-8 -*-

"""The math module contains classes and functions for mathematical operations."""

import math
import typing as ty


class Point3D:
    """The Point3D class represents a 3D vector."""

    def __init__(self, x: float, y: float, z: float) -> None:
        """Initialize the 3D point.

        :param x: The x-coordinate of the vector.
        :type x: float
        :param y: The y-coordinate of the vector.
        :type y: float
        :param z: The z-coordinate of the vector.
        :type z: float
        :raises TypeError: If the coordinates are not floats.
        """
        if not isinstance(x, (float)):
            raise TypeError(f"Point3D expects float for x, got {type(x)}")

        if not isinstance(y, (float)):
            raise TypeError(f"Point3D expects float for y, got {type(y)}")

        if not isinstance(z, (float)):
            raise TypeError(f"Point3D expects float for z, got {type(z)}")

        self.x = x
        self.y = y
        self.z = z

    def __eq__(self, other: ty.Any) -> bool:
        """Check if two points are equal.

        :param other: The other point to compare.
        :type other: ty.Any
        :return: True if the points are equal, otherwise False.
        :rtype: bool
        """
        if type(self) is not type(other):
            return False

        if (
            math.isclose(self.x, other.x, rel_tol=0.00005)
            and math.isclose(self.y, other.y, rel_tol=0.00005)
            and math.isclose(self.z, other.z, rel_tol=0.00005)
        ):
            return True
        else:
            return False

    def __hash__(self) -> int:
        """Return the hash value of the point.

        :return: The hash value of the point.
        :rtype: int
        """
        return hash((round(self.x, 4), round(self.y, 4), round(self.z, 4)))

    def __repr__(self) -> str:
        """Return the string representation of the point.

        :return: The string representation of the point.
        :rtype: str
        """
        return f"Point3D(x={(round(self.x, 1))}, y={(round(self.y, 1))}, z={(round(self.z, 1))})"

    def calc_euclidean_distance(self, other: "Point3D") -> float:
        """Calculate the Euclidean distance between two points.

        :param other: The other point to calculate the distance.
        :type other: Point3D
        :return: The Euclidean distance between the two points.
        :rtype: float
        :raises TypeError: If the other point is not a Point3D object.
        """
        if type(self) is not type(other):
            raise TypeError(f"Point3D expects Point3D for other, got {type(other)}")

        return math.sqrt((self.x - other.x) ** 2 + (self.y - other.y) ** 2 + (self.z - other.z) ** 2)

    def to_tuple(self) -> ty.Tuple[float, float, float]:
        """Return the point as a tuple.

        :return: The point as a tuple.
        :rtype: ty.Tuple[float, float, float]
        """
        return (self.x, self.y, self.z)


class Sphere:
    """The Sphere class represents a 3D sphere."""

    def __init__(self, center: Point3D, radius: float) -> None:
        """Initialize the sphere.

        :param center: The center of the sphere.
        :type center: Point3D
        :param radius: The radius of the sphere.
        :type radius: float
        """
        self.center = center
        self.radius = radius


class Cube:
    """The Cube class represents a 3D cube without defined rotation."""

    def __init__(self, center: Point3D, radius: float) -> None:
        """Initialize the cube.

        :param center: The center of the cube, equidistant from all faces.
        :type center: Point3D
        :param radius: The radius of the cube.
        :type radius: float
        """
        self.center = center
        self.radius = radius


def calc_if_cube_and_sphere_intersect(voxel: Cube, sphere: Sphere) -> bool:
    """Calculate if a cube and a sphere intersect.

    :param voxel: The cube to check for intersection.
    :type voxel: Cube
    :param sphere: The sphere to check for intersection.
    :type sphere: Sphere
    :return: True if the cube and sphere intersect, otherwise False.
    :rtype: bool
    :raises TypeError: If the voxel is not a Cube object or the sphere is not a Sphere object.
    """
    if type(voxel) is not Cube:
        raise TypeError(f"calc_if_cube_and_sphere_intersect expects Cube for voxel, got {type(voxel)}")

    if type(sphere) is not Sphere:
        raise TypeError(f"calc_if_cube_and_sphere_intersect expects Sphere for sphere, got {type(sphere)}")

    sphere_x = sphere.center.x - voxel.center.x
    sphere_y = sphere.center.y - voxel.center.y
    sphere_z = sphere.center.z - voxel.center.z

    scaling_factor = 1.0 / voxel.radius

    sphere_x *= scaling_factor
    sphere_y *= scaling_factor
    sphere_z *= scaling_factor
    sphere_radius = sphere.radius * scaling_factor

    distance = math.sqrt(
        max([0, abs(sphere_x) - 1]) ** 2 + max([0, abs(sphere_y) - 1]) ** 2 + max([0, abs(sphere_z) - 1]) ** 2
    )

    if distance < sphere_radius:
        return True
    else:
        return False

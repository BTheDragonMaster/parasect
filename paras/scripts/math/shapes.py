from dataclasses import dataclass
from math import isclose
from scipy.spatial import distance


@dataclass
class Vector3D:
    x: float
    y: float
    z: float

    def __eq__(self, other):
        if isclose(self.x, other.x, rel_tol=0.00005) and isclose(self.y, other.y, rel_tol=0.00005) and \
                isclose(self.z, other.z, rel_tol=0.00005):
            return True
        else:
            return False

    def __hash__(self):
        return hash((round(self.x, 5), round(self.y, 5), round(self.z, 5)))

    def __repr__(self):
        return f"{(round(self.x, 1))}_{(round(self.y, 1))}_{(round(self.z, 1))}"

    def euclidean_distance(self, other):
        return distance.euclidean([self.x, self.y, self.z], [other.x, other.y, other.z])

    def to_list(self):
        return [self.x, self.y, self.z]


@dataclass
class Sphere:
    center: Vector3D
    radius: float


@dataclass
class Cube:
    midpoint: Vector3D
    radius: float






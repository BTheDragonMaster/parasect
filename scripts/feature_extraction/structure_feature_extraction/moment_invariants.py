# Taken from https://github.com/TurtleTools/geometricus/blob/5d61047637ba980ca03a145fa755dc560f2e548b/geometricus/moment_utility.py#L1422

import typing as ty
from dataclasses import dataclass
from enum import Enum

import numba as nb
import numpy as np


@nb.njit
def nb_mean_axis_0(array: np.ndarray) -> np.ndarray:
    """
    Same as np.mean(array, axis=0) but njitted
    """
    mean_array = np.zeros(array.shape[1])
    for i in range(array.shape[1]):
        mean_array[i] = np.mean(array[:, i])
    return mean_array


@dataclass
class MomentInfo:
    moment_function: ty.Callable[[int, int, int, np.ndarray, np.ndarray], float]
    mu_arguments: ty.List[ty.Tuple[int, int, int]]


@nb.njit(cache=False)
def mu(p, q, r, coords, centroid):
    """
    Central moment
    """
    return np.sum(
        ((coords[:, 0] - centroid[0]) ** p) * ((coords[:, 1] - centroid[1]) ** q) * ((coords[:, 2] - centroid[2]) ** r)
    )


@nb.njit
def O_3(mu_200, mu_020, mu_002):
    return mu_200 + mu_020 + mu_002


@nb.njit
def O_4(mu_200, mu_020, mu_002, mu_110, mu_101, mu_011):
    return (
        mu_200 * mu_020 * mu_002
        + 2 * mu_110 * mu_101 * mu_011
        - mu_002 * mu_110**2
        - mu_020 * mu_101**2
        - mu_200 * mu_011**2
    )


@nb.njit
def O_5(mu_200, mu_020, mu_002, mu_110, mu_101, mu_011):
    return mu_200 * mu_020 + mu_200 * mu_002 + mu_020 * mu_002 - mu_110**2 - mu_101**2 - mu_011**2


@nb.njit
def F(
    mu_201,
    mu_021,
    mu_210,
    mu_300,
    mu_111,
    mu_012,
    mu_003,
    mu_030,
    mu_102,
    mu_120,
):
    return (
        mu_003**2
        + 6 * mu_012**2
        + 6 * mu_021**2
        + mu_030**2
        + 6 * mu_102**2
        + 15 * mu_111**2
        - 3 * mu_102 * mu_120
        + 6 * mu_120**2
        - 3 * mu_021 * mu_201
        + 6 * mu_201**2
        - 3 * mu_003 * (mu_021 + mu_201)
        - 3 * mu_030 * mu_210
        + 6 * mu_210**2
        - 3 * mu_012 * (mu_030 + mu_210)
        - 3 * mu_102 * mu_300
        - 3 * mu_120 * mu_300
        + mu_300**2
    )


def make_formula(name, formula_string):
    """
    Generate code from one of the formula in Appendix 4A of "2D and 3D Image Analysis by Moments"

    Parameters
    ----------
    name
        moment_name
    formula_string
        formula copy-pasted from PDF
    """
    formula = []
    mu_types = set()
    for x in formula_string.split("+"):
        f = ""
        parts = x.split("𝜇")
        if len(parts[0]):
            f += f"{parts[0]}"
        else:
            f += "1"
        for group in parts[1:]:
            f += " * mu_"
            if len(group) == 4:
                power = group[0]
                mul = group[1:]
                f += f"{mul} ** {power}"
            else:
                mul = group
                f += f"{mul}"
            mu_types.add(f"mu_{mul}")
        formula.append(f)
    formula = " + ".join(formula)
    mu_arguments = list(tuple(int(c) for c in m.split("_")[1]) for m in mu_types)
    mu_types = ", ".join(mu_types)
    print(f"{name} = MomentInfo({name}, {mu_arguments})")
    print(f"@nb.njit\ndef {name}({mu_types}):\n    return {formula}")


@nb.njit
def phi_2(mu_020, mu_011, mu_110, mu_200, mu_002, mu_101):
    return mu_200**2 + mu_020**2 + mu_002**2 + 2 * mu_110**2 + 2 * mu_101**2 + 2 * mu_011**2


@nb.njit
def phi_3(mu_020, mu_011, mu_110, mu_200, mu_002, mu_101):
    return (
        mu_200**3
        + 3 * mu_200 * mu_110**2
        + 3 * mu_200 * mu_101**2
        + 3 * mu_110**2 * mu_020
        + 3 * mu_101**2 * mu_002
        + mu_020**3
        + 3 * mu_020 * mu_011**2
        + 3 * mu_011**2 * mu_002
        + mu_002**3
        + 6 * mu_110 * mu_101 * mu_011
    )


@nb.njit
def phi_4(
    mu_030,
    mu_021,
    mu_120,
    mu_003,
    mu_111,
    mu_201,
    mu_102,
    mu_210,
    mu_012,
    mu_300,
):
    return (
        mu_300**2
        + mu_030**2
        + mu_003**2
        + 3 * mu_210**2
        + 3 * mu_201**2
        + 3 * mu_120**2
        + 3 * mu_102**2
        + 3 * mu_021**2
        + 3 * mu_012**2
        + 6 * mu_111**2
    )


@nb.njit
def phi_5(mu_030, mu_021, mu_120, mu_003, mu_201, mu_102, mu_210, mu_012, mu_300):
    return (
        mu_300**2
        + 2 * mu_300 * mu_120
        + 2 * mu_300 * mu_102
        + 2 * mu_210 * mu_030
        + 2 * mu_201 * mu_003
        + mu_030**2
        + 2 * mu_030 * mu_012
        + 2 * mu_021 * mu_003
        + mu_003**2
        + mu_210**2
        + 2 * mu_210 * mu_012
        + mu_201**2
        + 2 * mu_201 * mu_021
        + mu_120**2
        + 2 * mu_120 * mu_102
        + mu_102**2
        + mu_021**2
        + mu_012**2
    )


@nb.njit
def phi_6(
    mu_030,
    mu_021,
    mu_120,
    mu_003,
    mu_111,
    mu_201,
    mu_102,
    mu_210,
    mu_012,
    mu_300,
):
    return (
        1 * mu_300**4
        + 6 * mu_300**2 * mu_210**2
        + 6 * mu_300**2 * mu_201**2
        + 2 * mu_300**2 * mu_120**2
        + 4 * mu_300**2 * mu_111**2
        + 2 * mu_300**2 * mu_102**2
        + 8 * mu_300 * mu_210**2 * mu_120
        + 16 * mu_300 * mu_210 * mu_201 * mu_111
        + 4 * mu_300 * mu_210 * mu_120 * mu_030
        + 8 * mu_300 * mu_210 * mu_111 * mu_021
        + 4 * mu_300 * mu_210 * mu_102 * mu_012
        + 8 * mu_300 * mu_201**2 * mu_102
        + 4 * mu_300 * mu_201 * mu_120 * mu_021
        + 8 * mu_300 * mu_201 * mu_111 * mu_012
        + 4 * mu_300 * mu_201 * mu_102 * mu_003
        + 2 * mu_210**2 * mu_030**2
        + 4 * mu_210 * mu_201 * mu_030 * mu_021
        + 4 * mu_210 * mu_201 * mu_012 * mu_003
        + 8 * mu_210 * mu_120**2 * mu_030
        + 8 * mu_210 * mu_111 * mu_102 * mu_003
        + 2 * mu_201**2 * mu_003**2
        + 8 * mu_201 * mu_120 * mu_111 * mu_030
        + 8 * mu_201 * mu_102**2 * mu_003
        + 6 * mu_120**2 * mu_030**2
        + 16 * mu_120 * mu_111 * mu_030 * mu_021
        + 8 * mu_120 * mu_111 * mu_012 * mu_003
        + 4 * mu_120 * mu_102 * mu_030 * mu_012
        + 4 * mu_120 * mu_102 * mu_021 * mu_003
        + 4 * mu_111**2 * mu_030**2
        + 4 * mu_111**2 * mu_003**2
        + 8 * mu_111 * mu_102 * mu_030 * mu_021
        + 16 * mu_111 * mu_102 * mu_012 * mu_003
        + 6 * mu_102**2 * mu_003**2
        + 1 * mu_030**4
        + 6 * mu_030**2 * mu_021**2
        + 2 * mu_030**2 * mu_012**2
        + 8 * mu_030 * mu_021**2 * mu_012
        + 4 * mu_030 * mu_021 * mu_012 * mu_003
        + 2 * mu_021**2 * mu_003**2
        + 8 * mu_021 * mu_012**2 * mu_003
        + 6 * mu_012**2 * mu_003**2
        + 1 * mu_003**4
        + 5 * mu_210**4
        + 10 * mu_210**2 * mu_201**2
        + 16 * mu_210**2 * mu_120**2
        + 20 * mu_210**2 * mu_111**2
        + 4 * mu_210**2 * mu_102**2
        + 4 * mu_210**2 * mu_021**2
        + 2 * mu_210**2 * mu_012**2
        + 24 * mu_210 * mu_201 * mu_120 * mu_111
        + 24 * mu_210 * mu_201 * mu_111 * mu_102
        + 8 * mu_210 * mu_201 * mu_021 * mu_012
        + 24 * mu_210 * mu_120 * mu_111 * mu_021
        + 8 * mu_210 * mu_120 * mu_102 * mu_012
        + 16 * mu_210 * mu_111**2 * mu_012
        + 5 * mu_201**4
        + 4 * mu_201**2 * mu_120**2
        + 20 * mu_201**2 * mu_111**2
        + 16 * mu_201**2 * mu_102**2
        + 2 * mu_201**2 * mu_021**2
        + 4 * mu_201**2 * mu_012**2
        + 8 * mu_201 * mu_120 * mu_102 * mu_021
        + 16 * mu_201 * mu_111**2 * mu_021
        + 24 * mu_201 * mu_111 * mu_102 * mu_012
        + 5 * mu_120**4
        + 20 * mu_120**2 * mu_111**2
        + 2 * mu_120**2 * mu_102**2
        + 10 * mu_120**2 * mu_021**2
        + 4 * mu_120**2 * mu_012**2
        + 16 * mu_120 * mu_111**2 * mu_102
        + 24 * mu_120 * mu_111 * mu_021 * mu_012
        + 20 * mu_111**2 * mu_102**2
        + 20 * mu_111**2 * mu_021**2
        + 20 * mu_111**2 * mu_012**2
        + 24 * mu_111 * mu_102 * mu_021 * mu_012
        + 5 * mu_102**4
        + 4 * mu_102**2 * mu_021**2
        + 10 * mu_102**2 * mu_012**2
        + 5 * mu_021**4
        + 16 * mu_021**2 * mu_012**2
        + 5 * mu_012**4
        + 12 * mu_111**4
    )


@nb.njit
def phi_7(
    mu_030,
    mu_021,
    mu_120,
    mu_003,
    mu_111,
    mu_201,
    mu_102,
    mu_210,
    mu_012,
    mu_300,
):
    return (
        1 * mu_300**4
        + 1 * mu_300**3 * mu_120
        + 1 * mu_300**3 * mu_102
        + 5 * mu_300**2 * mu_210**2
        + 1 * mu_300**2 * mu_210 * mu_030
        + 1 * mu_300**2 * mu_210 * mu_012
        + 5 * mu_300**2 * mu_201**2
        + 1 * mu_300**2 * mu_201 * mu_021
        + 1 * mu_300**2 * mu_201 * mu_003
        + 1 * mu_300**2 * mu_120**2
        + 2 * mu_300**2 * mu_111**2
        + 1 * mu_300**2 * mu_102**2
        + 11 * mu_300 * mu_210**2 * mu_120
        + 4 * mu_300 * mu_210**2 * mu_102
        + 14 * mu_300 * mu_210 * mu_201 * mu_111
        + 4 * mu_300 * mu_210 * mu_120 * mu_030
        + 2 * mu_300 * mu_210 * mu_120 * mu_012
        + 6 * mu_300 * mu_210 * mu_111 * mu_021
        + 2 * mu_300 * mu_210 * mu_111 * mu_003
        + 2 * mu_300 * mu_210 * mu_102 * mu_012
        + 4 * mu_300 * mu_201**2 * mu_120
        + 11 * mu_300 * mu_201**2 * mu_102
        + 2 * mu_300 * mu_201 * mu_120 * mu_021
        + 2 * mu_300 * mu_201 * mu_111 * mu_030
        + 6 * mu_300 * mu_201 * mu_111 * mu_012
        + 2 * mu_300 * mu_201 * mu_102 * mu_021
        + 4 * mu_300 * mu_201 * mu_102 * mu_003
        + 3 * mu_300 * mu_120**3
        + 1 * mu_300 * mu_120**2 * mu_102
        + 8 * mu_300 * mu_120 * mu_111**2
        + 1 * mu_300 * mu_120 * mu_102**2
        + 1 * mu_300 * mu_120 * mu_030**2
        + 2 * mu_300 * mu_120 * mu_021**2
        + 1 * mu_300 * mu_120 * mu_012**2
        + 8 * mu_300 * mu_111**2 * mu_102
        + 2 * mu_300 * mu_111 * mu_030 * mu_021
        + 4 * mu_300 * mu_111 * mu_021 * mu_012
        + 2 * mu_300 * mu_111 * mu_012 * mu_003
        + 3 * mu_300 * mu_102**3
        + 1 * mu_300 * mu_102 * mu_021**2
        + 2 * mu_300 * mu_102 * mu_012**2
        + 1 * mu_300 * mu_102 * mu_003**2
        + 3 * mu_210**3 * mu_030
        + 2 * mu_210**2 * mu_201 * mu_003
        + 1 * mu_210**2 * mu_030**2
        + 1 * mu_210**2 * mu_030 * mu_012
        + 1 * mu_210**2 * mu_021 * mu_003
        + 2 * mu_210 * mu_201**2 * mu_030
        + 2 * mu_210 * mu_201 * mu_030 * mu_021
        + 2 * mu_210 * mu_201 * mu_012 * mu_003
        + 11 * mu_210 * mu_120**2 * mu_030
        + 4 * mu_210 * mu_120 * mu_111 * mu_003
        + 2 * mu_210 * mu_120 * mu_102 * mu_030
        + 8 * mu_210 * mu_111**2 * mu_030
        + 6 * mu_210 * mu_111 * mu_102 * mu_003
        + 1 * mu_210 * mu_102**2 * mu_030
        + 1 * mu_210 * mu_030**3
        + 4 * mu_210 * mu_030 * mu_021**2
        + 1 * mu_210 * mu_030 * mu_012**2
        + 2 * mu_210 * mu_021 * mu_012 * mu_003
        + 1 * mu_210 * mu_012 * mu_003**2
        + 3 * mu_201**3 * mu_003
        + 1 * mu_201**2 * mu_030 * mu_012
        + 1 * mu_201**2 * mu_021 * mu_003
        + 1 * mu_201**2 * mu_003**2
        + 1 * mu_201 * mu_120**2 * mu_003
        + 6 * mu_201 * mu_120 * mu_111 * mu_030
        + 2 * mu_201 * mu_120 * mu_102 * mu_003
        + 8 * mu_201 * mu_111**2 * mu_003
        + 4 * mu_201 * mu_111 * mu_102 * mu_030
        + 11 * mu_201 * mu_102**2 * mu_003
        + 1 * mu_201 * mu_030**2 * mu_021
        + 2 * mu_201 * mu_030 * mu_021 * mu_012
        + 1 * mu_201 * mu_021**2 * mu_003
        + 4 * mu_201 * mu_012**2 * mu_003
        + 1 * mu_201 * mu_003**3
        + 5 * mu_120**2 * mu_030**2
        + 4 * mu_120**2 * mu_030 * mu_012
        + 2 * mu_120**2 * mu_021 * mu_003
        + 14 * mu_120 * mu_111 * mu_030 * mu_021
        + 2 * mu_120 * mu_111 * mu_030 * mu_003
        + 6 * mu_120 * mu_111 * mu_012 * mu_003
        + 1 * mu_120 * mu_102 * mu_030**2
        + 2 * mu_120 * mu_102 * mu_030 * mu_012
        + 2 * mu_120 * mu_102 * mu_021 * mu_003
        + 1 * mu_120 * mu_102 * mu_003**2
        + 2 * mu_111**2 * mu_030**2
        + 8 * mu_111**2 * mu_030 * mu_012
        + 8 * mu_111**2 * mu_021 * mu_003
        + 2 * mu_111**2 * mu_003**2
        + 6 * mu_111 * mu_102 * mu_030 * mu_021
        + 2 * mu_111 * mu_102 * mu_030 * mu_003
        + 14 * mu_111 * mu_102 * mu_012 * mu_003
        + 2 * mu_102**2 * mu_030 * mu_012
        + 4 * mu_102**2 * mu_021 * mu_003
        + 5 * mu_102**2 * mu_003**2
        + 1 * mu_030**4
        + 1 * mu_030**3 * mu_012
        + 5 * mu_030**2 * mu_021**2
        + 1 * mu_030**2 * mu_021 * mu_003
        + 1 * mu_030**2 * mu_012**2
        + 11 * mu_030 * mu_021**2 * mu_012
        + 4 * mu_030 * mu_021 * mu_012 * mu_003
        + 3 * mu_030 * mu_012**3
        + 1 * mu_030 * mu_012 * mu_003**2
        + 3 * mu_021**3 * mu_003
        + 1 * mu_021**2 * mu_003**2
        + 11 * mu_021 * mu_012**2 * mu_003
        + 1 * mu_021 * mu_003**3
        + 5 * mu_012**2 * mu_003**2
        + 1 * mu_003**4
        + 2 * mu_210**4
        + 2 * mu_210**3 * mu_012
        + 4 * mu_210**2 * mu_201**2
        + 5 * mu_210**2 * mu_201 * mu_021
        + 10 * mu_210**2 * mu_120**2
        + 5 * mu_210**2 * mu_120 * mu_102
        + 6 * mu_210**2 * mu_111**2
        + 1 * mu_210**2 * mu_102**2
        + 1 * mu_210**2 * mu_021**2
        + 5 * mu_210 * mu_201**2 * mu_012
        + 18 * mu_210 * mu_201 * mu_120 * mu_111
        + 18 * mu_210 * mu_201 * mu_111 * mu_102
        + 4 * mu_210 * mu_201 * mu_021 * mu_012
        + 5 * mu_210 * mu_120**2 * mu_012
        + 18 * mu_210 * mu_120 * mu_111 * mu_021
        + 4 * mu_210 * mu_120 * mu_102 * mu_012
        + 12 * mu_210 * mu_111**2 * mu_012
        + 12 * mu_210 * mu_111 * mu_102 * mu_021
        + 5 * mu_210 * mu_102**2 * mu_012
        + 5 * mu_210 * mu_021**2 * mu_012
        + 2 * mu_210 * mu_012**3
        + 2 * mu_201**4
        + 2 * mu_201**3 * mu_021
        + 1 * mu_201**2 * mu_120**2
        + 5 * mu_201**2 * mu_120 * mu_102
        + 6 * mu_201**2 * mu_111**2
        + 10 * mu_201**2 * mu_102**2
        + 1 * mu_201**2 * mu_012**2
        + 5 * mu_201 * mu_120**2 * mu_021
        + 12 * mu_201 * mu_120 * mu_111 * mu_012
        + 4 * mu_201 * mu_120 * mu_102 * mu_021
        + 12 * mu_201 * mu_111**2 * mu_021
        + 18 * mu_201 * mu_111 * mu_102 * mu_012
        + 5 * mu_201 * mu_102**2 * mu_021
        + 2 * mu_201 * mu_021**3
        + 5 * mu_201 * mu_021 * mu_012**2
        + 2 * mu_120**4
        + 2 * mu_120**3 * mu_102
        + 6 * mu_120**2 * mu_111**2
        + 4 * mu_120**2 * mu_021**2
        + 1 * mu_120**2 * mu_012**2
        + 12 * mu_120 * mu_111**2 * mu_102
        + 18 * mu_120 * mu_111 * mu_021 * mu_012
        + 2 * mu_120 * mu_102**3
        + 5 * mu_120 * mu_102 * mu_021**2
        + 5 * mu_120 * mu_102 * mu_012**2
        + 6 * mu_111**2 * mu_102**2
        + 6 * mu_111**2 * mu_021**2
        + 6 * mu_111**2 * mu_012**2
        + 18 * mu_111 * mu_102 * mu_021 * mu_012
        + 2 * mu_102**4
        + 1 * mu_102**2 * mu_021**2
        + 4 * mu_102**2 * mu_012**2
        + 2 * mu_021**4
        + 10 * mu_021**2 * mu_012**2
        + 2 * mu_012**4
    )


@nb.njit
def phi_8(
    mu_030,
    mu_021,
    mu_120,
    mu_003,
    mu_111,
    mu_201,
    mu_102,
    mu_210,
    mu_012,
    mu_300,
):
    return (
        1 * mu_300**4
        + 2 * mu_300**3 * mu_120
        + 2 * mu_300**3 * mu_102
        + 4 * mu_300**2 * mu_210**2
        + 2 * mu_300**2 * mu_210 * mu_030
        + 2 * mu_300**2 * mu_210 * mu_012
        + 4 * mu_300**2 * mu_201**2
        + 2 * mu_300**2 * mu_201 * mu_021
        + 2 * mu_300**2 * mu_201 * mu_003
        + 2 * mu_300**2 * mu_120**2
        + 2 * mu_300**2 * mu_120 * mu_102
        + 2 * mu_300**2 * mu_111**2
        + 2 * mu_300**2 * mu_102**2
        + 10 * mu_300 * mu_210**2 * mu_120
        + 6 * mu_300 * mu_210**2 * mu_102
        + 8 * mu_300 * mu_210 * mu_201 * mu_111
        + 8 * mu_300 * mu_210 * mu_120 * mu_030
        + 6 * mu_300 * mu_210 * mu_120 * mu_012
        + 8 * mu_300 * mu_210 * mu_111 * mu_021
        + 4 * mu_300 * mu_210 * mu_111 * mu_003
        + 2 * mu_300 * mu_210 * mu_102 * mu_030
        + 4 * mu_300 * mu_210 * mu_102 * mu_012
        + 6 * mu_300 * mu_201**2 * mu_120
        + 10 * mu_300 * mu_201**2 * mu_102
        + 4 * mu_300 * mu_201 * mu_120 * mu_021
        + 2 * mu_300 * mu_201 * mu_120 * mu_003
        + 4 * mu_300 * mu_201 * mu_111 * mu_030
        + 8 * mu_300 * mu_201 * mu_111 * mu_012
        + 6 * mu_300 * mu_201 * mu_102 * mu_021
        + 8 * mu_300 * mu_201 * mu_102 * mu_003
        + 2 * mu_300 * mu_120**3
        + 2 * mu_300 * mu_120**2 * mu_102
        + 4 * mu_300 * mu_120 * mu_111**2
        + 2 * mu_300 * mu_120 * mu_102**2
        + 2 * mu_300 * mu_120 * mu_030**2
        + 2 * mu_300 * mu_120 * mu_030 * mu_012
        + 2 * mu_300 * mu_120 * mu_021**2
        + 2 * mu_300 * mu_120 * mu_021 * mu_003
        + 4 * mu_300 * mu_111**2 * mu_102
        + 4 * mu_300 * mu_111 * mu_030 * mu_021
        + 8 * mu_300 * mu_111 * mu_021 * mu_012
        + 4 * mu_300 * mu_111 * mu_012 * mu_003
        + 2 * mu_300 * mu_102**3
        + 2 * mu_300 * mu_102 * mu_030 * mu_012
        + 2 * mu_300 * mu_102 * mu_021 * mu_003
        + 2 * mu_300 * mu_102 * mu_012**2
        + 2 * mu_300 * mu_102 * mu_003**2
        + 2 * mu_210**3 * mu_030
        + 2 * mu_210**2 * mu_201 * mu_003
        + 2 * mu_210**2 * mu_030**2
        + 2 * mu_210**2 * mu_030 * mu_012
        + 2 * mu_210 * mu_201**2 * mu_030
        + 4 * mu_210 * mu_201 * mu_030 * mu_021
        + 2 * mu_210 * mu_201 * mu_030 * mu_003
        + 4 * mu_210 * mu_201 * mu_012 * mu_003
        + 10 * mu_210 * mu_120**2 * mu_030
        + 8 * mu_210 * mu_120 * mu_111 * mu_003
        + 6 * mu_210 * mu_120 * mu_102 * mu_030
        + 4 * mu_210 * mu_111**2 * mu_030
        + 8 * mu_210 * mu_111 * mu_102 * mu_003
        + 2 * mu_210 * mu_030**3
        + 2 * mu_210 * mu_030**2 * mu_012
        + 6 * mu_210 * mu_030 * mu_021**2
        + 2 * mu_210 * mu_030 * mu_021 * mu_003
        + 2 * mu_210 * mu_030 * mu_012**2
        + 6 * mu_210 * mu_021 * mu_012 * mu_003
        + 2 * mu_210 * mu_012 * mu_003**2
        + 2 * mu_201**3 * mu_003
        + 2 * mu_201**2 * mu_021 * mu_003
        + 2 * mu_201**2 * mu_003**2
        + 8 * mu_201 * mu_120 * mu_111 * mu_030
        + 6 * mu_201 * mu_120 * mu_102 * mu_003
        + 4 * mu_201 * mu_111**2 * mu_003
        + 8 * mu_201 * mu_111 * mu_102 * mu_030
        + 10 * mu_201 * mu_102**2 * mu_003
        + 2 * mu_201 * mu_030**2 * mu_021
        + 6 * mu_201 * mu_030 * mu_021 * mu_012
        + 2 * mu_201 * mu_030 * mu_012 * mu_003
        + 2 * mu_201 * mu_021**2 * mu_003
        + 2 * mu_201 * mu_021 * mu_003**2
        + 6 * mu_201 * mu_012**2 * mu_003
        + 2 * mu_201 * mu_003**3
        + 4 * mu_120**2 * mu_030**2
        + 6 * mu_120**2 * mu_030 * mu_012
        + 2 * mu_120**2 * mu_021 * mu_003
        + 8 * mu_120 * mu_111 * mu_030 * mu_021
        + 4 * mu_120 * mu_111 * mu_030 * mu_003
        + 8 * mu_120 * mu_111 * mu_012 * mu_003
        + 2 * mu_120 * mu_102 * mu_030**2
        + 4 * mu_120 * mu_102 * mu_030 * mu_012
        + 4 * mu_120 * mu_102 * mu_021 * mu_003
        + 2 * mu_120 * mu_102 * mu_003**2
        + 2 * mu_111**2 * mu_030**2
        + 4 * mu_111**2 * mu_030 * mu_012
        + 4 * mu_111**2 * mu_021 * mu_003
        + 2 * mu_111**2 * mu_003**2
        + 8 * mu_111 * mu_102 * mu_030 * mu_021
        + 4 * mu_111 * mu_102 * mu_030 * mu_003
        + 8 * mu_111 * mu_102 * mu_012 * mu_003
        + 2 * mu_102**2 * mu_030 * mu_012
        + 6 * mu_102**2 * mu_021 * mu_003
        + 4 * mu_102**2 * mu_003**2
        + 1 * mu_030**4
        + 2 * mu_030**3 * mu_012
        + 4 * mu_030**2 * mu_021**2
        + 2 * mu_030**2 * mu_021 * mu_003
        + 2 * mu_030**2 * mu_012**2
        + 10 * mu_030 * mu_021**2 * mu_012
        + 8 * mu_030 * mu_021 * mu_012 * mu_003
        + 2 * mu_030 * mu_012**3
        + 2 * mu_030 * mu_012 * mu_003**2
        + 2 * mu_021**3 * mu_003
        + 2 * mu_021**2 * mu_003**2
        + 10 * mu_021 * mu_012**2 * mu_003
        + 2 * mu_021 * mu_003**3
        + 4 * mu_012**2 * mu_003**2
        + 1 * mu_003**4
        + 1 * mu_210**4
        + 2 * mu_210**3 * mu_012
        + 2 * mu_210**2 * mu_201**2
        + 2 * mu_210**2 * mu_201 * mu_021
        + 8 * mu_210**2 * mu_120**2
        + 8 * mu_210**2 * mu_120 * mu_102
        + 2 * mu_210**2 * mu_111**2
        + 2 * mu_210**2 * mu_102**2
        + 2 * mu_210**2 * mu_021**2
        + 2 * mu_210**2 * mu_012**2
        + 2 * mu_210 * mu_201**2 * mu_012
        + 12 * mu_210 * mu_201 * mu_120 * mu_111
        + 12 * mu_210 * mu_201 * mu_111 * mu_102
        + 6 * mu_210 * mu_201 * mu_021 * mu_012
        + 8 * mu_210 * mu_120**2 * mu_012
        + 12 * mu_210 * mu_120 * mu_111 * mu_021
        + 6 * mu_210 * mu_120 * mu_102 * mu_012
        + 4 * mu_210 * mu_111**2 * mu_012
        + 12 * mu_210 * mu_111 * mu_102 * mu_021
        + 2 * mu_210 * mu_102**2 * mu_012
        + 8 * mu_210 * mu_021**2 * mu_012
        + 2 * mu_210 * mu_012**3
        + 1 * mu_201**4
        + 2 * mu_201**3 * mu_021
        + 2 * mu_201**2 * mu_120**2
        + 8 * mu_201**2 * mu_120 * mu_102
        + 2 * mu_201**2 * mu_111**2
        + 8 * mu_201**2 * mu_102**2
        + 2 * mu_201**2 * mu_021**2
        + 2 * mu_201**2 * mu_012**2
        + 2 * mu_201 * mu_120**2 * mu_021
        + 12 * mu_201 * mu_120 * mu_111 * mu_012
        + 6 * mu_201 * mu_120 * mu_102 * mu_021
        + 4 * mu_201 * mu_111**2 * mu_021
        + 12 * mu_201 * mu_111 * mu_102 * mu_012
        + 8 * mu_201 * mu_102**2 * mu_021
        + 2 * mu_201 * mu_021**3
        + 8 * mu_201 * mu_021 * mu_012**2
        + 1 * mu_120**4
        + 2 * mu_120**3 * mu_102
        + 2 * mu_120**2 * mu_111**2
        + 2 * mu_120**2 * mu_102**2
        + 2 * mu_120**2 * mu_021**2
        + 2 * mu_120**2 * mu_012**2
        + 4 * mu_120 * mu_111**2 * mu_102
        + 12 * mu_120 * mu_111 * mu_021 * mu_012
        + 2 * mu_120 * mu_102**3
        + 2 * mu_120 * mu_102 * mu_021**2
        + 2 * mu_120 * mu_102 * mu_012**2
        + 2 * mu_111**2 * mu_102**2
        + 2 * mu_111**2 * mu_021**2
        + 2 * mu_111**2 * mu_012**2
        + 12 * mu_111 * mu_102 * mu_021 * mu_012
        + 1 * mu_102**4
        + 2 * mu_102**2 * mu_021**2
        + 2 * mu_102**2 * mu_012**2
        + 1 * mu_021**4
        + 8 * mu_021**2 * mu_012**2
        + 1 * mu_012**4
    )


@nb.njit
def phi_9(
    mu_030,
    mu_021,
    mu_120,
    mu_101,
    mu_003,
    mu_200,
    mu_110,
    mu_201,
    mu_111,
    mu_102,
    mu_210,
    mu_020,
    mu_012,
    mu_002,
    mu_011,
    mu_300,
):
    return (
        1 * mu_200 * mu_300**2
        + 2 * mu_110 * mu_300 * mu_210
        + 2 * mu_110 * mu_120 * mu_030
        + 2 * mu_101 * mu_300 * mu_201
        + 2 * mu_101 * mu_102 * mu_003
        + 1 * mu_020 * mu_030**2
        + 2 * mu_011 * mu_030 * mu_021
        + 2 * mu_011 * mu_012 * mu_003
        + 1 * mu_002 * mu_003**2
        + 2 * mu_200 * mu_210**2
        + 2 * mu_200 * mu_201**2
        + 1 * mu_200 * mu_120**2
        + 2 * mu_200 * mu_111**2
        + 1 * mu_200 * mu_102**2
        + 4 * mu_110 * mu_210 * mu_120
        + 4 * mu_110 * mu_201 * mu_111
        + 4 * mu_110 * mu_111 * mu_021
        + 2 * mu_110 * mu_102 * mu_012
        + 4 * mu_101 * mu_210 * mu_111
        + 4 * mu_101 * mu_201 * mu_102
        + 2 * mu_101 * mu_120 * mu_021
        + 4 * mu_101 * mu_111 * mu_012
        + 1 * mu_020 * mu_210**2
        + 2 * mu_020 * mu_120**2
        + 2 * mu_020 * mu_111**2
        + 2 * mu_020 * mu_021**2
        + 1 * mu_020 * mu_012**2
        + 2 * mu_011 * mu_210 * mu_201
        + 4 * mu_011 * mu_120 * mu_111
        + 4 * mu_011 * mu_111 * mu_102
        + 4 * mu_011 * mu_021 * mu_012
        + 1 * mu_002 * mu_201**2
        + 2 * mu_002 * mu_111**2
        + 2 * mu_002 * mu_102**2
        + 1 * mu_002 * mu_021**2
        + 2 * mu_002 * mu_012**2
    )


@nb.njit
def phi_10(
    mu_030,
    mu_021,
    mu_120,
    mu_101,
    mu_003,
    mu_200,
    mu_110,
    mu_201,
    mu_111,
    mu_102,
    mu_210,
    mu_020,
    mu_012,
    mu_002,
    mu_011,
    mu_300,
):
    return (
        1 * mu_200 * mu_300**2
        + 1 * mu_200 * mu_300 * mu_120
        + 1 * mu_200 * mu_300 * mu_102
        + 1 * mu_200 * mu_210 * mu_030
        + 1 * mu_200 * mu_201 * mu_003
        + 2 * mu_110 * mu_300 * mu_210
        + 2 * mu_110 * mu_120 * mu_030
        + 2 * mu_110 * mu_111 * mu_003
        + 2 * mu_101 * mu_300 * mu_201
        + 2 * mu_101 * mu_111 * mu_030
        + 2 * mu_101 * mu_102 * mu_003
        + 1 * mu_020 * mu_300 * mu_120
        + 1 * mu_020 * mu_210 * mu_030
        + 1 * mu_020 * mu_030**2
        + 1 * mu_020 * mu_030 * mu_012
        + 1 * mu_020 * mu_021 * mu_003
        + 2 * mu_011 * mu_300 * mu_111
        + 2 * mu_011 * mu_030 * mu_021
        + 2 * mu_011 * mu_012 * mu_003
        + 1 * mu_002 * mu_300 * mu_102
        + 1 * mu_002 * mu_201 * mu_003
        + 1 * mu_002 * mu_030 * mu_012
        + 1 * mu_002 * mu_021 * mu_003
        + 1 * mu_002 * mu_003**2
        + 1 * mu_200 * mu_210**2
        + 1 * mu_200 * mu_210 * mu_012
        + 1 * mu_200 * mu_201**2
        + 1 * mu_200 * mu_201 * mu_021
        + 4 * mu_110 * mu_210 * mu_120
        + 2 * mu_110 * mu_210 * mu_102
        + 2 * mu_110 * mu_201 * mu_111
        + 2 * mu_110 * mu_120 * mu_012
        + 2 * mu_110 * mu_111 * mu_021
        + 2 * mu_101 * mu_210 * mu_111
        + 2 * mu_101 * mu_201 * mu_120
        + 4 * mu_101 * mu_201 * mu_102
        + 2 * mu_101 * mu_111 * mu_012
        + 2 * mu_101 * mu_102 * mu_021
        + 1 * mu_020 * mu_201 * mu_021
        + 1 * mu_020 * mu_120**2
        + 1 * mu_020 * mu_120 * mu_102
        + 1 * mu_020 * mu_021**2
        + 2 * mu_011 * mu_210 * mu_021
        + 2 * mu_011 * mu_201 * mu_012
        + 2 * mu_011 * mu_120 * mu_111
        + 2 * mu_011 * mu_111 * mu_102
        + 4 * mu_011 * mu_021 * mu_012
        + 1 * mu_002 * mu_210 * mu_012
        + 1 * mu_002 * mu_120 * mu_102
        + 1 * mu_002 * mu_102**2
        + 1 * mu_002 * mu_012**2
    )


@nb.njit
def phi_11(
    mu_030,
    mu_021,
    mu_120,
    mu_101,
    mu_003,
    mu_200,
    mu_110,
    mu_201,
    mu_102,
    mu_210,
    mu_012,
    mu_020,
    mu_002,
    mu_011,
    mu_300,
):
    return (
        1 * mu_200 * mu_300**2
        + 2 * mu_200 * mu_300 * mu_120
        + 2 * mu_200 * mu_300 * mu_102
        + 2 * mu_110 * mu_300 * mu_210
        + 2 * mu_110 * mu_300 * mu_030
        + 2 * mu_110 * mu_300 * mu_012
        + 2 * mu_110 * mu_120 * mu_030
        + 2 * mu_110 * mu_102 * mu_030
        + 2 * mu_101 * mu_300 * mu_201
        + 2 * mu_101 * mu_300 * mu_021
        + 2 * mu_101 * mu_300 * mu_003
        + 2 * mu_101 * mu_120 * mu_003
        + 2 * mu_101 * mu_102 * mu_003
        + 2 * mu_020 * mu_210 * mu_030
        + 1 * mu_020 * mu_030**2
        + 2 * mu_020 * mu_030 * mu_012
        + 2 * mu_011 * mu_210 * mu_003
        + 2 * mu_011 * mu_201 * mu_030
        + 2 * mu_011 * mu_030 * mu_021
        + 2 * mu_011 * mu_030 * mu_003
        + 2 * mu_011 * mu_012 * mu_003
        + 2 * mu_002 * mu_201 * mu_003
        + 2 * mu_002 * mu_021 * mu_003
        + 1 * mu_002 * mu_003**2
        + 1 * mu_200 * mu_120**2
        + 2 * mu_200 * mu_120 * mu_102
        + 1 * mu_200 * mu_102**2
        + 2 * mu_110 * mu_210 * mu_120
        + 2 * mu_110 * mu_210 * mu_102
        + 2 * mu_110 * mu_120 * mu_012
        + 2 * mu_110 * mu_102 * mu_012
        + 2 * mu_101 * mu_201 * mu_120
        + 2 * mu_101 * mu_201 * mu_102
        + 2 * mu_101 * mu_120 * mu_021
        + 2 * mu_101 * mu_102 * mu_021
        + 1 * mu_020 * mu_210**2
        + 2 * mu_020 * mu_210 * mu_012
        + 1 * mu_020 * mu_012**2
        + 2 * mu_011 * mu_210 * mu_201
        + 2 * mu_011 * mu_210 * mu_021
        + 2 * mu_011 * mu_201 * mu_012
        + 2 * mu_011 * mu_021 * mu_012
        + 1 * mu_002 * mu_201**2
        + 2 * mu_002 * mu_201 * mu_021
        + 1 * mu_002 * mu_021**2
    )


@nb.njit
def phi_12(
    mu_030,
    mu_021,
    mu_120,
    mu_101,
    mu_003,
    mu_200,
    mu_110,
    mu_201,
    mu_111,
    mu_102,
    mu_210,
    mu_020,
    mu_012,
    mu_002,
    mu_011,
    mu_300,
):
    return (
        1 * mu_200**2 * mu_300**2
        + 4 * mu_200 * mu_110 * mu_300 * mu_210
        + 4 * mu_200 * mu_101 * mu_300 * mu_201
        + 2 * mu_200 * mu_020 * mu_300 * mu_120
        + 2 * mu_200 * mu_020 * mu_210 * mu_030
        + 4 * mu_200 * mu_011 * mu_300 * mu_111
        + 2 * mu_200 * mu_002 * mu_300 * mu_102
        + 2 * mu_200 * mu_002 * mu_201 * mu_003
        + 4 * mu_110 * mu_020 * mu_120 * mu_030
        + 4 * mu_110 * mu_002 * mu_111 * mu_003
        + 4 * mu_101 * mu_020 * mu_111 * mu_030
        + 4 * mu_101 * mu_002 * mu_102 * mu_003
        + 1 * mu_020**2 * mu_030**2
        + 4 * mu_020 * mu_011 * mu_030 * mu_021
        + 2 * mu_020 * mu_002 * mu_030 * mu_012
        + 2 * mu_020 * mu_002 * mu_021 * mu_003
        + 4 * mu_011 * mu_002 * mu_012 * mu_003
        + 1 * mu_002**2 * mu_003**2
        + 1 * mu_200**2 * mu_210**2
        + 1 * mu_200**2 * mu_201**2
        + 4 * mu_200 * mu_110 * mu_210 * mu_120
        + 4 * mu_200 * mu_110 * mu_201 * mu_111
        + 4 * mu_200 * mu_101 * mu_210 * mu_111
        + 4 * mu_200 * mu_101 * mu_201 * mu_102
        + 2 * mu_200 * mu_020 * mu_201 * mu_021
        + 4 * mu_200 * mu_011 * mu_210 * mu_021
        + 4 * mu_200 * mu_011 * mu_201 * mu_012
        + 2 * mu_200 * mu_002 * mu_210 * mu_012
        + 4 * mu_110**2 * mu_210**2
        + 4 * mu_110**2 * mu_120**2
        + 8 * mu_110 * mu_101 * mu_210 * mu_201
        + 8 * mu_110 * mu_101 * mu_120 * mu_111
        + 8 * mu_110 * mu_101 * mu_111 * mu_102
        + 4 * mu_110 * mu_020 * mu_210 * mu_120
        + 4 * mu_110 * mu_020 * mu_111 * mu_021
        + 8 * mu_110 * mu_011 * mu_210 * mu_111
        + 8 * mu_110 * mu_011 * mu_120 * mu_021
        + 8 * mu_110 * mu_011 * mu_111 * mu_012
        + 4 * mu_110 * mu_002 * mu_210 * mu_102
        + 4 * mu_110 * mu_002 * mu_120 * mu_012
        + 4 * mu_101**2 * mu_201**2
        + 4 * mu_101**2 * mu_102**2
        + 4 * mu_101 * mu_020 * mu_201 * mu_120
        + 4 * mu_101 * mu_020 * mu_102 * mu_021
        + 8 * mu_101 * mu_011 * mu_201 * mu_111
        + 8 * mu_101 * mu_011 * mu_111 * mu_021
        + 8 * mu_101 * mu_011 * mu_102 * mu_012
        + 4 * mu_101 * mu_002 * mu_201 * mu_102
        + 4 * mu_101 * mu_002 * mu_111 * mu_012
        + 1 * mu_020**2 * mu_120**2
        + 1 * mu_020**2 * mu_021**2
        + 4 * mu_020 * mu_011 * mu_120 * mu_111
        + 4 * mu_020 * mu_011 * mu_021 * mu_012
        + 2 * mu_020 * mu_002 * mu_120 * mu_102
        + 4 * mu_011**2 * mu_021**2
        + 4 * mu_011**2 * mu_012**2
        + 4 * mu_011 * mu_002 * mu_111 * mu_102
        + 4 * mu_011 * mu_002 * mu_021 * mu_012
        + 1 * mu_002**2 * mu_102**2
        + 1 * mu_002**2 * mu_012**2
        + 4 * mu_110**2 * mu_111**2
        + 4 * mu_101**2 * mu_111**2
        + 4 * mu_011**2 * mu_111**2
    )


@nb.njit
def phi_13(
    mu_030,
    mu_021,
    mu_120,
    mu_101,
    mu_003,
    mu_200,
    mu_110,
    mu_201,
    mu_111,
    mu_102,
    mu_210,
    mu_012,
    mu_020,
    mu_002,
    mu_011,
    mu_300,
):
    return (
        1 * mu_200**2 * mu_300**2
        + 2 * mu_200 * mu_110 * mu_300 * mu_210
        + 2 * mu_200 * mu_110 * mu_120 * mu_030
        + 2 * mu_200 * mu_101 * mu_300 * mu_201
        + 2 * mu_200 * mu_101 * mu_102 * mu_003
        + 1 * mu_110**2 * mu_300**2
        + 1 * mu_110**2 * mu_030**2
        + 2 * mu_110 * mu_101 * mu_030 * mu_021
        + 2 * mu_110 * mu_101 * mu_012 * mu_003
        + 2 * mu_110 * mu_020 * mu_300 * mu_210
        + 2 * mu_110 * mu_020 * mu_120 * mu_030
        + 2 * mu_110 * mu_011 * mu_300 * mu_201
        + 2 * mu_110 * mu_011 * mu_102 * mu_003
        + 1 * mu_101**2 * mu_300**2
        + 1 * mu_101**2 * mu_003**2
        + 2 * mu_101 * mu_011 * mu_300 * mu_210
        + 2 * mu_101 * mu_011 * mu_120 * mu_030
        + 2 * mu_101 * mu_002 * mu_300 * mu_201
        + 2 * mu_101 * mu_002 * mu_102 * mu_003
        + 1 * mu_020**2 * mu_030**2
        + 2 * mu_020 * mu_011 * mu_030 * mu_021
        + 2 * mu_020 * mu_011 * mu_012 * mu_003
        + 1 * mu_011**2 * mu_030**2
        + 1 * mu_011**2 * mu_003**2
        + 2 * mu_011 * mu_002 * mu_030 * mu_021
        + 2 * mu_011 * mu_002 * mu_012 * mu_003
        + 1 * mu_002**2 * mu_003**2
        + 2 * mu_200**2 * mu_210**2
        + 2 * mu_200**2 * mu_201**2
        + 1 * mu_200**2 * mu_120**2
        + 2 * mu_200**2 * mu_111**2
        + 1 * mu_200**2 * mu_102**2
        + 4 * mu_200 * mu_110 * mu_210 * mu_120
        + 4 * mu_200 * mu_110 * mu_201 * mu_111
        + 4 * mu_200 * mu_110 * mu_111 * mu_021
        + 2 * mu_200 * mu_110 * mu_102 * mu_012
        + 4 * mu_200 * mu_101 * mu_210 * mu_111
        + 4 * mu_200 * mu_101 * mu_201 * mu_102
        + 2 * mu_200 * mu_101 * mu_120 * mu_021
        + 4 * mu_200 * mu_101 * mu_111 * mu_012
        + 3 * mu_110**2 * mu_210**2
        + 2 * mu_110**2 * mu_201**2
        + 3 * mu_110**2 * mu_120**2
        + 1 * mu_110**2 * mu_102**2
        + 2 * mu_110**2 * mu_021**2
        + 1 * mu_110**2 * mu_012**2
        + 2 * mu_110 * mu_101 * mu_210 * mu_201
        + 4 * mu_110 * mu_101 * mu_120 * mu_111
        + 4 * mu_110 * mu_101 * mu_111 * mu_102
        + 4 * mu_110 * mu_101 * mu_021 * mu_012
        + 4 * mu_110 * mu_020 * mu_210 * mu_120
        + 4 * mu_110 * mu_020 * mu_201 * mu_111
        + 4 * mu_110 * mu_020 * mu_111 * mu_021
        + 2 * mu_110 * mu_020 * mu_102 * mu_012
        + 4 * mu_110 * mu_011 * mu_210 * mu_111
        + 4 * mu_110 * mu_011 * mu_201 * mu_102
        + 2 * mu_110 * mu_011 * mu_120 * mu_021
        + 4 * mu_110 * mu_011 * mu_111 * mu_012
        + 2 * mu_101**2 * mu_210**2
        + 3 * mu_101**2 * mu_201**2
        + 1 * mu_101**2 * mu_120**2
        + 3 * mu_101**2 * mu_102**2
        + 1 * mu_101**2 * mu_021**2
        + 2 * mu_101**2 * mu_012**2
        + 4 * mu_101 * mu_011 * mu_210 * mu_120
        + 4 * mu_101 * mu_011 * mu_201 * mu_111
        + 4 * mu_101 * mu_011 * mu_111 * mu_021
        + 2 * mu_101 * mu_011 * mu_102 * mu_012
        + 4 * mu_101 * mu_002 * mu_210 * mu_111
        + 4 * mu_101 * mu_002 * mu_201 * mu_102
        + 2 * mu_101 * mu_002 * mu_120 * mu_021
        + 4 * mu_101 * mu_002 * mu_111 * mu_012
        + 1 * mu_020**2 * mu_210**2
        + 2 * mu_020**2 * mu_120**2
        + 2 * mu_020**2 * mu_111**2
        + 2 * mu_020**2 * mu_021**2
        + 1 * mu_020**2 * mu_012**2
        + 2 * mu_020 * mu_011 * mu_210 * mu_201
        + 4 * mu_020 * mu_011 * mu_120 * mu_111
        + 4 * mu_020 * mu_011 * mu_111 * mu_102
        + 4 * mu_020 * mu_011 * mu_021 * mu_012
        + 1 * mu_011**2 * mu_210**2
        + 1 * mu_011**2 * mu_201**2
        + 2 * mu_011**2 * mu_120**2
        + 2 * mu_011**2 * mu_102**2
        + 3 * mu_011**2 * mu_021**2
        + 3 * mu_011**2 * mu_012**2
        + 2 * mu_011 * mu_002 * mu_210 * mu_201
        + 4 * mu_011 * mu_002 * mu_120 * mu_111
        + 4 * mu_011 * mu_002 * mu_111 * mu_102
        + 4 * mu_011 * mu_002 * mu_021 * mu_012
        + 1 * mu_002**2 * mu_201**2
        + 2 * mu_002**2 * mu_111**2
        + 2 * mu_002**2 * mu_102**2
        + 1 * mu_002**2 * mu_021**2
        + 2 * mu_002**2 * mu_012**2
        + 4 * mu_110**2 * mu_111**2
        + 4 * mu_101**2 * mu_111**2
        + 4 * mu_011**2 * mu_111**2
    )


@nb.njit
def CI(
    mu_000,
    mu_200,
    mu_020,
    mu_002,
    mu_110,
    mu_101,
    mu_011,
    mu_111,
    mu_210,
    mu_201,
    mu_120,
    mu_021,
    mu_012,
    mu_102,
    mu_003,
    mu_030,
    mu_300,
    mu_013,
    mu_103,
    mu_130,
    mu_310,
    mu_301,
    mu_031,
    mu_112,
    mu_121,
    mu_211,
    mu_022,
    mu_202,
    mu_220,
    mu_400,
    mu_040,
    mu_004,
):
    alpha_1 = mu_002 - mu_020
    alpha_2 = mu_020 - mu_200
    alpha_3 = mu_200 - mu_002
    beta_1 = mu_021 - mu_201
    beta_2 = mu_102 - mu_120
    beta_3 = mu_210 - mu_012
    beta_4 = mu_003 - mu_201 - 2 * mu_021
    beta_5 = mu_030 - mu_012 - 2 * mu_210
    beta_6 = mu_300 - mu_120 - 2 * mu_102
    beta_7 = mu_021 - mu_003 + 2 * mu_201
    beta_8 = mu_102 - mu_300 + 2 * mu_120
    beta_9 = mu_210 - mu_030 + 2 * mu_012
    beta_10 = mu_021 + mu_201 - 3 * mu_003
    beta_11 = mu_012 + mu_210 - 3 * mu_030
    beta_12 = mu_102 + mu_120 - 3 * mu_300
    beta_13 = mu_021 + mu_003 + 3 * mu_201
    beta_14 = mu_102 + mu_300 + 3 * mu_120
    beta_15 = mu_210 + mu_030 + 3 * mu_012
    beta_16 = mu_012 + mu_030 + 3 * mu_210
    beta_17 = mu_201 + mu_003 + 3 * mu_021
    beta_18 = mu_120 + mu_300 + 3 * mu_102
    gamma_1 = mu_022 - mu_400
    gamma_2 = mu_202 - mu_040
    gamma_3 = mu_220 - mu_004
    gamma_4 = mu_112 + mu_130 + mu_310
    gamma_5 = mu_121 + mu_103 + mu_301
    gamma_6 = mu_211 + mu_013 + mu_031
    gamma_7 = mu_022 - mu_220 + mu_004 - mu_400
    gamma_8 = mu_202 - mu_022 + mu_400 - mu_040
    gamma_9 = mu_220 - mu_202 + mu_040 - mu_004
    r_gyr = np.sqrt((mu_200 + mu_020 + mu_002) / (3 * mu_000))
    if r_gyr == 0:
        return np.nan
    s_3 = 1 / (mu_000**3 * r_gyr**9)
    s_4 = 1 / (mu_000**4 * r_gyr**9)
    return 4 * s_3 * (
        mu_110
        * (
            mu_021 * (3 * gamma_2 - 2 * gamma_3 - gamma_1)
            - mu_201 * (3 * gamma_1 - 2 * gamma_3 - gamma_2)
            + beta_12 * gamma_5
            - beta_11 * gamma_6
            + mu_003 * gamma_8
        )
        + mu_101
        * (
            mu_210 * (3 * gamma_1 - 2 * gamma_2 - gamma_3)
            - mu_012 * (3 * gamma_3 - 2 * gamma_2 - gamma_1)
            + beta_10 * gamma_6
            - beta_12 * gamma_4
            + mu_030 * gamma_7
        )
        + mu_011
        * (
            mu_102 * (3 * gamma_3 - 2 * gamma_1 - gamma_2)
            - mu_120 * (3 * gamma_2 - 2 * gamma_1 - gamma_3)
            + beta_11 * gamma_4
            - beta_10 * gamma_5
            + mu_300 * gamma_9
        )
        + mu_002 * (beta_18 * gamma_6 - beta_15 * gamma_5 - 2 * (mu_111 * gamma_8 + beta_1 * gamma_4))
        + mu_020 * (beta_17 * gamma_4 - beta_14 * gamma_6 - 2 * (mu_111 * gamma_7 + beta_3 * gamma_5))
        + mu_200 * (beta_16 * gamma_5 - beta_13 * gamma_4 - 2 * (mu_111 * gamma_9 + beta_2 * gamma_6))
    ) - 16 * s_4 * (
        mu_011 * alpha_2 * alpha_3 * beta_2
        + mu_101 * alpha_1 * alpha_2 * beta_3
        + mu_110 * alpha_1 * alpha_3 * beta_1
        - mu_111 * alpha_1 * alpha_2 * alpha_3
        - mu_011 * mu_011 * (mu_111 * alpha_1 - mu_011 * beta_2 - mu_101 * beta_5 - mu_110 * beta_7)
        - mu_101 * mu_101 * (mu_111 * alpha_3 - mu_101 * beta_3 - mu_110 * beta_4 - mu_011 * beta_8)
        - mu_110 * mu_110 * (mu_111 * alpha_2 - mu_110 * beta_1 - mu_011 * beta_6 - mu_101 * beta_9)
        + mu_011 * mu_101 * (mu_002 * beta_1 + mu_020 * beta_4 + mu_200 * beta_7)
        + mu_011 * mu_110 * (mu_020 * beta_3 + mu_200 * beta_5 + mu_002 * beta_9)
        + mu_101 * mu_110 * (mu_200 * beta_2 + mu_002 * beta_6 + mu_020 * beta_8)
    )


class MomentType(Enum):
    """
    Different rotation invariant moments (order 2 and order 3)

    Choose from ['O_3', 'O_4', 'O_5', 'F', 'phi_2', 'phi_3', 'phi_4', 'phi_5', 'phi_6', 'phi_7', 'phi_8', 'phi_9', 'phi_10', 'phi_11', 'phi_12', 'phi_13']

    O_3, O_4, and O_5 are second order moments from [1] and F is a third order moment from [2].
    These four moments are used in the original Geometricus manuscript [3].

    phi_{2-13} are independent third order moments from [4].

    CI is the chiral invariant moment from [5].


    [1] Mamistvalov, Alexander G. "N-dimensional moment invariants and conceptual mathematical theory of recognition n-dimensional solids."
    IEEE Transactions on pattern analysis and machine intelligence 20.8 (1998): 819-831.

    [2] Flusser, Jan, Jirí Boldys, and Barbara Zitová. "Moment forms invariant to rotation and blur in arbitrary number of dimensions."
    IEEE Transactions on Pattern Analysis and Machine Intelligence 25.2 (2003): 234-246.

    [3] Durairaj, Janani, et al. "Geometricus represents protein structures as shape-mers derived from moment invariants." Bioinformatics 36.Supplement_2 (2020): i718-i725.

    [4] Flusser, Jan, Tomas Suk, and Barbara Zitová. 2D and 3D image analysis by moments. John Wiley & Sons, 2016.

    [5] Hattne, Johan, and Victor S. Lamzin. "A moment invariant for evaluating the chirality of three-dimensional objects." Journal of The Royal Society Interface 8.54 (2011): 144-151.
    """

    O_3 = MomentInfo(O_3, [(2, 0, 0), (0, 2, 0), (0, 0, 2)])
    O_4 = MomentInfo(O_4, [(2, 0, 0), (0, 2, 0), (0, 0, 2), (1, 1, 0), (1, 0, 1), (0, 1, 1)])
    O_5 = MomentInfo(O_5, [(2, 0, 0), (0, 2, 0), (0, 0, 2), (1, 1, 0), (1, 0, 1), (0, 1, 1)])
    F = MomentInfo(
        F,
        [
            (2, 0, 1),
            (0, 2, 1),
            (2, 1, 0),
            (3, 0, 0),
            (1, 1, 1),
            (0, 1, 2),
            (0, 0, 3),
            (0, 3, 0),
            (1, 0, 2),
            (1, 2, 0),
        ],
    )
    phi_2 = MomentInfo(phi_2, [(0, 2, 0), (0, 1, 1), (1, 1, 0), (2, 0, 0), (0, 0, 2), (1, 0, 1)])
    phi_3 = MomentInfo(phi_3, [(0, 2, 0), (0, 1, 1), (1, 1, 0), (2, 0, 0), (0, 0, 2), (1, 0, 1)])
    phi_4 = MomentInfo(
        phi_4,
        [
            (0, 3, 0),
            (0, 2, 1),
            (1, 2, 0),
            (0, 0, 3),
            (1, 1, 1),
            (2, 0, 1),
            (1, 0, 2),
            (2, 1, 0),
            (0, 1, 2),
            (3, 0, 0),
        ],
    )
    phi_5 = MomentInfo(
        phi_5,
        [
            (0, 3, 0),
            (0, 2, 1),
            (1, 2, 0),
            (0, 0, 3),
            (2, 0, 1),
            (1, 0, 2),
            (2, 1, 0),
            (0, 1, 2),
            (3, 0, 0),
        ],
    )
    phi_6 = MomentInfo(
        phi_6,
        [
            (0, 3, 0),
            (0, 2, 1),
            (1, 2, 0),
            (0, 0, 3),
            (1, 1, 1),
            (2, 0, 1),
            (1, 0, 2),
            (2, 1, 0),
            (0, 1, 2),
            (3, 0, 0),
        ],
    )
    phi_7 = MomentInfo(
        phi_7,
        [
            (0, 3, 0),
            (0, 2, 1),
            (1, 2, 0),
            (0, 0, 3),
            (1, 1, 1),
            (2, 0, 1),
            (1, 0, 2),
            (2, 1, 0),
            (0, 1, 2),
            (3, 0, 0),
        ],
    )
    phi_8 = MomentInfo(
        phi_8,
        [
            (0, 3, 0),
            (0, 2, 1),
            (1, 2, 0),
            (0, 0, 3),
            (1, 1, 1),
            (2, 0, 1),
            (1, 0, 2),
            (2, 1, 0),
            (0, 1, 2),
            (3, 0, 0),
        ],
    )
    phi_9 = MomentInfo(
        phi_9,
        [
            (0, 3, 0),
            (0, 2, 1),
            (1, 2, 0),
            (1, 0, 1),
            (0, 0, 3),
            (2, 0, 0),
            (1, 1, 0),
            (2, 0, 1),
            (1, 1, 1),
            (1, 0, 2),
            (2, 1, 0),
            (0, 2, 0),
            (0, 1, 2),
            (0, 0, 2),
            (0, 1, 1),
            (3, 0, 0),
        ],
    )
    phi_10 = MomentInfo(
        phi_10,
        [
            (0, 3, 0),
            (0, 2, 1),
            (1, 2, 0),
            (1, 0, 1),
            (0, 0, 3),
            (2, 0, 0),
            (1, 1, 0),
            (2, 0, 1),
            (1, 1, 1),
            (1, 0, 2),
            (2, 1, 0),
            (0, 2, 0),
            (0, 1, 2),
            (0, 0, 2),
            (0, 1, 1),
            (3, 0, 0),
        ],
    )
    phi_11 = MomentInfo(
        phi_11,
        [
            (0, 3, 0),
            (0, 2, 1),
            (1, 2, 0),
            (1, 0, 1),
            (0, 0, 3),
            (2, 0, 0),
            (1, 1, 0),
            (2, 0, 1),
            (1, 0, 2),
            (2, 1, 0),
            (0, 1, 2),
            (0, 2, 0),
            (0, 0, 2),
            (0, 1, 1),
            (3, 0, 0),
        ],
    )
    phi_12 = MomentInfo(
        phi_12,
        [
            (0, 3, 0),
            (0, 2, 1),
            (1, 2, 0),
            (1, 0, 1),
            (0, 0, 3),
            (2, 0, 0),
            (1, 1, 0),
            (2, 0, 1),
            (1, 1, 1),
            (1, 0, 2),
            (2, 1, 0),
            (0, 2, 0),
            (0, 1, 2),
            (0, 0, 2),
            (0, 1, 1),
            (3, 0, 0),
        ],
    )
    phi_13 = MomentInfo(
        phi_13,
        [
            (0, 3, 0),
            (0, 2, 1),
            (1, 2, 0),
            (1, 0, 1),
            (0, 0, 3),
            (2, 0, 0),
            (1, 1, 0),
            (2, 0, 1),
            (1, 1, 1),
            (1, 0, 2),
            (2, 1, 0),
            (0, 1, 2),
            (0, 2, 0),
            (0, 0, 2),
            (0, 1, 1),
            (3, 0, 0),
        ],
    )
    CI = MomentInfo(
        CI,
        [
            (0, 0, 0),
            (2, 0, 0),
            (0, 2, 0),
            (0, 0, 2),
            (1, 1, 0),
            (1, 0, 1),
            (0, 1, 1),
            (1, 1, 1),
            (2, 1, 0),
            (2, 0, 1),
            (1, 2, 0),
            (0, 2, 1),
            (0, 1, 2),
            (1, 0, 2),
            (0, 0, 3),
            (0, 3, 0),
            (3, 0, 0),
            (0, 1, 3),
            (1, 0, 3),
            (1, 3, 0),
            (3, 1, 0),
            (3, 0, 1),
            (0, 3, 1),
            (1, 1, 2),
            (1, 2, 1),
            (2, 1, 1),
            (0, 2, 2),
            (2, 0, 2),
            (2, 2, 0),
            (4, 0, 0),
            (0, 4, 0),
            (0, 0, 4),
        ],
    )

    def get_moments_from_coordinates(self, mus: ty.List[float]):
        return self.value.moment_function(*mus)


def get_moments_from_coordinates(
    coordinates: np.ndarray,
    moment_types: ty.List[MomentType] = (
        MomentType.O_3,
        MomentType.O_4,
        MomentType.O_5,
        MomentType.F,
    ),
) -> ty.List[float]:
    """
    Gets rotation-invariant moments for a set of coordinates

    Parameters
    ----------
    coordinates
    moment_types
        Which moments to calculate
        Choose from ['O_3', 'O_4', 'O_5', 'F', 'phi_2', 'phi_3', 'phi_4', 'phi_5', 'phi_6', 'phi_7', 'phi_8', 'phi_9', 'phi_10', 'phi_11', 'phi_12', 'phi_13', 'CI']

    Returns
    -------
    list of moments
    """
    if coordinates.shape[0] < 2:
        return [np.nan] * len(moment_types)
    all_moment_mu_types: ty.Set[ty.Tuple[int, int, int]] = set(
        m for moment_type in moment_types for m in moment_type.value.mu_arguments
    )
    centroid = nb_mean_axis_0(coordinates)
    mus = {(x, y, z): mu(float(x), float(y), float(z), coordinates, centroid) for (x, y, z) in all_moment_mu_types}
    moments = [
        moment_type.get_moments_from_coordinates([mus[m] for m in moment_type.value.mu_arguments])
        for moment_type in moment_types
    ]
    return moments


def normalize(moments):
    return (np.sign(moments) * np.log1p(np.abs(moments))).astype("float32")


def get_moments(coords):
    return normalize(
        get_moments_from_coordinates(
            np.array(coords),
            moment_types=[
                MomentType.O_3,
                MomentType.O_4,
                MomentType.O_5,
                MomentType.F,
                MomentType.phi_2,
                MomentType.phi_3,
                MomentType.phi_4,
                MomentType.phi_5,
                MomentType.phi_6,
                MomentType.phi_7,
                MomentType.phi_8,
                MomentType.phi_9,
                MomentType.phi_10,
                MomentType.phi_11,
                MomentType.phi_12,
                MomentType.phi_13,
                MomentType.CI,
            ],
        )
    )
# -*- coding: utf-8 -*-
#
"""
Stefanos-Aldo Papanicolopulos,
New fully symmetric and rotationally symmetric cubature rules on the triangle using
minimal orthonormal bases,
Journal of Computational and Applied Mathematics,
Volume 294, 1 March 2016, Pages 39–48,
<https://doi.org/10.1016/j.cam.2015.08.001>,
<https://arxiv.org/abs/1411.5631>.
"""
import json
import os

from ..helpers import untangle2, TriangleScheme


def _read(filename):
    this_dir = os.path.dirname(os.path.realpath(__file__))
    with open(os.path.join(this_dir, filename), "r") as f:
        data = json.load(f)

    degree = data.pop("degree")
    bary, weights = untangle2(data)
    return degree, weights, bary


def papanicolopulos_sym_0():
    return TriangleScheme("Papanicolopulos 0 (full symmetry)", *_read("full00.json"))


def papanicolopulos_sym_1():
    return TriangleScheme("Papanicolopulos 1 (full symmetry)", *_read("full01.json"))


def papanicolopulos_sym_2():
    return TriangleScheme("Papanicolopulos 2 (full symmetry)", *_read("full02.json"))


def papanicolopulos_sym_3():
    return TriangleScheme("Papanicolopulos 3 (full symmetry)", *_read("full03.json"))


def papanicolopulos_sym_4():
    return TriangleScheme("Papanicolopulos 4 (full symmetry)", *_read("full04.json"))


def papanicolopulos_sym_5():
    return TriangleScheme("Papanicolopulos 5 (full symmetry)", *_read("full05.json"))


def papanicolopulos_sym_6():
    return TriangleScheme("Papanicolopulos 6 (full symmetry)", *_read("full06.json"))


def papanicolopulos_sym_7():
    return TriangleScheme("Papanicolopulos 7 (full symmetry)", *_read("full07.json"))


def papanicolopulos_sym_8():
    return TriangleScheme("Papanicolopulos 8 (full symmetry)", *_read("full08.json"))


PapanicolopulosSym = {
    0: papanicolopulos_sym_0,
    1: papanicolopulos_sym_1,
    2: papanicolopulos_sym_2,
    3: papanicolopulos_sym_3,
    4: papanicolopulos_sym_4,
    5: papanicolopulos_sym_5,
    6: papanicolopulos_sym_6,
    7: papanicolopulos_sym_7,
    8: papanicolopulos_sym_8,
}

# TODO ERR the first 8 schemes are flawed by round-off error
# def papanicolopulos_rot_0():
#     return TriangleScheme("Papanicolopulos 0 (rotational symmetry)", *_read("rot00.json"))
#
#
# def papanicolopulos_rot_1():
#     return TriangleScheme("Papanicolopulos 1 (rotational symmetry)", *_read("rot01.json"))
#
#
# def papanicolopulos_rot_2():
#     return TriangleScheme("Papanicolopulos 2 (rotational symmetry)", *_read("rot02.json"))
#
#
# def papanicolopulos_rot_3():
#     return TriangleScheme("Papanicolopulos 3 (rotational symmetry)", *_read("rot03.json"))
#
#
# def papanicolopulos_rot_4():
#     return TriangleScheme("Papanicolopulos 4 (rotational symmetry)", *_read("rot04.json"))
#
#
# def papanicolopulos_rot_5():
#     return TriangleScheme("Papanicolopulos 5 (rotational symmetry)", *_read("rot05.json"))
#
#
# def papanicolopulos_rot_6():
#     return TriangleScheme("Papanicolopulos 6 (rotational symmetry)", *_read("rot06.json"))
#
#
# def papanicolopulos_rot_7():
#     return TriangleScheme("Papanicolopulos 7 (rotational symmetry)", *_read("rot07.json"))


def papanicolopulos_rot_8():
    return TriangleScheme(
        "Papanicolopulos 8 (rotational symmetry)", *_read("rot08.json")
    )


def papanicolopulos_rot_9():
    return TriangleScheme(
        "Papanicolopulos 9 (rotational symmetry)", *_read("rot09.json")
    )


def papanicolopulos_rot_10():
    return TriangleScheme(
        "Papanicolopulos 10 (rotational symmetry)", *_read("rot10.json")
    )


def papanicolopulos_rot_11():
    return TriangleScheme(
        "Papanicolopulos 11 (rotational symmetry)", *_read("rot11.json")
    )


def papanicolopulos_rot_12():
    return TriangleScheme(
        "Papanicolopulos 12 (rotational symmetry)", *_read("rot12.json")
    )


def papanicolopulos_rot_13():
    return TriangleScheme(
        "Papanicolopulos 13 (rotational symmetry)", *_read("rot13.json")
    )


def papanicolopulos_rot_14():
    return TriangleScheme(
        "Papanicolopulos 14 (rotational symmetry)", *_read("rot14.json")
    )


def papanicolopulos_rot_15():
    return TriangleScheme(
        "Papanicolopulos 15 (rotational symmetry)", *_read("rot15.json")
    )


def papanicolopulos_rot_16():
    return TriangleScheme(
        "Papanicolopulos 16 (rotational symmetry)", *_read("rot16.json")
    )


def papanicolopulos_rot_17():
    return TriangleScheme(
        "Papanicolopulos 17 (rotational symmetry)", *_read("rot17.json")
    )


PapanicolopulosRot = {
    8: papanicolopulos_rot_8,
    9: papanicolopulos_rot_9,
    10: papanicolopulos_rot_10,
    11: papanicolopulos_rot_11,
    12: papanicolopulos_rot_12,
    13: papanicolopulos_rot_13,
    14: papanicolopulos_rot_14,
    15: papanicolopulos_rot_15,
    16: papanicolopulos_rot_16,
    17: papanicolopulos_rot_17,
}

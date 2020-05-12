import json
import os

from ...helpers import article
from .._helpers import T2Scheme, untangle2

citation = article(
    authors=["Stefanos-Aldo Papanicolopulos"],
    title="New fully symmetric and rotationally symmetric cubature rules on the triangle using minimal orthonormal bases",
    journal="Journal of Computational and Applied Mathematics",
    volume="294",
    month="mar",
    year="2016",
    pages="39–48",
    url="https://doi.org/10.1016/j.cam.2015.08.001",
)
# https://arxiv.org/abs/1411.5631


def _read(filename):
    this_dir = os.path.dirname(os.path.realpath(__file__))
    with open(os.path.join(this_dir, filename), "r") as f:
        data = json.load(f)

    degree = data.pop("degree")
    points, weights = untangle2(data)
    return weights, points, degree, citation


def papanicolopulos_sym_0():
    return T2Scheme("Papanicolopulos 0 (full symmetry)", *_read("full00.json"))


def papanicolopulos_sym_1():
    return T2Scheme("Papanicolopulos 1 (full symmetry)", *_read("full01.json"))


def papanicolopulos_sym_2():
    return T2Scheme("Papanicolopulos 2 (full symmetry)", *_read("full02.json"))


def papanicolopulos_sym_3():
    return T2Scheme("Papanicolopulos 3 (full symmetry)", *_read("full03.json"))


def papanicolopulos_sym_4():
    return T2Scheme("Papanicolopulos 4 (full symmetry)", *_read("full04.json"))


def papanicolopulos_sym_5():
    return T2Scheme("Papanicolopulos 5 (full symmetry)", *_read("full05.json"))


def papanicolopulos_sym_6():
    return T2Scheme("Papanicolopulos 6 (full symmetry)", *_read("full06.json"))


def papanicolopulos_sym_7():
    return T2Scheme("Papanicolopulos 7 (full symmetry)", *_read("full07.json"))


def papanicolopulos_sym_8():
    return T2Scheme("Papanicolopulos 8 (full symmetry)", *_read("full08.json"))


# TODO ERR the first 8 schemes are flawed by round-off error
# def papanicolopulos_rot_00():
#     return T2Scheme("Papanicolopulos 0 (rotational symmetry)", *_read("rot00.json"))
#
#
# def papanicolopulos_rot_01():
#     return T2Scheme("Papanicolopulos 1 (rotational symmetry)", *_read("rot01.json"))
#
#
# def papanicolopulos_rot_02():
#     return T2Scheme("Papanicolopulos 2 (rotational symmetry)", *_read("rot02.json"))
#
#
# def papanicolopulos_rot_03():
#     return T2Scheme("Papanicolopulos 3 (rotational symmetry)", *_read("rot03.json"))
#
#
# def papanicolopulos_rot_04():
#     return T2Scheme("Papanicolopulos 4 (rotational symmetry)", *_read("rot04.json"))
#
#
# def papanicolopulos_rot_05():
#     return T2Scheme("Papanicolopulos 5 (rotational symmetry)", *_read("rot05.json"))
#
#
# def papanicolopulos_rot_06():
#     return T2Scheme("Papanicolopulos 6 (rotational symmetry)", *_read("rot06.json"))
#
#
# def papanicolopulos_rot_07():
#     return T2Scheme("Papanicolopulos 7 (rotational symmetry)", *_read("rot07.json"))


def papanicolopulos_rot_08():
    return T2Scheme("Papanicolopulos 8 (rotational symmetry)", *_read("rot08.json"))


def papanicolopulos_rot_09():
    return T2Scheme("Papanicolopulos 9 (rotational symmetry)", *_read("rot09.json"))


def papanicolopulos_rot_10():
    return T2Scheme("Papanicolopulos 10 (rotational symmetry)", *_read("rot10.json"))


def papanicolopulos_rot_11():
    return T2Scheme("Papanicolopulos 11 (rotational symmetry)", *_read("rot11.json"))


def papanicolopulos_rot_12():
    return T2Scheme("Papanicolopulos 12 (rotational symmetry)", *_read("rot12.json"))


def papanicolopulos_rot_13():
    return T2Scheme("Papanicolopulos 13 (rotational symmetry)", *_read("rot13.json"))


def papanicolopulos_rot_14():
    return T2Scheme("Papanicolopulos 14 (rotational symmetry)", *_read("rot14.json"))


def papanicolopulos_rot_15():
    return T2Scheme("Papanicolopulos 15 (rotational symmetry)", *_read("rot15.json"))


def papanicolopulos_rot_16():
    return T2Scheme("Papanicolopulos 16 (rotational symmetry)", *_read("rot16.json"))


def papanicolopulos_rot_17():
    return T2Scheme("Papanicolopulos 17 (rotational symmetry)", *_read("rot17.json"))

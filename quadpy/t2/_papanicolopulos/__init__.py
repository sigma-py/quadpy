import pathlib

from ...helpers import article
from .._helpers import _read, register

source = article(
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

this_dir = pathlib.Path(__file__).resolve().parent


def papanicolopulos_sym_0():
    return _read(this_dir / "full00.json", source)


def papanicolopulos_sym_1():
    return _read(this_dir / "full01.json", source)


def papanicolopulos_sym_2():
    return _read(this_dir / "full01.json", source)


def papanicolopulos_sym_3():
    return _read(this_dir / "full01.json", source)


def papanicolopulos_sym_4():
    return _read(this_dir / "full01.json", source)


def papanicolopulos_sym_5():
    return _read(this_dir / "full01.json", source)


def papanicolopulos_sym_6():
    return _read(this_dir / "full01.json", source)


def papanicolopulos_sym_7():
    return _read(this_dir / "full01.json", source)


def papanicolopulos_sym_8():
    return _read(this_dir / "full08.json", source)


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
    return _read(this_dir / "rot08.json", source)


def papanicolopulos_rot_09():
    return _read(this_dir / "rot09.json", source)


def papanicolopulos_rot_10():
    return _read(this_dir / "rot10.json", source)


def papanicolopulos_rot_11():
    return _read(this_dir / "rot11.json", source)


def papanicolopulos_rot_12():
    return _read(this_dir / "rot12.json", source)


def papanicolopulos_rot_13():
    return _read(this_dir / "rot13.json", source)


def papanicolopulos_rot_14():
    return _read(this_dir / "rot14.json", source)


def papanicolopulos_rot_15():
    return _read(this_dir / "rot15.json", source)


def papanicolopulos_rot_16():
    return _read(this_dir / "rot16.json", source)


def papanicolopulos_rot_17():
    return _read(this_dir / "rot17.json", source)


register(
    [
        papanicolopulos_sym_0,
        papanicolopulos_sym_1,
        papanicolopulos_sym_2,
        papanicolopulos_sym_3,
        papanicolopulos_sym_4,
        papanicolopulos_sym_5,
        papanicolopulos_sym_6,
        papanicolopulos_sym_7,
        papanicolopulos_sym_8,
        #
        papanicolopulos_rot_08,
        papanicolopulos_rot_09,
        papanicolopulos_rot_10,
        papanicolopulos_rot_11,
        papanicolopulos_rot_12,
        papanicolopulos_rot_13,
        papanicolopulos_rot_14,
        papanicolopulos_rot_15,
        papanicolopulos_rot_16,
        papanicolopulos_rot_17,
    ]
)

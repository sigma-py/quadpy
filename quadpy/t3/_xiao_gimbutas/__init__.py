import json
import os

import numpy

from ...helpers import article
from .._helpers import T3Scheme

citation = article(
    authors=["Hong Xiao", "Zydrunas Gimbutas"],
    title="A numerical algorithm for the construction of efficient quadrature rules in two and higher dimensions",
    journal="Computers & Mathematics with Applications",
    volume="59",
    number="2",
    month="jan",
    year="2010",
    pages="663–676",
    url="https://doi.org/10.1016/j.camwa.2009.10.027",
)
# Data extracted from
# https://people.sc.fsu.edu/~jburkardt/f_src/triangle_symq_rule/triangle_symq_rule.f90


def _read(degree):
    this_dir = os.path.dirname(os.path.realpath(__file__))
    filename = f"xg{degree:02d}.json"
    with open(os.path.join(this_dir, filename), "r") as f:
        data = json.load(f)

    degree = data.pop("degree")

    points = numpy.array(data["bary"])
    weights = numpy.array(data["weights"])
    return T3Scheme(f"Xiao-Gimbutas {degree}", weights, points, degree, citation)


def xiao_gimbutas_01():
    return _read(1)


def xiao_gimbutas_02():
    return _read(2)


def xiao_gimbutas_03():
    return _read(3)


def xiao_gimbutas_04():
    return _read(4)


def xiao_gimbutas_05():
    return _read(5)


def xiao_gimbutas_06():
    return _read(6)


def xiao_gimbutas_07():
    return _read(7)


def xiao_gimbutas_08():
    return _read(8)


def xiao_gimbutas_09():
    return _read(9)


def xiao_gimbutas_10():
    return _read(10)


def xiao_gimbutas_11():
    return _read(11)


def xiao_gimbutas_12():
    return _read(12)


def xiao_gimbutas_13():
    return _read(13)


def xiao_gimbutas_14():
    return _read(14)


def xiao_gimbutas_15():
    return _read(15)

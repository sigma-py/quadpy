# -*- coding: utf-8 -*-
#
import json
import os

from ...helpers import article
from .._helpers import TetrahedronScheme, untangle2

citation = article(
    authors=["B. Vioreanu", "V. Rokhlin"],
    title="Spectra of Multiplication Operators as a Numerical Tool",
    year="2014",
    journal="SIAM J. Sci. Comput.",
    volume="36",
    number="1",
    url="https://doi.org/10.1137/110860082",
)
# http://www.cs.yale.edu/publications/techreports/tr1443.pdf

# Data adapted from modepy
# https://github.com/inducer/modepy/blob/master/modepy/quadrature/vr_quad_data_tri.py


def _read(index):
    this_dir = os.path.dirname(os.path.realpath(__file__))
    filename = "vr{:02d}.json".format(index)
    with open(os.path.join(this_dir, filename), "r") as f:
        data = json.load(f)

    degree = data.pop("degree")
    points, weights = untangle2(data)
    weights *= 3.0 / 4.0
    return TetrahedronScheme(
        "Vioreanu-Rokhlin {}".format(index), weights, points, degree, citation
    )


def vioreanu_rokhlin_0():
    return _read(0)


def vioreanu_rokhlin_1():
    return _read(1)


def vioreanu_rokhlin_2():
    return _read(2)


def vioreanu_rokhlin_3():
    return _read(3)


def vioreanu_rokhlin_4():
    return _read(4)


def vioreanu_rokhlin_5():
    return _read(5)


def vioreanu_rokhlin_6():
    return _read(6)


def vioreanu_rokhlin_7():
    return _read(7)


def vioreanu_rokhlin_8():
    return _read(8)


def vioreanu_rokhlin_9():
    return _read(9)

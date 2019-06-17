# -*- coding: utf-8 -*-
#
"""
F.D. Witherden, P.E. Vincent,
On the identification of symmetric quadrature rules for finite
element methods,
Computers & Mathematics with Applications,
Volume 69, Issue 10, May 2015, Pages 1232–1241,
<https://doi.org/10.1016/j.camwa.2015.03.017>.

Abstract:
In this paper we describe a methodology for the identification of symmetric quadrature
rules inside of quadrilaterals, triangles, tetrahedra, prisms, pyramids, and hexahedra.
The methodology is free from manual intervention and is capable of identifying a set of
rules with a given strength and a given number of points. We also present polyquad which
is an implementation of our methodology. Using polyquad v1.0 we proceed to derive a
complete set of symmetric rules on the aforementioned domains. All rules possess purely
positive weights and have all points inside the domain. Many of the rules appear to be
new, and an improvement over those tabulated in the literature.
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


def witherden_vincent_1():
    return TriangleScheme("Witherden-Vincent 1", *_read("wv01.json"))


def witherden_vincent_2():
    return TriangleScheme("Witherden-Vincent 2", *_read("wv02.json"))


def witherden_vincent_4():
    return TriangleScheme("Witherden-Vincent 4", *_read("wv04.json"))


def witherden_vincent_5():
    return TriangleScheme("Witherden-Vincent 5", *_read("wv05.json"))


def witherden_vincent_6():
    return TriangleScheme("Witherden-Vincent 6", *_read("wv06.json"))


def witherden_vincent_7():
    return TriangleScheme("Witherden-Vincent 7", *_read("wv07.json"))


def witherden_vincent_8():
    return TriangleScheme("Witherden-Vincent 8", *_read("wv08.json"))


def witherden_vincent_9():
    return TriangleScheme("Witherden-Vincent 9", *_read("wv09.json"))


def witherden_vincent_10():
    return TriangleScheme("Witherden-Vincent 10", *_read("wv10.json"))


def witherden_vincent_11():
    return TriangleScheme("Witherden-Vincent 11", *_read("wv11.json"))


def witherden_vincent_12():
    return TriangleScheme("Witherden-Vincent 12", *_read("wv12.json"))


def witherden_vincent_13():
    return TriangleScheme("Witherden-Vincent 13", *_read("wv13.json"))


def witherden_vincent_14():
    return TriangleScheme("Witherden-Vincent 14", *_read("wv14.json"))


def witherden_vincent_15():
    return TriangleScheme("Witherden-Vincent 15", *_read("wv15.json"))


def witherden_vincent_16():
    return TriangleScheme("Witherden-Vincent 16", *_read("wv16.json"))


def witherden_vincent_17():
    return TriangleScheme("Witherden-Vincent 17", *_read("wv17.json"))


def witherden_vincent_18():
    return TriangleScheme("Witherden-Vincent 18", *_read("wv18.json"))


def witherden_vincent_19():
    return TriangleScheme("Witherden-Vincent 19", *_read("wv19.json"))


def witherden_vincent_20():
    return TriangleScheme("Witherden-Vincent 20", *_read("wv20.json"))


WitherdenVincent = {
    1: witherden_vincent_1,
    2: witherden_vincent_2,
    4: witherden_vincent_4,
    5: witherden_vincent_5,
    6: witherden_vincent_6,
    7: witherden_vincent_7,
    8: witherden_vincent_8,
    9: witherden_vincent_9,
    10: witherden_vincent_10,
    11: witherden_vincent_11,
    12: witherden_vincent_12,
    13: witherden_vincent_13,
    14: witherden_vincent_14,
    15: witherden_vincent_15,
    16: witherden_vincent_16,
    17: witherden_vincent_17,
    18: witherden_vincent_18,
    19: witherden_vincent_19,
    20: witherden_vincent_20,
}

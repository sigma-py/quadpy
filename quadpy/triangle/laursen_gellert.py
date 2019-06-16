# -*- coding: utf-8 -*-
#
"""
M.E. Laursen, M. Gellert,
Some criteria for numerically integrated matrices and quadrature formulas for
triangles,
International Journal for Numerical Methods in Engineering,
Volume 12, Issue 1, 1978, Pages 67–76.
DOI: 10.1002/nme.1620120107,
<https://doi.org/10.1002/nme.1620120107>.

Abstract:
For a wide class of finite element matrices integrated numerically rather than
exactly, a definable number of sampling points is found to be sufficient for
keeping their theoretical properties unchanged. A systematic criterion limiting
the number of possible point configurations for numerical quadrature formulas
on triangles is established. Some new high order formulas are presented. Tables
containing optimal formulas with respect to minimum number of sampling points
and required degrees of accuracy are given. They are arranged so as to assist
with selection of suitable quadrature formulas for finite element computer
programming.
"""
from __future__ import division

import sympy

from .helpers import untangle2


def _gen1(symbolic):
    data = {"s3": [1]}
    return 1, data


def _gen2a(symbolic):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    data = {"s2": [frac(1, 3), frac(1, 6)]}
    return 2, data


def _gen2b(symbolic):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    data = {"s2": [frac(1, 3), frac(1, 2)]}
    return 2, data


def _gen3(symbolic):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    data = {"s3": [-frac(9, 16)], "s2": [frac(25, 48), frac(1, 5)]}
    return 3, data


def _gen4(symbolic):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    data = {"s1": [frac(1, 6), 0.659027622374092, 0.231933368553031]}
    return 3, data


def _gen5(symbolic):
    data = {
        "s2": [
            [0.109951743655322, 0.091576213509771],
            [0.223381589678011, 0.445948490915965],
        ]
    }
    return 4, data


def _gen6(symbolic):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    data = {
        "s3": [frac(3, 8)],
        "s1": [frac(5, 48), 0.736712498968435, 0.237932366472434],
    }
    return 4, data


def _gen7(symbolic):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    data = {
        "s3": [frac(9, 40)],
        "s2": [
            [0.125939180544827, 0.101286507323456],
            [0.132394152788506, 0.470142064105115],
        ],
    }
    return 5, data


def _gen8(symbolic):
    data = {
        "s2": [[0.205950504760887, 0.437525248383384]],
        "s1": [[0.063691414286223, 0.797112651860071, 0.165409927389841]],
    }
    return 5, data


def _gen9(symbolic):
    data = {
        "s2": [
            [0.050844906370207, 0.063089014491502],
            [0.116786275726379, 0.249286745170910],
        ],
        "s1": [[0.082851075618374, 0.636502499121399, 0.310352451033785]],
    }
    return 6, data


def _gen10(symbolic):
    data = {
        "s3": [[-0.149570044467670]],
        "s2": [
            [+0.175615257433204, 0.260345966079038],
            [+0.053347235608839, 0.065130102902216],
        ],
        "s1": [[+0.077113760890257, 0.638444188569809, 0.312865496004875]],
    }
    return 7, data


def _gen11(symbolic):
    data = {
        "s2": [[0.053077801790233, 0.064930513159165]],
        "s1": [
            [0.070853083692136, 0.284575584249173, 0.517039939069325],
            [0.069274682079415, 0.313559184384932, 0.043863471792371],
        ],
    }
    return 7, data


def _gen12(symbolic):
    data = {
        "s3": [[0.144315607677787]],
        "s2": [
            [0.103217370534718, 0.170569307751761],
            [0.032458497623198, 0.050547228317031],
            [0.095091634267284, 0.459292588292723],
        ],
        "s1": [[0.027230314174435, 0.008394777409958, 0.263112829634638]],
    }
    return 8, data


def _gen13(symbolic):
    data = {
        "s3": [[0.097135796282799]],
        "s2": [
            [0.031334700227139, 0.489682519198738],
            [0.077827541004774, 0.437089591492937],
            [0.079647738927210, 0.188203535619033],
            [0.025577675658698, 0.044729513394453],
        ],
        "s1": [[0.043283539377289, 0.036838412054736, 0.221962989160766]],
    }
    return 9, data


def _gen14(symbolic):
    data = {
        "s2": [
            [0.051617202569021, 0.481519834783311],
            [0.094080073458356, 0.403603979817940],
            [0.025993571032320, 0.045189009784377],
        ],
        "s1": [
            [0.045469538047619, 0.136991201264904, 0.218290070971381],
            [0.035351705089199, 0.030424361728820, 0.222063165537318],
        ],
    }
    return 9, data


def _gen15a(symbolic):
    data = {
        "s3": [[0.079894504741240]],
        "s2": [
            [0.071123802232377, 0.425086210602091],
            [0.008223818690464, 0.023308867510000],
        ],
        "s1": [
            [0.045430592296170, 0.147925626209534, 0.223766973576973],
            [0.037359856234305, 0.029946031954171, 0.358740141864431],
            [0.030886656884564, 0.035632559587504, 0.143295370426867],
        ],
    }
    return 10, data


def _gen15b(symbolic):
    data = {
        "s3": [[0.081743329146286]],
        "s2": [
            [0.045957963604745, 0.142161101056564],
            [0.013352968813150, 0.032055373216944],
        ],
        "s1": [
            [0.063904906396424, 0.148132885783821, 0.321812995288835],
            [0.034184648162959, 0.029619889488730, 0.369146781827811],
            [0.025297757707288, 0.028367665339938, 0.163701733737182],
        ],
    }
    return 10, data


_gen = {
    "1": _gen1,
    "2a": _gen2a,
    "2b": _gen2b,
    "3": _gen3,
    "4": _gen4,
    "5": _gen5,
    "6": _gen6,
    "7": _gen7,
    "8": _gen8,
    "9": _gen9,
    "10": _gen10,
    "11": _gen11,
    "12": _gen12,
    "13": _gen13,
    "14": _gen14,
    "15a": _gen15a,
    "15b": _gen15b,
}


class LaursenGellert(object):
    keys = _gen.keys()

    def __init__(self, key, symbolic=False):
        self.name = "Laursen-Gellert({})".format(key)
        self.degree, data = _gen[key](symbolic)
        self.bary, self.weights = untangle2(data)
        return

# -*- coding: utf-8 -*-
#
from __future__ import division

import sympy

from ._helpers import TriangleScheme, concat, s1, s2, s3
from ..helpers import book, article

citation = book(
    authors=["Gilbert Strang", "George Fix"],
    title="An Analysis of the Finite Element Method",
    publisher="Wellesley-Cambridge Press",
    year="1973",
    isbn="096140888X",
    url="https://bookstore.siam.org/wc08/",
)

citation = article(
    authors=["G.R. Cowper"],
    title="Gaussian quadrature formulas for triangles",
    journal="Numerical Methods in Engineering",
    volume="7",
    number="3",
    year="1973",
    pages="405–408",
    url="https://doi.org/10.1002/nme.1620070316",
)
# See also
# https://people.sc.fsu.edu/~jburkardt/datasets/quadrature_rules_tri/quadrature_rules_tri.html


def strang_fix_cowper_01(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y

    weights, points = s2([frac(1, 3), frac(1, 6)])
    return TriangleScheme("Strang 1", weights, points, 2, citation)


def strang_fix_cowper_02(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y

    weights, points = s2([frac(1, 3), frac(1, 2)])
    return TriangleScheme("Strang 2", weights, points, 2, citation)


def strang_fix_cowper_03(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y

    weights, points = concat(s3(-frac(9, 16)), s2([frac(25, 48), frac(1, 5)]))
    return TriangleScheme("Strang 3", weights, points, 3, citation)


def strang_fix_cowper_04(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y

    weights, points = s1([frac(1, 6), 0.659027622374092, 0.231933368553031])
    return TriangleScheme("Strang 4", weights, points, 3, citation)


def strang_fix_cowper_05():
    weights, points = s2(
        [0.109951743655322, 0.091576213509771], [0.223381589678011, 0.445948490915965]
    )
    return TriangleScheme("Strang 5", weights, points, 4, citation)


def strang_fix_cowper_06(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y

    weights, points = concat(
        s3(frac(3, 8)), s1([frac(5, 48), 0.736712498968435, 0.237932366472434])
    )
    return TriangleScheme("Strang 6", weights, points, 4, citation)


def strang_fix_cowper_07(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y

    weights, points = concat(
        s3(frac(9, 40)),
        s2(
            [0.12593918054482717, 0.10128650732345633],
            [0.13239415278850616, 0.47014206410511505],
        ),
    )
    return TriangleScheme("Strang 7", weights, points, 5, citation)


def strang_fix_cowper_08(symbolic=False):
    weights, points = concat(
        s2([0.205950504760887, 0.437525248383384]),
        s1([0.063691414286223, 0.797112651860071, 0.165409927389841]),
    )
    return TriangleScheme("Strang 8", weights, points, 5, citation)


def strang_fix_cowper_09(symbolic=False):
    weights, points = concat(
        s2(
            [0.050844906370207, 0.063089014491502],
            [0.116786275726379, 0.249286745170910],
        ),
        s1([0.082851075618374, 0.636502499121399, 0.310352451033785]),
    )
    return TriangleScheme("Strang 9", weights, points, 6, citation)


def strang_fix_cowper_10():
    weights, points = concat(
        s3(-0.149570044467670),
        s2(
            [+0.175615257433204, 0.260345966079038],
            [+0.053347235608839, 0.065130102902216],
        ),
        s1([+0.077113760890257, 0.638444188569809, 0.312865496004875]),
    )
    return TriangleScheme("Strang 10", weights, points, 7, citation)

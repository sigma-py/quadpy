# -*- coding: utf-8 -*-
#
"""
Ann Haegemans and Robert Piessens,
Construction of Cubature Formulas of Degree Seven and Nine Symmetric Planar Regions
Using Orthogonal Polynomials,
SIAM Journal on Numerical Analysis, Vol. 14, No. 3 (Jun., 1977), pp. 492-508,
<https://www.jstor.org/stable/2156699>.

Abstract:
A method of constructing twelve-point cubature formulas with polynomial precision seven
and nineteen-point and eighteen-point cubature formulas with polynomial precision nine
is given for planar regions and weight functions, which are symmetric in each variable.
The nodes are computed as common zeros of a set of linearly independent orthogonal
polynomials.
"""
from .helpers import QuadrilateralScheme, concat, pm, pm2


def HaegemansPiessens():
    weights, points = concat(
        pm2(
            [0.213057211620949126, 0.9171178223127705862, 0.547931206828092323],
            [0.17400948894689560610, 0.61126876646532841440, 0.93884325665885830459],
        ),
        pm(
            [0.63585388344327977182, +0.52942280204265532589, 0.0],
            [0.59001271542103076297, 0.0, +0.62704137378039531763],
        ),
    )
    return QuadrilateralScheme("Haegemans-Piessens", 7, weights, points)

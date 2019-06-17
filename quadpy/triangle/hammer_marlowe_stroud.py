# -*- coding: utf-8 -*-
#
"""
P.C. Hammer, O.J. Marlowe and A.H. Stroud,
Numerical Integration Over Simplexes and Cones,
Mathematical Tables and Other Aids to Computation,
Vol. 10, No. 55, Jul. 1956, pp. 130-137,
<https://doi.org/10.1090/S0025-5718-1956-0086389-6>.

Abstract:
In this paper we develop numerical integration formulas for simplexes and cones in
n-space for n>=2. While several papers have been written on numerical integration in
higher spaces, most of these have dealt with hyperrectangular regions. For certain
exceptions see [3]. Hammer and Wymore [1] have given a first general type theory
designed through systematic use of cartesian product regions and affine transformations
to extend the possible usefulness of formulas for each region.

Two of the schemes also appear in

P.C. Hammer, Arthur H. Stroud,
Numerical Evaluation of Multiple Integrals II,
Mathematical Tables and Other Aids to Computation.
Vol. 12, No. 64 (Oct., 1958), pp. 272-280,
<https://www.jstor.org/stable/2002370>
"""
from __future__ import division

import numpy
import sympy

from .helpers import TriangleScheme, s3, r, concat


def hammer_marlowe_stroud_1(symbolic=False):
    weights, bary = s3(1)
    return TriangleScheme("Hammer-Marlowe-Stroud 1", 1, weights, bary)


def hammer_marlowe_stroud_2(symbolic=False):
    frac = sympy.frac if symbolic else lambda x, y: x / y

    weights, bary = r([frac(1, 3), frac(1, 2)])
    return TriangleScheme("Hammer-Marlowe-Stroud 1", 2, weights, bary)


def hammer_marlowe_stroud_3(symbolic=False):
    frac = sympy.frac if symbolic else lambda x, y: x / y

    weights, bary = r([frac(1, 3), -frac(1, 2)])
    return TriangleScheme("Hammer-Marlowe-Stroud 3", 2, weights, bary)


def hammer_marlowe_stroud_4(symbolic=False):
    frac = sympy.frac if symbolic else lambda x, y: x / y

    weights, bary = concat(s3(-frac(9, 16)), r([frac(25, 48), frac(2, 5)]))
    return TriangleScheme("Hammer-Marlowe-Stroud 4", 3, weights, bary)


def hammer_marlowe_stroud_5(symbolic=False):
    frac = sympy.frac if symbolic else lambda x, y: x / y
    sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt

    w1, w2 = [(155 - i * sqrt(15)) / 1200 for i in [+1, -1]]
    x1, x2 = [(1 + i * sqrt(15)) / 7 for i in [+1, -1]]
    weights, bary = concat(s3(frac(9, 40)), r([w1, x1], [w2, x2]))
    return TriangleScheme("Hammer-Marlowe-Stroud 5", 5, weights, bary)


HammerMarloweStroud = {
    1: hammer_marlowe_stroud_1,
    2: hammer_marlowe_stroud_2,
    3: hammer_marlowe_stroud_3,
    4: hammer_marlowe_stroud_4,
    5: hammer_marlowe_stroud_5,
}

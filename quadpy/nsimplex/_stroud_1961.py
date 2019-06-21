# -*- coding: utf-8 -*-
#
from __future__ import division

import sympy

from ..helpers import untangle, rd, article
from ._helpers import NSimplexScheme


citation = article(
    authors=["A.H. Stroud"],
    title="Numerical Integration Formulas of Degree 3 for Product Regions and Cones",
    journal="Mathematics of Computation",
    volume="15",
    number="74",
    month="apr",
    year="1961",
    pages="143-150",
    url="https://doi.org/10.2307/2004220",
)


def stroud_1961(n, symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y

    degree = 3

    r = frac(1, n + 1)
    s = frac(1, n)

    prod = (n + 1) * (n + 2) * (n + 3)
    A = frac((3 - n) * (n + 1) ** 3, prod)
    B = frac(3, prod)
    C = frac(n ** 3, prod)

    data = [(A, [(n + 1) * [r]]), (B, rd(n + 1, [(1, 1)])), (C, rd(n + 1, [(s, n)]))]

    bary, weights = untangle(data)
    return NSimplexScheme("Stroud 1961", n, weights, bary, degree, citation)

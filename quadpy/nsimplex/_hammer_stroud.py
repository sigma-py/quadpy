# -*- coding: utf-8 -*-
#
from __future__ import division

import numpy
import sympy

from ..helpers import untangle, rd, article
from ._helpers import NSimplexScheme

citation = article(
    authors=["Preston C. Hammer", "Arthur H. Stroud"],
    title="Numerical Integration Over Simplexes",
    journal="Mathematical Tables and Other Aids to Computation",
    volume="10",
    number="55",
    month="jul",
    year="1956",
    pages="137-139",
    url="https://doi.org/10.2307/2002484",
)


def hammer_stroud_1a(n, symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt

    degree = 2
    r = (n + 2 - sqrt(n + 2)) / (n + 1) / (n + 2)
    s = (n + 2 + n * sqrt(n + 2)) / (n + 1) / (n + 2)
    data = [(frac(1, n + 1), rd(n + 1, [(r, n), (s, 1)]))]

    bary, weights = untangle(data)
    return NSimplexScheme("Hammer-Stround 1a", n, weights, bary, degree, citation)


def hammer_stroud_1b(n, symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt

    degree = 2
    r = (n + 2 + sqrt(n + 2)) / (n + 1) / (n + 2)
    s = (n + 2 - n * sqrt(n + 2)) / (n + 1) / (n + 2)
    data = [(frac(1, n + 1), rd(n + 1, [(r, n), (s, 1)]))]

    bary, weights = untangle(data)
    return NSimplexScheme("Hammer-Stround 1b", n, weights, bary, degree, citation)


def hammer_stroud_2(n, symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y

    degree = 3

    B = -frac((n + 1) ** 2, 4 * (n + 2))
    C = frac((n + 3) ** 2, 4 * (n + 1) * (n + 2))

    r = frac(1, n + 1)
    s = frac(1, n + 3)
    t = frac(3, n + 3)

    data = [(B, [(n + 1) * [r]]), (C, rd(n + 1, [(t, 1), (s, n)]))]

    bary, weights = untangle(data)
    return NSimplexScheme("Hammer-Stround 2", n, weights, bary, degree, citation)

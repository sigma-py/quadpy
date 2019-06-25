# -*- coding: utf-8 -*-
#

from mpmath import mp
import numpy
import sympy

from ..helpers import untangle, rd, article
from ._helpers import NSimplexScheme

citation = article(
    authors=["A.H. Stroud"],
    title="Approximate Integration Formulas of Degree 3 for Simplexes",
    journal="Mathematics of Computation",
    volume="18",
    number="88",
    month="oct",
    year="1964",
    pages="590-597",
    url="https://doi.org/10.2307/2002945",
)


def _stroud_1964(variant_a, n, symbolic=False):

    frac = sympy.Rational if symbolic else lambda x, y: x / y
    roots = mp.polyroots if symbolic else numpy.roots

    degree = 3

    if n == 2:
        # The roots sum up to 1
        r, s, t = mp.polyroots([1, -1, frac(1, 4), -frac(1, 60)])
        data = [(frac(1, n * (n + 1)), rd(n + 1, [(r, 1), (s, 1), (t, 1)]))]
    else:
        assert n > 2

        # Stroud's book only gives numerical values for certain n; the article explains
        # it in more detail, namely: r is a root of a polynomial of degree 3.
        rts = numpy.sort(
            roots([n + 1, -3, frac(3, n + 2), -frac(1, (n + 2) * (n + 3))])
        )

        # all roots are real-valued
        if n > 8:
            assert not variant_a, "Choose variant b for n >= 9."

        r = rts[0] if variant_a else rts[1]

        # s and t are zeros of a polynomial of degree 2
        s, t = numpy.sort(
            roots(
                [
                    1,
                    -(1 - (n - 1) * r),
                    frac(n, 2 * (n + 2)) - (n - 1) * r + frac(n * (n - 1), 2) * r ** 2,
                ]
            )
        )

        data = [(frac(1, n * (n + 1)), rd(n + 1, [(r, n - 1), (s, 1), (t, 1)]))]

    points, weights = untangle(data)

    name = "Stroud 1964{}".format("a" if variant_a else "b")
    return NSimplexScheme(name, n, weights, points, degree, citation)


def stroud_1964a(n, symbolic=False):
    return _stroud_1964(True, n, symbolic)


def stroud_1964b(n, symbolic=False):
    return _stroud_1964(False, n, symbolic)

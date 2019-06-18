# -*- coding: utf-8 -*-
#
from __future__ import division

import numpy
import sympy

from ._helpers import E2rScheme
from ..helpers import untangle, pm_array, pm, fsd, article


_citation = article(
    authors=["A.H. Stroud", "D. Secrest"],
    title="Approximate integration formulas for certain spherically symmetric regions",
    journal="Math. Comp.",
    volume="17",
    year="1963",
    pages="105-135",
    url="https://doi.org/10.1090/S0025-5718-1963-0161473-0",
)


def stroud_secrest_v(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt
    pi = sympy.pi if symbolic else numpy.pi

    nu = 2 * sqrt(5)
    xi = sqrt(5)
    eta = sqrt(15)

    data = [
        (frac(7, 10), numpy.array([[0, 0]])),
        (frac(1, 20), numpy.array([[+nu, 0], [-nu, 0]])),
        (frac(1, 20), pm_array([xi, eta])),
    ]
    points, weights = untangle(data)
    weights *= 2 * pi
    return E2rScheme("Stroud-Secrest V", weights, points, 5, _citation)


def stroud_secrest_vi(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt
    pi = sympy.pi if symbolic else numpy.pi

    sqrt74255 = sqrt(74255)

    nu = sqrt(42)
    xi, eta = [sqrt((6615 - p_m * 21 * sqrt74255) / 454) for p_m in [+1, -1]]
    A = frac(5, 588)
    B, C = [(5272105 + p_m * 18733 * sqrt74255) / 43661940 for p_m in [+1, -1]]

    data = [(A, fsd(2, (nu, 1))), (B, pm(2, xi)), (C, pm(2, eta))]

    points, weights = untangle(data)
    weights *= 2 * pi
    return E2rScheme("Stroud-Secrest VI", weights, points, 7, _citation)

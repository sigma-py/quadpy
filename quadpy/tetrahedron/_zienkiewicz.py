# -*- coding: utf-8 -*-
#
from __future__ import division

import sympy

from ._helpers import untangle2, TetrahedronScheme
from ..helpers import book

citation = book(
    authors=["Olgierd Zienkiewicz"],
    title="The Finite Element Method, Sixth Edition",
    publisher="Butterworth-Heinemann",
    year="2005",
    isbn="0750663200",
    url="http://www.sciencedirect.com/science/book/9780750664318",
)
# https://people.sc.fsu.edu/~jburkardt/datasets/quadrature_rules_tet/quadrature_rules_tet.html


def zienkiewicz_4(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y

    degree = 2
    data = {"s31": [[frac(1, 4), 0.1381966011250105]]}
    bary, weights = untangle2(data)
    return TetrahedronScheme("Zienkiewicz 4", weights, bary, degree, citation)


def zienkiewicz_5(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y

    degree = 3
    data = {"s4": [[-frac(4, 5)]], "s31": [[frac(9, 20), frac(1, 6)]]}
    bary, weights = untangle2(data)
    return TetrahedronScheme("Zienkiewicz 5", weights, bary, degree, citation)

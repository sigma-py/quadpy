# -*- coding: utf-8 -*-
#
from __future__ import division

import sympy

from .helpers import untangle2


class HammerStroud(object):
    """
    Preston C. Hammer and Arthur H. Stroud,
    Numerical Evaluation of Multiple Integrals II,
    Math. Comp. 12 (1958), 272-280,
    <https://doi.org/10.1090/S0025-5718-1958-0102176-6>.
    """

    def __init__(self, degree, symbolic=False):
        frac = sympy.Rational if symbolic else lambda x, y: x / y

        self.name = "Hammer-Stroud"
        self.degree = degree
        if degree == 2:
            data = {"s2": [[frac(1, 3), frac(1, 6)]]}
        else:
            assert degree == 3
            data = {"s3": [[-frac(27, 48)]], "s2": [[+frac(25, 48), frac(1, 5)]]}

        self.bary, self.weights = untangle2(data)
        self.points = self.bary[:, 1:]
        return

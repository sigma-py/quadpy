# -*- coding: utf-8 -*-
#
from __future__ import division

import sympy

from .helpers import _symm_r_0, _symm_s, _z, _symm_s_t
from ..helpers import untangle


class Meister(object):
    """
    Bernd Meister,
    On a Family of Cubature Formulae,
    Comput J (1966) 8 (4): 368-371,
    <https://doi.org/10.1093/comjnl/8.4.368>.
    """

    def __init__(self, symbolic=False):
        frac = sympy.Rational if symbolic else lambda x, y: x / y

        self.name = "Meister"
        self.degree = 7

        r = frac(2, 3)
        s = frac(1, 3)

        data = [
            (frac(1024, 6720), _z()),
            (frac(576, 6720), _symm_s(r)),
            (frac(576, 6720), _symm_r_0(r)),
            (-frac(9, 6720), _symm_s(s)),
            (frac(117, 6720), _symm_s_t(1, s)),
            (frac(47, 6720), _symm_s(1)),
        ]

        self.points, self.weights = untangle(data)
        self.weights *= 4
        return

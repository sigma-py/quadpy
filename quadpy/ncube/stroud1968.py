# -*- coding: utf-8 -*-
#
from __future__ import division

import numpy
import sympy

from .helpers import _s2, _s11
from ..helpers import untangle, fsd, z


class Stroud1968(object):
    """
    A.H. Stroud,
    Extensions of Symmetric Integration Formulas,
    Mathematics of Computation,
    Vol. 22, No. 102 (Apr., 1968), pp. 271-274,
    Published by: American Mathematical Society,
    <https://doi.org/10.2307/2004655>.
    """

    def __init__(self, n, symbolic=False):
        frac = sympy.Rational if symbolic else lambda x, y: x / y
        sqrt = sympy.sqrt if symbolic else numpy.sqrt

        self.degree = 5

        r = sqrt(frac(7, 15))
        s, t = [sqrt((7 + i * sqrt(24)) / 15) for i in [+1, -1]]
        data = [
            (frac(5 * n ** 2 - 15 * n + 14, 14), z(n)),
            (frac(25, 168), _s2(n, +r)),
            (frac(25, 168), _s2(n, -r)),
            (frac(-25 * (n - 2), 168), fsd(n, (r, 1))),
            (frac(5, 48), _s11(n, +s, -t)),
            (frac(5, 48), _s11(n, -s, +t)),
            (frac(-5 * (n - 2), 48), fsd(n, (s, 1))),
            (frac(-5 * (n - 2), 48), fsd(n, (t, 1))),
        ]

        self.points, self.weights = untangle(data)
        reference_volume = 2 ** n
        self.weights *= reference_volume
        return

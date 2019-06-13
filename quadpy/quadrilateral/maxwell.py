# -*- coding: utf-8 -*-
#
from __future__ import division

import numpy
import sympy

from .helpers import concat, zero, symm_r0, symm_s_t


class Maxwell(object):
    """
    J.C. Maxwell,
    On Approximate Multiple Integration between Limits by Summation.
    In W. Niven (Ed.), The Scientific Papers of James Clerk Maxwell,
    Cambridge Library Collection - Physical Sciences, pp. 604-611.
    Cambridge: Cambridge University Press.
    First published in 1890.
    <https://doi.org/10.1017/CBO9780511710377.061>.
    """

    def __init__(self, symbolic=False):
        frac = sympy.Rational if symbolic else lambda x, y: x / y
        sqrt = sympy.sqrt if symbolic else numpy.sqrt

        self.name = "Maxwell"
        self.degree = 7

        r = sqrt(frac(12, 35))
        s, t = [sqrt((93 + i * 3 * sqrt(186)) / 155) for i in [+1, -1]]

        self.weights, self.points = concat(
            zero(frac(1, 81)),
            symm_r0([frac(49, 324), r]),
            # ERR typo in Stroud: 648 vs 649
            symm_s_t([frac(31, 648), s, t]),
        )
        self.weights *= 4
        return

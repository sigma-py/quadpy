# -*- coding: utf-8 -*-
#
from __future__ import division

import sympy

from .helpers import untangle2


class Zienkiewicz(object):
    """
    Olgierd Zienkiewicz,
    The Finite Element Method,
    Sixth Edition,
    Butterworth-Heinemann, 2005,
    ISBN: 0750663200,
    <http://www.sciencedirect.com/science/book/9780750664318>,
    <https://people.sc.fsu.edu/~jburkardt/datasets/quadrature_rules_tet/quadrature_rules_tet.html>.
    """

    def __init__(self, index, symbolic=False):
        frac = sympy.Rational if symbolic else lambda x, y: x / y

        self.name = "Zienkiewicz({})".format(index)

        data = {
            4: {"degree": 2, "s31": [[frac(1, 4), 0.1381966011250105]]},
            5: {"degree": 3, "s4": [[-frac(4, 5)]], "s31": [[frac(9, 20), frac(1, 6)]]},
        }[index]

        self.degree = data.pop("degree")

        self.bary, self.weights = untangle2(data)
        self.points = self.bary[:, 1:]
        return

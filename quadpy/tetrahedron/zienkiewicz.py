# -*- coding: utf-8 -*-
#
from __future__ import division

import sympy

from .helpers import _s4, _s31
from ..helpers import untangle


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
        if index == 4:
            self.degree = 2
            data = [(frac(1, 4), _s31(0.1381966011250105))]
        else:
            assert index == 5
            self.degree = 3
            data = [(-frac(4, 5), _s4()), (frac(9, 20), _s31(frac(1, 6)))]

        self.bary, self.weights = untangle(data)
        self.points = self.bary[:, 1:]
        return

# -*- coding: utf-8 -*-
#
from sympy import Rational as fr

from .helpers import _s4, _s31
from ..helpers import untangle


class Zienkiewicz(object):
    '''
    Olgierd Zienkiewicz,
    The Finite Element Method,
    Sixth Edition,
    Butterworth-Heinemann, 2005,
    ISBN: 0750663200,
    <http://www.sciencedirect.com/science/book/9780750664318>,
    <https://people.sc.fsu.edu/~jburkardt/datasets/quadrature_rules_tet/quadrature_rules_tet.html>.
    '''
    def __init__(self, index):
        if index == 4:
            self.degree = 2
            data = [
                (fr(1, 4), _s31(0.1381966011250105))
                ]
        else:
            assert index == 5
            self.degree = 3
            data = [
                (-fr(4, 5), _s4()),
                (fr(9, 20), _s31(fr(1, 6))),
                ]

        self.bary, self.weights = untangle(data)
        self.points = self.bary[:, 1:]
        return

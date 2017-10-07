# -*- coding: utf-8 -*-
#
from sympy import sqrt, Rational as fr

from .helpers import _s4, _s31
from ..helpers import untangle


class HammerStroud(object):
    '''
    Preston C. Hammer and Arthur H. Stroud,
    Numerical Evaluation of Multiple Integrals II,
    Math. Comp. 12 (1958), 272-280,
    <https://doi.org/10.1090/S0025-5718-1958-0102176-6>.
    '''
    def __init__(self, degree):
        self.degree = degree
        if degree == 2:
            data = [
                (fr(1, 4), _s31((5 - sqrt(5))/20)),
                ]
        else:
            assert degree == 3
            data = [
                (-fr(4, 5), _s4()),
                (+fr(9, 20), _s31(fr(1, 6))),
                ]

        self.bary, self.weights = untangle(data)
        self.points = self.bary[:, 1:]
        return

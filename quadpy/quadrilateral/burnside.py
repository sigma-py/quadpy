# -*- coding: utf-8 -*-
#
from sympy import sqrt, Rational as fr

from .helpers import _symm_r_0, _symm_s
from ..helpers import untangle


class Burnside(object):
    '''
    W. Burnside,
    An approximate quadrature formula,
    Messenger of Math., v. 37, 1908, pp. 166-167.
    '''
    def __init__(self):
        self.name = 'Burnside'
        self.degree = 5
        r = sqrt(fr(7, 15))
        s = sqrt(fr(7, 9))
        data = [
            (fr(10, 49), _symm_r_0(r)),
            (fr(9, 196), _symm_s(s)),
            ]

        self.points, self.weights = untangle(data)
        self.weights *= 4
        return

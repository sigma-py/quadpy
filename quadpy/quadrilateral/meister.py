# -*- coding: utf-8 -*-
#
from sympy import Rational as fr

from .helpers import _symm_r_0, _symm_s, _z, _symm_s_t
from ..helpers import untangle


class Meister(object):
    '''
    Bernd Meister,
    On a Family of Cubature Formulae,
    Comput J (1966) 8 (4): 368-371,
    <https://doi.org/10.1093/comjnl/8.4.368>.
    '''
    def __init__(self):
        self.name = 'Meister'
        self.degree = 7

        r = fr(2, 3)
        s = fr(1, 3)

        data = [
            (fr(1024, 6720), _z()),
            (fr(576, 6720), _symm_s(r)),
            (fr(576, 6720), _symm_r_0(r)),
            (-fr(9, 6720), _symm_s(s)),
            (fr(117, 6720), _symm_s_t(1, s)),
            (fr(47, 6720), _symm_s(1)),
            ]

        self.points, self.weights = untangle(data)
        self.weights *= 4
        return

# -*- coding: utf-8 -*-
#
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

        r = 2.0/3.0
        s = 1.0/3.0

        data = [
            (1024.0/6720.0, _z()),
            (576.0/6720.0, _symm_s(r)),
            (576.0/6720.0, _symm_r_0(r)),
            (-9.0/6720.0, _symm_s(s)),
            (117.0/6720.0, _symm_s_t(1.0, s)),
            (47.0/6720.0, _symm_s(1.0)),
            ]

        self.points, self.weights = untangle(data)
        self.weights *= 4.0
        return

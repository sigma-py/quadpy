# -*- coding: utf-8 -*-
#
from .helpers import fs_rr0, z

from ..helpers import untangle


class MustardLynessBlatt(object):
    '''
    D. Mustard, J.N. Lyness, J.M. Blatt,
    Numerical quadrature in n dimensions,
    Comput J (1963) 6 (1): 75-87,
    <https://doi.org/10.1093/comjnl/6.1.75>.
    '''
    def __init__(self):
        self.degree = 3
        data = [
            (0.5, z()),
            (1.0/24.0, fs_rr0(1.0)),
            ]

        self.points, self.weights = untangle(data)
        self.weights *= 8.0
        return

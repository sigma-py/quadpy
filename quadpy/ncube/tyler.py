# -*- coding: utf-8 -*-
#
from ..helpers import untangle, fsd, z


class Tyler(object):
    '''
    G.W. Tyler,
    Numerical integration of functions of several variables,
    Canad. J. Math. 5(1953), 393-412,
    <https://dx.doi.org/10.4153/CJM-1953-044-1>.
    '''
    def __init__(self, n):
        reference_volume = 2.0**n
        self.degree = 3
        data = [
            ((3.0 - n)/3.0, z(n)),
            (1.0/6.0, fsd(n, 1.0, 1)),
            ]

        self.points, self.weights = untangle(data)
        self.weights *= reference_volume
        return

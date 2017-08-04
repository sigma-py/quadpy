# -*- coding: utf-8 -*-
#
from .helpers import fs00

from ..helpers import untangle


class Tyler(object):
    '''
    G.W. Tyler,
    Numerical integration of functions of several variables,
    Canad. J. Math. 5(1953), 393-412,
    <https://dx.doi.org/10.4153/CJM-1953-044-1>.
    '''
    def __init__(self):
        self.degree = 3
        data = [
            (1.0/6.0, fs00(1.0)),
            ]

        self.points, self.weights = untangle(data)
        self.weights *= 8.0
        return

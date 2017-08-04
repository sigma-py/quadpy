# -*- coding: utf-8 -*-
#
from .helpers import fs_r00, pm_rrr, z

from ..helpers import untangle


class Tyler(object):
    '''
    G.W. Tyler,
    Numerical integration of functions of several variables,
    Canad. J. Math. 5(1953), 393-412,
    <https://dx.doi.org/10.4153/CJM-1953-044-1>.
    '''
    def __init__(self, index):
        if index == 1:
            self.degree = 3
            data = [
                (1.0/6.0, fs_r00(1.0)),
                ]
        else:
            assert index == 2
            self.degree = 5
            data = [
                (-62.0/45.0, z()),
                (16.0/45.0, fs_r00(0.5)),
                (1.0/45.0, fs_r00(1.0)),
                (1.0/72.0, pm_rrr(1.0)),
                ]

        self.points, self.weights = untangle(data)
        self.weights *= 8.0
        return

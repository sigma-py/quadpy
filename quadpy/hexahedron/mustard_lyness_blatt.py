# -*- coding: utf-8 -*-
#
from .helpers import fs_rr0, fs_r00, pm_rrr, z

from ..helpers import untangle


class MustardLynessBlatt(object):
    '''
    D. Mustard, J.N. Lyness, J.M. Blatt,
    Numerical quadrature in n dimensions,
    Comput J (1963) 6 (1): 75-87,
    <https://doi.org/10.1093/comjnl/6.1.75>.
    '''
    def __init__(self, index):
        if index == 1:
            self.degree = 3
            data = [
                (0.5, z()),
                (1.0/24.0, fs_rr0(1.0)),
                ]
        elif index == 2:
            self.degree = 3
            data = [
                (2.0/9.0, z()),
                (1.0/9.0, fs_r00(1.0)),
                (1.0/72.0, pm_rrr(1.0))
                ]
        elif index == 3:
            self.degree = 3
            data = [
                (+1.0/6.0, fs_rr0(1.0)),
                (-1.0/8.0, pm_rrr(1.0))
                ]
        elif index == 4:
            self.degree = 5
            data = [
                (-2.0/45.0, z()),
                (+2.0/45.0, fs_r00(1.0)),
                (+4.0/45.0, pm_rrr(0.5)),
                (1.0/120.0, pm_rrr(1.0)),
                ]
        elif index == 5:
            self.degree = 5
            data = [
                (-19.0/15.0, z()),
                (+16.0/45.0, fs_r00(0.5)),
                (-1.0/30.0, fs_r00(1.0)),
                (+1.0/36.0, fs_rr0(1.0)),
                ]
        elif index == 6:
            self.degree = 5
            data = [
                (-4.0/3.0, z()),
                (+16.0/45.0, fs_r00(0.5)),
                (1.0/90.0, fs_rr0(1.0)),
                (1.0/120.0, pm_rrr(1.0)),
                ]
        else:
            assert index == 7
            self.degree = 5
            data = [
                (2.0/45.0, z()),
                (1.0/45.0, fs_rr0(1.0)),
                (4.0/45.0, pm_rrr(0.5)),
                (-1.0/360.0, pm_rrr(1.0)),
                ]

        self.points, self.weights = untangle(data)
        self.weights *= 8.0
        return

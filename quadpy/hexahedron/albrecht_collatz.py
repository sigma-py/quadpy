# -*- coding: utf-8 -*-
#
from .helpers import fs_r00, fs_rr0, z

from ..helpers import untangle


class AlbrechtCollatz(object):
    '''
    J. Albrecht, L. Collatz,
    Zur numerischen Auswertung mehrdimensionaler Integrale,
    ZAMM, Volume 38, Issue 1-2, 1958, Pages 1–15,
    <https://dx.doi.org/10.1002/zamm.19580380102>
    '''
    def __init__(self):
        self.degree = 3
        data = [
            (1.0/4.0, z()),
            (1.0/12.0, fs_r00(1.0)),
            (1.0/48.0, fs_rr0(1.0)),
            ]

        self.points, self.weights = untangle(data)
        self.weights *= 8.0
        return

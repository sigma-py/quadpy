# -*- coding: utf-8 -*-
#
import numpy

from .helpers import fs_r00, pm_rrr

from ..helpers import untangle


class HammerStroud(object):
    '''
    Preston C. Hammer and Arthur H. Stroud,
    Numerical Evaluation of Multiple Integrals II,
    Math. Comp. 12 (1958), 272-280,
    <https://doi.org/10.1090/S0025-5718-1958-0102176-6>.
    '''
    def __init__(self):
        self.degree = 5
        data = [
            (40.0/361.0, fs_r00(numpy.sqrt(19.0/30.0))),
            (121.0/2888.0, pm_rrr(numpy.sqrt(19.0/33.0)))
            ]

        self.points, self.weights = untangle(data)
        self.weights *= 8.0
        return

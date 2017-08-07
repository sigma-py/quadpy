# -*- coding: utf-8 -*-
#
import numpy

from ..helpers import untangle, fsd
from .helpers import volume_unit_ball


class HammerStroud(object):
    '''
    Preston C. Hammer and Arthur H. Stroud,
    Numerical Evaluation of Multiple Integrals II,
    Math. Comp. 12 (1958), 272-280,
    <https://doi.org/10.1090/S0025-5718-1958-0102176-6>.
    '''
    def __init__(self, n, alpha):
        self.degree = 3
        self.dim = n

        r = numpy.sqrt((n + alpha) / (n + alpha + 2))

        data = [
            (1.0/(2*n), fsd(n, r, 1))
            ]

        self.points, self.weights = untangle(data)
        self.weights *= volume_unit_ball(n)
        return

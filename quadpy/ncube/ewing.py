# -*- coding: utf-8 -*-
#
from .helpers import _z, _pm

from ..helpers import untangle


class Ewing(object):
    '''
    G.M. Ewing,
    On Approximate Cubature,
    The American Mathematical Monthly,
    Vol. 48, No. 2 (Feb., 1941), pp. 134-136,
    <https://dx.doi.org/dx.doi.org/10.2307/2303604>.
    '''
    def __init__(self, n):
        reference_volume = 2.0**n
        self.degree = 3
        data = [
            (2.0/3.0, _z(n)),
            (1.0/3.0 / 2**n, _pm(n, 1.0)),
            ]

        self.points, self.weights = untangle(data)
        self.weights *= reference_volume
        return

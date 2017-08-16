# -*- coding: utf-8 -*-
#
import math

from ..helpers import untangle, fsd, z


class HammerStroud(object):
    '''
    Preston C. Hammer and Arthur H. Stroud,
    Numerical Evaluation of Multiple Integrals II,
    Math. Comp. 12 (1958), 272-280,
    <https://doi.org/10.1090/S0025-5718-1958-0102176-6>.
    '''
    def __init__(self, n, index):
        self.dim = n

        if index == '1-n':
            self.degree = 3
            data = [
                (0.5/n, fsd(n, math.sqrt(n/3.0), 1)),
                ]
        else:
            assert index == '2-n'
            self.degree = 5
            r = math.sqrt(3.0 / 5.0)
            data = [
                ((25*n**2 - 115*n + 162)/162.0, z(n)),
                ((70 - 25*n)/162.0, fsd(n, r, 1)),
                (25.0/324.0, fsd(n, r, 2)),
                ]

        self.points, self.weights = untangle(data)
        reference_volume = 2.0**n
        self.weights *= reference_volume
        return

# -*- coding: utf-8 -*-
#
from sympy import sqrt, Rational as fr

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
                (fr(1, 2*n), fsd(n, (sqrt(fr(n, 3)), 1))),
                ]
        else:
            assert index == '2-n'
            self.degree = 5
            r = sqrt(fr(3, 5))
            data = [
                (fr(25*n**2 - 115*n + 162, 162), z(n)),
                (fr(70 - 25*n, 162), fsd(n, (r, 1))),
                (fr(25, 324), fsd(n, (r, 2))),
                ]

        self.points, self.weights = untangle(data)
        reference_volume = 2**n
        self.weights *= reference_volume
        return

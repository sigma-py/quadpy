# -*- coding: utf-8 -*-
#
from sympy import Rational as fr

from ..helpers import untangle, rd


class Stroud1961(object):
    '''
    A.H. Stroud,
    Numerical Integration Formulas of Degree 3 for Product Regions and Cones
    Mathematics of Computation, Vol. 15, No. 74 (Apr., 1961), pp. 143-150,
    <https://doi.org/10.2307/2004220>.
    '''
    def __init__(self, n):
        self.dim = n
        self.degree = 3

        r = fr(1, n+1)
        s = fr(1, n)

        prod = (n+1) * (n+2) * (n+3)
        A = fr((3-n) * (n+1)**3, prod)
        B = fr(3, prod)
        C = fr(n**3, prod)

        data = [
            (A, [(n+1) * [r]]),
            (B, rd(n+1, [(1, 1)])),
            (C, rd(n+1, [(s, n)])),
            ]

        self.bary, self.weights = untangle(data)
        self.points = self.bary[:, 1:]
        return

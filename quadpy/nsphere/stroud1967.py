# -*- coding: utf-8 -*-
#
from sympy import sqrt, Rational as fr

from ..helpers import untangle, fsd, pm
from .helpers import integrate_monomial_over_unit_nsphere


class Stroud1967(object):
    '''
    A.H. Stroud,
    Some Seventh Degree Integration Formulas for Symmetric Regions,
    SIAM J. Numer. Anal., 4(1), 37–44. (8 pages),
    <https://doi.org/10.1137/0704004>.
    '''
    def __init__(self, n):
        self.dim = n
        self.degree = 7

        r = 1
        s = sqrt(fr(1, n))
        t = sqrt(fr(1, 2))
        B = fr(8-n, n * (n+2) * (n+4))
        C = fr(n**3, 2**n * n * (n+2) * (n+4))
        D = fr(4, n * (n+2) * (n+4))

        data = [
            (B, fsd(n, (r, 1))),
            (C, pm(n, s)),
            # ERR Stroud's book wrongly states (t, t,..., t)_FS instead of
            # (t, t, 0, ..., 0)_FS.
            (D, fsd(n, (t, 2))),
            ]

        self.points, self.weights = untangle(data)
        self.weights *= integrate_monomial_over_unit_nsphere(n * [0])
        return

# -*- coding: utf-8 -*-
#
from __future__ import division
import numpy

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

        r = 1.0
        s = numpy.sqrt(1.0 / n)
        t = numpy.sqrt(0.5)
        B = (8 - n) / n / (n+2) / (n+4)
        C = n**3 / 2**n / n / (n+2) / (n+4)
        D = 4 / n / (n+2) / (n+4)

        data = [
            (B, fsd(n, r, 1)),
            (C, pm(n, s)),
            # ERR Stroud's book wrongly states (t, t,..., t)_FS instead of
            # (t, t, 0, ..., 0)_FS.
            (D, fsd(n, t, 2)),
            ]

        self.points, self.weights = untangle(data)
        self.weights *= integrate_monomial_over_unit_nsphere(n * [0])
        return

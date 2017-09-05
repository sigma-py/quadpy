# -*- coding: utf-8 -*-
#
from __future__ import division

import math

import numpy

from ..helpers import untangle, fsd, pm
from .. import nsphere


class Stroud1967b(object):
    '''
    A.H. Stroud,
    Some Seventh Degree Integration Formulas for Symmetric Regions,
    SIAM J. Numer. Anal., 4(1), 37–44. (8 pages),
    <https://doi.org/10.1137/0704004>.
    '''
    # pylint: disable=too-many-locals
    def __init__(self, index, n):
        self.dim = n
        self.degree = 7

        if index in ['2a', '2b']:
            # the points/weights are complex for n >= 9; one could permit that
            if index == '2a':
                assert n in [2, 3, 4, 6, 7]
                p_m = +1
            else:
                assert index == '2b'
                assert n in [3, 4]
                p_m = -1

            p_m = +1 if index == '2a' else -1

            sqrt38n = math.sqrt(3*(8-n))

            r2 = (3 * (8-n) - p_m * (n-2) * sqrt38n) / 2 / (5-n)
            s2 = (3*n - p_m * 2 * sqrt38n) / 2 / (3*n - 8)
            t2 = (6 + p_m * sqrt38n) / 2
            B = (8 - n) / 8 / r2**3
            C = 1 / 2**(n+3) / s2**3
            D = 1 / 16 / t2**3
            A = 1 - 2*n*B - 2**n*C - 2*n*(n-1)*D

            r = math.sqrt(r2)
            s = math.sqrt(s2)
            t = math.sqrt(t2)

            data = [
                (A, numpy.full((1, n), 0.0)),
                (B, fsd(n, (r, 1))),
                (C, pm(n, s)),
                (D, fsd(n, (t, 2))),
                ]

            self.points, self.weights = untangle(data)
            self.weights *= math.sqrt(math.pi)**n
        else:
            assert index == '4'
            assert n >= 3
            p_m = numpy.array([+1, -1])

            sqrt2n2 = math.sqrt(2*(n+2))
            r1, r2 = numpy.sqrt((n + 2 - p_m * sqrt2n2) / 2)
            g = math.gamma(n / 2)
            A1, A2 = (n + 2 + p_m * sqrt2n2) / 4 / (n+2) * g

            s = nsphere.Stroud1967(n)

            self.points = numpy.concatenate([
                r1 * s.points,
                r2 * s.points,
                ])

            self.weights = numpy.concatenate([
                A1 * s.weights,
                A2 * s.weights,
                ])

        return

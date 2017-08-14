# -*- coding: utf-8 -*-
#
from __future__ import division
import math
import numpy

from ..helpers import untangle, fsd, pm, z
from .helpers import volume_unit_ball

from .. import nsphere


class Stroud1967b(object):
    '''
    A.H. Stroud,
    Some Seventh Degree Integration Formulas for Symmetric Regions,
    SIAM J. Numer. Anal., 4(1), 37–44. (8 pages),
    <https://doi.org/10.1137/0704004>.
    '''
    def __init__(self, n, variant):
        self.degree = 7
        self.dim = n

        if variant in ['a', 'b']:
            if variant == 'a':
                assert 3 <= n <= 7
                t = 1
            else:
                # Stroud mentions nothing of variant b being only valid until
                # dimension 6, but that's the way it is.
                assert variant == 'b'
                assert 3 <= n <= 6
                t = -1

            alpha = math.sqrt(6*(n+6)*(8-n))

            r2 = (3*(n+6)*(8-n) - t*(n-2) * alpha) / ((n+6) * (34-5*n))
            s2 = (3*n*(n+6) - t*2*alpha) / ((n+6) * (3*n**2+6*n-16))
            t2 = (6*(n+6) + t*alpha) / (14*(n+6))

            B = (8-n) / (n+2) / (n+4) / (n+6) / r2**3
            C = 1 / (n+2) / (n+4) / (n+6) / s2**3 / 2**n
            D = 1 / (n+2) / (n+4) / (n+6) / t2**3 / 2
            A = 1 - 2*n*B - 2**n*C - 2*n*(n-1)*D

            r = math.sqrt(r2)
            s = math.sqrt(s2)
            t = math.sqrt(t2)

            data = [
                (A, z(n)),
                (B, fsd(n, r, 1)),
                (C, pm(n, s)),
                (D, fsd(n, t, 2)),
                ]
            self.points, self.weights = untangle(data)
        else:
            assert variant == 'c'
            assert n >= 3

            alpha = math.sqrt(2*(n+2)*(n+4))

            def r(n, t):
                return math.sqrt(((n+2)*(n+4) + t*2*alpha) / (n+4) / (n+6))

            def A(n, t):
                return (2*(n+2)**2 + t*(n-2)*alpha) / (4*n*(n+2)**2)

            r1 = r(n, -1)
            r2 = r(n, +1)

            A1 = A(n, -1)
            A2 = A(n, +1)

            s = nsphere.Stroud1967(n)
            u = s.points
            B = s.weights / math.fsum(s.weights) * n

            self.points = numpy.concatenate([
                r1 * u,
                r2 * u,
                ])

            self.weights = numpy.concatenate([
                A1 * B,
                A2 * B,
                ])

        self.weights *= volume_unit_ball(n)
        return

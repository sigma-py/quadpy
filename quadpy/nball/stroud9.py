# -*- coding: utf-8 -*-
#
from __future__ import division
import math

from ..helpers import untangle, fsd, pm, z
from .helpers import volume_unit_ball


class Stroud9(object):
    # TODO look up citation
    '''
    '''
    def __init__(self, n, variant):
        self.degree = 7
        self.dim = n

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
        self.weights *= volume_unit_ball(n)
        return

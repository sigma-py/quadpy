# -*- coding: utf-8 -*-
#
from __future__ import division
import math
import numpy

from ..helpers import untangle, pm, fsd2
from .helpers import integrate_monomial_over_unit_nsphere


class Stroud1969(object):
    # The reference doesn't really seem to fit here, but this is what Stroud's
    # book says. Perhaps a typo?
    '''
    A.H. Stroud,
    A Fifth Degree Integration Formula for the n-Simplex,
    SIAM J. Numer. Anal., 6(1), 90–98. (9 pages),
    <https://doi.org/10.1137/0706009>.
    '''
    # pylint: disable=too-many-locals
    def __init__(self, n):
        assert 5 <= n <= 16

        self.dim = n
        self.degree = 11

        plus_minus = numpy.array([+1, -1])
        sqrt3 = math.sqrt(3.0)

        t = math.sqrt(1.0 / n)
        r1, r2 = numpy.sqrt(
                (n + 6 - plus_minus*4*sqrt3) / (n**2 + 12*n - 12)
                )
        s1, s2 = numpy.sqrt(
                (7*n - 6 + plus_minus*4*(n-1)*sqrt3) / (n**2 + 12*n - 12)
                )
        u1, u2 = numpy.sqrt(
                (n + 12 + plus_minus*8*sqrt3) / (n**2 + 24*n - 48)
                )
        v1, v2 = numpy.sqrt(
                (7*n - 12 - plus_minus*4*(n-2)*sqrt3) / (n**2 + 24*n - 48)
                )

        A = [
            0.211979378646e-1,
            0.281250000000,
            0.111934731935e+1,
            0.282751322751e+1,
            0.568266145619e+1,
            0.993785824515e+1,
            0.158196616478e+2,
            0.235285714285e+2,
            0.332409299392e+2,
            0.451113811729e+2,
            0.592754264177e+2,
            0.758518518518e+2,
            ]
        B1 = [
            0.638253880175e-1,
            0.452340041459e-1,
            0.337329118818e-1,
            0.261275095270e-1,
            0.208331595340e-1,
            0.169937111647e-1,
            0.141147212492e-1,
            0.118949128383e-1,
            0.101424250926e-1,
            0.873046796644e-2,
            0.757257014768e-2,
            0.660813369775e-2,
            ]
        B2 = [
            +0.213579471658e-1,
            -0.108726067638,
            -0.371589499738,
            -0.786048144448,
            +0.136034060198e+1,
            +0.209547695631e+1,
            +0.298784764467e+1,
            +0.403107480702e+1,
            +0.521726499521e+1,
            +0.653783099707e+1,
            +0.798401677102e+1,
            +0.954722261180e+1,
            ]
        C1 = [
            0.236639091329e-1,
            0.525940190875e-1,
            0.925052768546e-1,
            0.141316953438,
            0.196818580052,
            0.257027634179,
            0.320299222258,
            0.385326226441,
            0.451098131789,
            0.516849445559,
            0.582010515746,
            0.646165210110,
            ]
        C2 = [
            0.316246294890e-1,
            0.207194729760e-1,
            0.144303800811e-1,
            0.105348984135e-1,
            0.798435122193e-2,
            0.623845929545e-2,
            0.499896882962e-2,
            0.409176297655e-2,
            0.341037426698e-2,
            0.288710646943e-2,
            0.247745182907e-2,
            0.215128820597e-2,
            ]

        data = [
            (A[n-5], pm(n, t)),
            (B1[n-5], fsd2(n, s1, r1, 1, n-1)),
            (B2[n-5], fsd2(n, s2, r2, 1, n-1)),
            (C1[n-5], fsd2(n, v1, u1, 2, n-2)),
            (C2[n-5], fsd2(n, v2, u2, 2, n-2)),
            ]

        self.points, self.weights = untangle(data)
        self.weights *= integrate_monomial_over_unit_nsphere(n * [0])
        self.weights /= 2.0**n

        return

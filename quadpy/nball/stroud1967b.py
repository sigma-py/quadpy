# -*- coding: utf-8 -*-
#
from __future__ import division

import numpy
import sympy

from ..helpers import untangle, fsd, pm, z
from .helpers import volume_unit_ball

from .. import nsphere


class Stroud1967b(object):
    """
    A.H. Stroud,
    Some Seventh Degree Integration Formulas for Symmetric Regions,
    SIAM J. Numer. Anal., 4(1), 37–44. (8 pages),
    <https://doi.org/10.1137/0704004>.
    """

    def __init__(self, n, variant, symbolic=False):
        sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt

        self.degree = 7
        self.dim = n

        if variant in ["a", "b"]:
            if variant == "a":
                assert 3 <= n <= 7
                t = 1
            else:
                # Stroud mentions nothing of variant b being only valid until
                # dimension 6, but that's the way it is.
                assert variant == "b"
                assert 3 <= n <= 6
                t = -1

            alpha = sqrt(6 * (n + 6) * (8 - n))

            r2 = (3 * (n + 6) * (8 - n) - t * (n - 2) * alpha) / (
                (n + 6) * (34 - 5 * n)
            )
            s2 = (3 * n * (n + 6) - t * 2 * alpha) / (
                (n + 6) * (3 * n ** 2 + 6 * n - 16)
            )
            t2 = (6 * (n + 6) + t * alpha) / (14 * (n + 6))

            B = (8 - n) / (n + 2) / (n + 4) / (n + 6) / r2 ** 3
            C = 1 / (n + 2) / (n + 4) / (n + 6) / s2 ** 3 / 2 ** n
            D = 1 / (n + 2) / (n + 4) / (n + 6) / t2 ** 3 / 2
            A = 1 - 2 * n * B - 2 ** n * C - 2 * n * (n - 1) * D

            r = sqrt(r2)
            s = sqrt(s2)
            t = sqrt(t2)

            data = [(A, z(n)), (B, fsd(n, (r, 1))), (C, pm(n, s)), (D, fsd(n, (t, 2)))]
            self.points, self.weights = untangle(data)

            self.weights *= volume_unit_ball(n, symbolic=symbolic)
        else:
            assert variant == "c"
            assert n >= 3

            alpha = sqrt(2 * (n + 2) * (n + 4))

            pm_ = numpy.array([+1, -1])
            r1, r2 = sqrt(((n + 2) * (n + 4) + pm_ * 2 * alpha) / (n + 4) / (n + 6))
            A1, A2 = (2 * (n + 2) ** 2 + pm_ * (n - 2) * alpha) / (4 * n * (n + 2) ** 2)

            s = nsphere.Stroud1967(n, symbolic=symbolic)

            self.points = numpy.concatenate([r1 * s.points, r2 * s.points])

            self.weights = numpy.concatenate([A1 * s.weights, A2 * s.weights])
        return

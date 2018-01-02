# -*- coding: utf-8 -*-
#
from __future__ import division

import numpy
import sympy

from ..helpers import untangle
from .helpers import volume_unit_ball


class Stroud1957(object):
    '''
    A. H. Stroud,
    Remarks on the Disposition of Points in Numerical Integration Formulas,
    Mathematical Tables and Other Aids to Computation,
    Vol. 11, No. 60 (Oct., 1957), pp. 257-261,
    <https://doi.org/10.2307/2001945>.
    '''
    def __init__(self, n, symbolic=False):
        frac = sympy.Rational if symbolic else lambda x, y: x/y
        pi = sympy.pi if symbolic else numpy.pi
        sqrt = sympy.sqrt if symbolic else numpy.sqrt
        cos = sympy.cos if symbolic else numpy.cos
        sin = sympy.sin if symbolic else numpy.sin
        pi = sympy.pi if symbolic else numpy.pi

        self.degree = 2
        self.dim = n

        pts = [[
            [sqrt(frac(2, n+2)) * cos(2*k*i*pi / (n+1)) for i in range(n+1)],
            [sqrt(frac(2, n+2)) * sin(2*k*i*pi / (n+1)) for i in range(n+1)],
            ] for k in range(1, n//2 + 1)]
        if n % 2 == 1:
            sqrt3pm = numpy.full(n+1, 1/sqrt(n+2))
            sqrt3pm[1::2] *= -1
            pts.append(sqrt3pm)
        pts = numpy.vstack(pts).T

        data = [(frac(1, n+1), pts)]

        self.points, self.weights = untangle(data)

        self.weights *= volume_unit_ball(n, symbolic=symbolic)
        return

# -*- coding: utf-8 -*-
#
from __future__ import division

import numpy
import sympy

from .helpers import _s

from ..helpers import untangle


class Stroud1957(object):
    '''
    A. H. Stroud,
    Remarks on the Disposition of Points in Numerical Integration Formulas,
    Mathematical Tables and Other Aids to Computation,
    Vol. 11, No. 60 (Oct., 1957), pp. 257-261,
    <https://doi.org/10.2307/2001945>.
    '''
    # pylint: disable=too-many-locals
    def __init__(self, n, index, symbolic=True):
        frac = sympy.Rational if symbolic else lambda x, y: x/y
        sqrt = sympy.sqrt if symbolic else numpy.sqrt
        pi = sympy.pi if symbolic else numpy.pi
        sin = sympy.sin if symbolic else numpy.sin
        cos = sympy.cos if symbolic else numpy.cos

        self.dim = n
        if index == 2:
            self.degree = 2
            r = sqrt(3) / 6
            data = [
                (1.0, numpy.array([numpy.full(n, 2*r)])),
                (+r, _s(n, -1, r)),
                (-r, _s(n, +1, r)),
                ]
        else:
            assert index == 3
            self.degree = 3
            n2 = n // 2 if n % 2 == 0 else (n-1)//2
            i_range = range(1, 2*n+1)
            pts = [[
                [sqrt(frac(2, 3)) * cos((2*k-1)*i*pi / n) for i in i_range],
                [sqrt(frac(2, 3)) * sin((2*k-1)*i*pi / n) for i in i_range],
                ] for k in range(1, n2+1)]
            if n % 2 == 1:
                sqrt3pm = numpy.full(2*n, 1/sqrt(3))
                sqrt3pm[1::2] *= -1
                pts.append(sqrt3pm)
            pts = numpy.vstack(pts).T

            data = [(frac(1, 2*n), pts)]

        self.points, self.weights = untangle(data)
        reference_volume = 2**n
        self.weights *= reference_volume
        return

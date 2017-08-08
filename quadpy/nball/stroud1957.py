# -*- coding: utf-8 -*-
#
import numpy

from ..helpers import untangle
from .helpers import volume_unit_ball


class Stroud1957(object):
    '''
    A. H. Stroud,
    Remarks on the Disposition of Points in Numerical Integration Formulas,
    Mathematical Tables and Other Aids to Computation,
    Vol. 11, No. 60 (Oct., 1957), pp. 257-261,
    <https://dx.doi.org/10.2307/2001945>.
    '''
    def __init__(self, n):
        self.degree = 2
        self.dim = n

        i = numpy.arange(n+1)
        n2 = n // 2 if n % 2 == 0 else (n-1)//2
        pts = [[
            numpy.sqrt(2.0/(n+2)) * numpy.cos(2*k*i*numpy.pi / (n+1)),
            numpy.sqrt(2.0/(n+2)) * numpy.sin(2*k*i*numpy.pi / (n+1)),
            ] for k in range(1, n2+1)]
        if n % 2 == 1:
            sqrt3pm = numpy.full(n+1, 1.0 / numpy.sqrt(n+2))
            sqrt3pm[1::2] *= -1
            pts.append(sqrt3pm)
        pts = numpy.vstack(pts).T

        data = [(1.0/(n+1), pts)]

        self.points, self.weights = untangle(data)
        self.weights *= volume_unit_ball(n)
        return

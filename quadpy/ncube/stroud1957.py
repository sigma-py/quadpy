# -*- coding: utf-8 -*-
#
import numpy

from .helpers import _s

from ..helpers import untangle


class Stroud1957(object):
    '''
    A. H. Stroud,
    Remarks on the Disposition of Points in Numerical Integration Formulas,
    Mathematical Tables and Other Aids to Computation,
    Vol. 11, No. 60 (Oct., 1957), pp. 257-261,
    <https://dx.doi.org/10.2307/2001945>.
    '''
    def __init__(self, n, index):
        reference_volume = 2.0**n
        self.dim = n
        if index == 2:
            self.degree = 2
            r = numpy.sqrt(3.0) / 6.0
            data = [
                (1.0, numpy.array([numpy.full(n, 2*r)])),
                (+r, _s(n, -1.0, r)),
                (-r, _s(n, +1.0, r)),
                ]
        else:
            assert index == 3
            self.degree = 3
            i = numpy.arange(1, 2*n+1)
            n2 = n // 2 if n % 2 == 0 else (n-1)//2
            pts = [[
                numpy.sqrt(2.0/3.0) * numpy.cos((2*k-1)*i*numpy.pi / n),
                numpy.sqrt(2.0/3.0) * numpy.sin((2*k-1)*i*numpy.pi / n),
                ] for k in range(1, n2+1)]
            if n % 2 == 1:
                sqrt3pm = numpy.full(2*n, 1.0 / numpy.sqrt(3.0))
                sqrt3pm[1::2] *= -1
                pts.append(sqrt3pm)
            pts = numpy.vstack(pts).T

            data = [(1.0/(2*n), pts)]

        self.points, self.weights = untangle(data)
        self.weights *= reference_volume
        return

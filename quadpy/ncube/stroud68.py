# -*- coding: utf-8 -*-
#
import numpy

from .helpers import _s2, _s11, _z, _fsd


class Stroud68(object):
    '''
    A. H. Stroud,
    Extensions of Symmetric Integration Formulas,
    Mathematics of Computation,
    Vol. 22, No. 102 (Apr., 1968), pp. 271-274,
    Published by: American Mathematical Society,
    <https://dx.doi.org/10.2307/2004655>.
    '''
    def __init__(self, n):
        reference_volume = 2.0**n
        self.degree = 5
        r = numpy.sqrt(7.0 / 15.0)
        s = numpy.sqrt((7.0 + numpy.sqrt(24.0)) / 15.0)
        t = numpy.sqrt((7.0 - numpy.sqrt(24.0)) / 15.0)
        self.points = numpy.concatenate([
            _z(n),
            _s2(n, r),
            _s2(n, -r),
            _fsd(n, r, 1),
            _s11(n, +s, -t),
            _s11(n, -s, +t),
            _fsd(n, s, 1),
            _fsd(n, t, 1),
            ])
        self.weights = numpy.concatenate([
            numpy.full(1, (5*n**2 - 15*n+14)/14.0 * reference_volume),
            numpy.full(n*(n-1), 25.0/168.0 * reference_volume),
            numpy.full(2*n, -25*(n-2)/168.0 * reference_volume),
            numpy.full(2*n*(n-1), 5.0/48.0 * reference_volume),
            numpy.full(4*n, -5*(n-2)/48.0 * reference_volume),
            ])
        return

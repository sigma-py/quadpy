# -*- coding: utf-8 -*-
#
import numpy

from .helpers import _s2, _s11, _z, _fsd

from ..helpers import untangle


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
        data = [
            ((5*n**2 - 15*n+14)/14.0, _z(n)),
            (25.0/168.0, _s2(n, +r)),
            (25.0/168.0, _s2(n, -r)),
            (-25*(n-2)/168.0, _fsd(n, r, 1)),
            (5.0/48.0, _s11(n, +s, -t)),
            (5.0/48.0, _s11(n, -s, +t)),
            (-5*(n-2)/48.0, _fsd(n, s, 1)),
            (-5*(n-2)/48.0, _fsd(n, t, 1)),
            ]

        self.points, self.weights = untangle(data)
        self.weights *= reference_volume
        return

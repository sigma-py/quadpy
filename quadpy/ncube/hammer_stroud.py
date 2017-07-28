# -*- coding: utf-8 -*-
#
import numpy

from .helpers import _z, _fsd


class HammerStroud(object):
    '''
    Preston C. Hammer and Arthur H. Stroud,
    Numerical Evaluation of Multiple Integrals II,
    Math. Comp. 12 (1958), 272-280,
    <https://doi.org/10.1090/S0025-5718-1958-0102176-6>.
    '''
    def __init__(self, n):
        reference_volume = 2.0**n
        self.degree = 5
        r = numpy.sqrt(3.0 / 5.0)
        self.points = numpy.concatenate([
            _z(n),
            _fsd(n, r, 1),
            _fsd(n, r, 2),
            ])
        self.weights = numpy.concatenate([
            numpy.full(
                1, (25*n**2 - 115*n + 162)/162.0 * reference_volume
                ),
            numpy.full(2*n, (70 - 25*n)/162.0 * reference_volume),
            numpy.full(2*n*(n-1), 25.0/324.0 * reference_volume),
            ])
        return

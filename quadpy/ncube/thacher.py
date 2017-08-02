# -*- coding: utf-8 -*-
#
import numpy

from .helpers import _s


class Thacher(object):
    '''
    Henry C. Thacher,
    An efficient composite formula for multidimensional quadrature,
    Communications of the ACM CACM,
    Volume 7 Issue 1, Jan. 1964, Pages 23-25,
    ACM New York, NY, USA,
    <https://dx.doi.org/10.1145/363872.363897>.

    Abstract:
    A (2s+1)-point, second-degree quadrature formula for integration over an
    s-dimensional hyper-rectangle is presented. All but one of the points lie
    on the surface with weights of opposite sign attached to points on opposite
    faces. When a large volume is subdivided into congruent rectangular
    subdivisions, only one point is required in each interior subdivision to
    achieve second-degree accuracy.
    '''
    def __init__(self, n):
        reference_volume = 2.0**n
        self.degree = 2
        r = numpy.sqrt(3.0) / 6.0
        self.weights = numpy.concatenate([
            numpy.full(1, reference_volume),
            numpy.full(n, r*reference_volume),
            numpy.full(n, -r*reference_volume),
            ])
        self.points = numpy.concatenate([
            numpy.array([numpy.full(n, 2*r)]),
            _s(n, -1.0, r),
            _s(n, +1.0, r),
            ])
        return

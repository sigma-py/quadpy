# -*- coding: utf-8 -*-
#
import numpy

from ..helpers import untangle, z


class Albrecht(object):
    # Look up citation
    '''
    Albrecht [1]
    '''
    def __init__(self):
        # Stroud claims degree 3, but really it's only order 1.
        self.degree = 9

        rho1 = numpy.sqrt((96.0 - 4*numpy.sqrt(111.0)) / 155.0)
        rho2 = numpy.sqrt((96.0 + 4*numpy.sqrt(111.0)) / 155.0)

        k = numpy.arange(1, 7)
        x = numpy.array([
            numpy.cos(k*numpy.pi/3.0),
            numpy.sin(k*numpy.pi/3.0),
            ]).T

        uv = numpy.sqrt(0.8) * numpy.array([
            numpy.cos((2*k-1)*numpy.pi / 6.0),
            numpy.sin((2*k-1)*numpy.pi / 6.0),
            ]).T

        B0 = 251.0 / 2304.0
        B1 = (110297.0 + 5713.0*numpy.sqrt(111.0)) / 2045952.0
        B2 = (110297.0 - 5713.0*numpy.sqrt(111.0)) / 2045952.0
        C = 125.0 / 3072.0

        data = [
            (B0, z(2)),
            (B1, rho1 * x),
            (B2, rho2 * x),
            (C, uv),
            ]
        self.points, self.weights = untangle(data)
        self.weights *= numpy.pi
        return

# -*- coding: utf-8 -*-
#
from math import pi
import numpy

from .. import helpers


def show(
        scheme,
        backend='mpl'
        ):
    '''Displays scheme for 3D ball quadrature.
    '''
    helpers.backend_to_function[backend](
            scheme.points,
            scheme.weights,
            volume=4.0/3.0*pi,
            edges=[],
            balls=[((0.0, 0.0, 0.0), 1.0)],
            )
    return


def integrate(f, center, radius, rule, dot=numpy.dot):
    flt = numpy.vectorize(float)

    center = numpy.array(center)
    rr = numpy.multiply.outer(radius, flt(rule.points))
    rr = numpy.swapaxes(rr, 0, -2)
    ff = numpy.array(f((rr + center).T))

    return numpy.array(radius)**3 * dot(ff, flt(rule.weights))

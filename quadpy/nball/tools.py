# -*- coding: utf-8 -*-
#
import numpy

from .. import helpers


def integrate(f, center, radius, rule, sumfun=helpers.kahan_sum):
    center = numpy.array(center)
    rr = numpy.multiply.outer(radius, rule.points)
    rr = numpy.swapaxes(rr, 0, -2)
    ff = numpy.array(f((rr + center).T))
    out = sumfun(rule.weights * ff, axis=-1)
    return numpy.array(radius)**rule.dim * out

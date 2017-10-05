# -*- coding: utf-8 -*-
#
import numpy

from .. import helpers


def integrate(f, center, radius, rule, sumfun=helpers.kahan_sum):
    flt = numpy.vectorize(float)

    center = numpy.array(center)
    rr = numpy.multiply.outer(radius, flt(rule.points))
    rr = numpy.swapaxes(rr, 0, -2)
    ff = numpy.array(f((rr + center).T))
    out = sumfun(flt(rule.weights) * ff, axis=-1)

    return numpy.array(radius)**(rule.dim-1) * out

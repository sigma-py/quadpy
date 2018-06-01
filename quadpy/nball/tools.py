# -*- coding: utf-8 -*-
#
import numpy


def integrate(f, center, radius, rule, dot=numpy.dot):
    center = numpy.array(center)
    rr = numpy.multiply.outer(radius, rule.points)
    rr = numpy.swapaxes(rr, 0, -2)
    ff = numpy.array(f((rr + center).T))
    return numpy.array(radius) ** rule.dim * dot(ff, rule.weights)

# -*- coding: utf-8 -*-
#
import numpy


def integrate(f, rule, dot=numpy.dot):
    flt = numpy.vectorize(float)
    return dot(f(flt(rule.points).T), flt(rule.weights))

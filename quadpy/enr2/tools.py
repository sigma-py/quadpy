# -*- coding: utf-8 -*-
#
import numpy

from .. import helpers


def integrate(f, rule, sumfun=helpers.kahan_sum):
    flt = numpy.vectorize(float)
    ff = numpy.array(f(flt(rule.points).T))
    return sumfun(flt(rule.weights) * ff, axis=-1)

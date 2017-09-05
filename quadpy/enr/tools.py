# -*- coding: utf-8 -*-
#
import numpy

from .. import helpers


def integrate(f, rule, sumfun=helpers.kahan_sum):
    ff = numpy.array(f(rule.points.T))
    return sumfun(rule.weights * ff, axis=-1)

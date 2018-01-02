# -*- coding: utf-8 -*-
#
import numpy


def integrate(f, rule, dot=numpy.dot):
    return dot(f(rule.points.T), rule.weights)

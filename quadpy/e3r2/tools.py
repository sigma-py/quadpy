# -*- coding: utf-8 -*-
#
import numpy

from .. import helpers


def integrate(f, rule, dot=numpy.dot):
    flt = numpy.vectorize(float)
    return dot(f(flt(rule.points).T), flt(rule.weights))


def show(scheme, backend="mpl"):
    """Displays scheme for E_3^{r^2} quadrature.
    """
    helpers.backend_to_function[backend](
        scheme.points, scheme.weights, volume=numpy.pi ** 1.5, edges=[]
    )
    return

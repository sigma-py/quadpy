# -*- coding: utf-8 -*-
#
import numpy

from .. import helpers


def integrate(f, rule, sumfun=helpers.kahan_sum):
    ff = numpy.array(f(rule.points.T))
    return sumfun(rule.weights * ff, axis=-1)


def show(
        scheme,
        backend='mpl'
        ):
    '''Displays scheme for E_3^r quadrature.
    '''
    helpers.backend_to_function[backend](
            scheme.points,
            scheme.weights,
            volume=8*numpy.pi,
            edges=[]
            )
    return

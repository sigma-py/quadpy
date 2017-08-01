# -*- coding: utf-8 -*-
#
import numpy


def transform(xi, simplex):
    '''Transform the points `xi` from the reference simplex onto `simplex`.
    '''
    # For n == 2:
    # x = (
    #     + outer(triangle[0].T, 1.0 - xi[0] - xi[1])
    #     + outer(triangle[1].T, xi[0])
    #     + outer(triangle[2].T, xi[1])
    #     )
    shape_funs = numpy.vstack([1.0 - numpy.sum(xi, axis=0), xi])
    return numpy.dot(simplex, shape_funs)


def get_vol(simplex):
    n = simplex.shape[1]
    return abs(get_detJ(simplex)) / numpy.math.factorial(n)


def get_detJ(simplex):
    # det is the signed volume of the tetrahedron
    # <https://en.wikipedia.org/wiki/Simplex#Volume>
    J = simplex[1:] - simplex[0]
    return numpy.linalg.det(numpy.moveaxis(J, 0, 1))

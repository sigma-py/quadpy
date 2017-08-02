# -*- coding: utf-8 -*-
#
import math
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
    # Compute the volume via the Cayley-Menger determinant
    # <http://mathworld.wolfram.com/Cayley-MengerDeterminant.html>. One
    # advantage is that it can compute the volume of the simplex indenpendent
    # of the dimension of the space where it's embedded.

    # compute all edge lengths
    edges = numpy.subtract(simplex[:, None], simplex[None, :])
    ei_dot_ej = numpy.einsum('...k,...k->...', edges, edges)

    j = simplex.shape[0] - 1
    a = numpy.empty((j+2, j+2) + ei_dot_ej.shape[2:])
    a[1:, 1:] = ei_dot_ej
    a[0, 1:] = 1.0
    a[1:, 0] = 1.0
    a[0, 0] = 0.0

    a = numpy.moveaxis(a, (0, 1), (-2, -1))
    det = numpy.linalg.det(a)

    vol = numpy.sqrt((-1.0)**(j+1) / 2**j / math.factorial(j)**2 * det)
    return vol

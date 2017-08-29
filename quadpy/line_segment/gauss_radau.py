# -*- coding: utf-8 -*-
#
import numpy
from orthopy import jacobi_recurrence_coefficients, gauss_from_coefficients


class GaussRadau(object):
    '''
    Gauss-Radau quadrature.
    '''
    def __init__(self, n, a=0.0, b=0.0):
        assert n >= 2
        self.degree = 2*n - 1
        alpha, beta = jacobi_recurrence_coefficients(n, a, b)
        self.points, self.weights = _radau(alpha, beta, -1.0)
        return


def _radau(alpha, beta, xr):
    '''From <http://www.scientificpython.net/pyblog/radau-quadrature>:
    Compute the Radau nodes and weights with the preassigned node xr.

    Based on the section 7 of the paper

        Some modified matrix eigenvalue problems,
        Gene Golub,
        SIAM Review Vol 15, No. 2, April 1973, pp.318--334.
    '''
    from scipy.linalg import solve_banded

    n = len(alpha)-1
    f = numpy.zeros(n)
    f[-1] = beta[-1]
    A = numpy.vstack((numpy.sqrt(beta), alpha-xr))
    J = numpy.vstack((A[:, 0:-1], A[0, 1:]))
    delta = solve_banded((1, 1), J, f)
    alphar = alpha.copy()
    alphar[-1] = xr + delta[-1]
    x, w = gauss_from_coefficients(alphar, beta)
    return x, w

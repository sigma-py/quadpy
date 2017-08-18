# -*- coding: utf-8 -*-
#
import numpy

from .helpers import _gauss


class Gauss(object):
    '''Given moments

    mu_k = int_a^b omega(x) x^k dx,  k = {0, 1,...,2N}

    (with omega being a nonnegative weight function), this class creates the
    Gauss scheme corresponding to the above integral. It uses the mechanism
    from

    Gene H. Golub and John H. Welsch,
    Calculation of Gauss Quadrature Rules,
    Mathematics of Computation,
    Vol. 23, No. 106 (Apr., 1969), pp. 221-230+s1-s10,
    <https://dx.doi.org/10.2307/2004418>,
    <https://pdfs.semanticscholar.org/c715/119d5464f614fd8ec590b732ccfea53e72c4.pdf>.
    '''
    def __init__(self, n, moments):
        self.degree = 2*n-1

        M = numpy.array([[
            moments[i+j] for j in range(n+1)
            ] for i in range(n+1)])
        R = numpy.linalg.cholesky(M).T

        # (upper) diagonal
        Rd = R.diagonal()
        q = R.diagonal(1) / Rd[:-1]

        alpha = numpy.zeros(n)
        alpha = q.copy()
        alpha[+1:] -= q[:-1]

        # TODO don't square here, but adapt _gauss to accept squared values at
        #      input
        beta = numpy.hstack([
            Rd[0], Rd[1:-1] / Rd[:-2]
            ])**2

        self.points, self.weights = _gauss(alpha, beta)
        return

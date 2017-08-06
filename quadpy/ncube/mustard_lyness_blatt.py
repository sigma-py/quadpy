# -*- coding: utf-8 -*-
#
import numpy

from .helpers import _pm

from ..helpers import untangle, fsd, z


class MustardLynessBlatt(object):
    '''
    D. Mustard, J.N. Lyness, J.M. Blatt,
    Numerical quadrature in n dimensions,
    Comput J (1963) 6 (1): 75-87,
    <https://doi.org/10.1093/comjnl/6.1.75>.

    Abstract:
    We investigate a selection of integration rules based on the combination of
    third degree rules for elementary three-dimensional sub-domains. The
    practical problems associated with the application of these rules to
    domains bounded by planes of a certain type are discussed in detail. These
    included integration rules for less symmetrical domains which may occur
    near the boundary of the volume of integration, methods for combining
    integration coefficients from adjacent sub-domains, and methods for
    changing net size within the volume of integration.
    Particular attention is paid to minimizing the number of points at which
    the function has to be evaluated, and error estimates in terms of
    computation time are given. A list of integration coefficients of general
    interest for three-dimensional integrations is presented. The discussion is
    generalized to n dimensions for hyper-cubic domains.
    '''
    def __init__(self, n):
        reference_volume = 2.0**n
        self.degree = 5
        r = numpy.sqrt(2.0 / 5.0)
        data = [
            ((8 - 5*n)/9.0, z(n)),
            (5.0/18.0, fsd(n, r, 1)),
            (1.0/9.0 / 2**n, _pm(n, 1.0)),
            ]

        self.points, self.weights = untangle(data)
        self.weights *= reference_volume
        return

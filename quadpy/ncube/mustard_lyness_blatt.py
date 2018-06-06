# -*- coding: utf-8 -*-
#
from __future__ import division

import numpy
import sympy

from ..helpers import untangle, fsd, z, pm


class MustardLynessBlatt(object):
    """
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
    """

    def __init__(self, n, symbolic=False):
        frac = sympy.Rational if symbolic else lambda x, y: x / y
        sqrt = sympy.sqrt if symbolic else numpy.sqrt

        self.degree = 5
        r = sqrt(frac(2, 5))
        data = [
            (frac(8 - 5 * n, 9), z(n)),
            (frac(5, 18), fsd(n, (r, 1))),
            (frac(1, 9 * 2 ** n), pm(n, 1)),
        ]

        self.points, self.weights = untangle(data)
        reference_volume = 2 ** n
        self.weights *= reference_volume
        return

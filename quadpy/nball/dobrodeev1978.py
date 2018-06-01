# -*- coding: utf-8 -*-
#
from __future__ import division

from ..helpers import untangle, fsd, compute_dobrodeev
from .helpers import integrate_monomial_over_unit_nball


class Dobrodeev1978(object):
    '''
    L.N. Dobrodeev,
    Cubature rules with equal coefficients for integrating functions with
    respect to symmetric domains,
    USSR Computational Mathematics and Mathematical Physics,
    Volume 18, Issue 4, 1978, Pages 27-34,
    <https://doi.org/10.1016/0041-5553(78)90064-2>.

    Abstract:
    Three-parameter cubature sums, whose values have to be the same as the
    values of the computed integral on all polynomials of degree p <= 5, are
    investigated. The parameters are selected on the basis of a specific
    optimality criterion for a number of symmetric domains of integration, of
    dimensionalities 2 <= n <= 20.
    '''
    def __init__(self, n, symbolic=False):
        assert 2 <= n <= 20

        self.name = 'Dobrodeev1978'
        self.degree = 5
        self.dim = n

        dim_config = {
            2: ('II', 1, 1, 1),
            3: ('II', 1, 1, 1),
            4: ('II', 4, 1, 1),
            5: ('I', 3, 1, 1),
            6: ('I', 3, 1, 1),
            7: ('I', 7, 1, 1),
            8: ('I', 8, 1, 1),
            9: ('I', 9, 1, 1),
            10: ('I', 9, 1, 2),
            11: ('I', 5, 1, 2),
            12: ('I', 12, 1, 2),
            13: ('I', 13, 1, 2),
            14: ('I', 14, 1, 2),
            15: ('I', 15, 1, 2),
            16: ('II', 15, 1, 4),
            17: ('I', 17, 1, 3),
            18: ('I', 18, 1, 3),
            19: ('I', 19, 1, 3),
            20: ('I', 20, 1, 3),
            }

        I0 = integrate_monomial_over_unit_nball(n * [0], symbolic=symbolic)
        I2 = integrate_monomial_over_unit_nball(
                [2] + (n-1) * [0],
                symbolic=symbolic
                )
        I22 = integrate_monomial_over_unit_nball(
                [2, 2] + (n-2) * [0],
                symbolic=symbolic
                )
        I4 = integrate_monomial_over_unit_nball(
                [4] + (n-1) * [0],
                symbolic=symbolic
                )
        pm_type, i, j, k = dim_config[n]
        G, a, b, c = compute_dobrodeev(
            n, I0, I2, I22, I4, pm_type, i, j, k, symbolic=symbolic
            )

        data = [
            (G, fsd(n, (a, i))),
            (G, fsd(n, (b, j), (c, k))),
            ]

        self.points, self.weights = untangle(data)
        self.weights *= I0
        return

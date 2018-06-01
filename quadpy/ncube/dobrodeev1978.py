# -*- coding: utf-8 -*-
#
from __future__ import division

from ..helpers import untangle, fsd, compute_dobrodeev
from .helpers import integrate_monomial_over_ncube


class Dobrodeev1978(object):
    """
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
    """

    def __init__(self, n, symbolic):
        assert 2 <= n <= 20

        self.name = "Dobrodeev1978"
        self.degree = 5
        self.dim = n

        dim_config = {
            2: ("II", 1, 1, 1),
            3: ("I", 3, 1, 1),
            4: ("I", 4, 1, 1),
            5: ("I", 5, 1, 1),
            6: ("I", 5, 1, 1),
            7: ("I", 6, 1, 2),
            8: ("I", 7, 1, 2),
            9: ("I", 6, 1, 2),
            10: ("I", 6, 1, 9),
            11: ("I", 6, 1, 10),
            12: ("I", 7, 1, 11),
            13: ("I", 8, 2, 2),
            14: ("I", 8, 1, 13),
            15: ("I", 9, 1, 14),
            16: ("I", 10, 2, 3),
            17: ("I", 10, 2, 3),
            18: ("I", 15, 3, 3),
            19: ("II", 9, 1, 18),
            20: ("I", 9, 1, 18),
        }

        ncube_limits = [[-1.0, 1.0]] * n
        I0 = integrate_monomial_over_ncube(ncube_limits, n * [0])
        I2 = integrate_monomial_over_ncube(ncube_limits, [2] + (n - 1) * [0])
        I22 = integrate_monomial_over_ncube(ncube_limits, [2, 2] + (n - 2) * [0])
        I4 = integrate_monomial_over_ncube(ncube_limits, [4] + (n - 1) * [0])

        pm_type, i, j, k = dim_config[n]

        G, a, b, c = compute_dobrodeev(
            n, I0, I2, I22, I4, pm_type, i, j, k, symbolic=symbolic
        )

        data = [(G, fsd(n, (a, i))), (G, fsd(n, (b, j), (c, k)))]

        self.points, self.weights = untangle(data)
        self.weights *= I0
        return

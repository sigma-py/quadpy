# -*- coding: utf-8 -*-
#
import math

import scipy.special
import sympy


def integrate_monomial_over_unit_simplex(k, symbolic=False):
    """The integrals of monomials over the standard triangle and tetrahedron are
    given by

    \\int_T x_0^k0 * x1^k1 = (k0!*k1!) / (2+k0+k1)!,
    \\int_T x_0^k0 * x1^k1 * x2^k2 = (k0!*k1!*k2!) / (4+k0+k1+k2)!,

    see, e.g.,
    A set of symmetric quadrature rules on triangles and tetrahedra,
    Linbo Zhang, Tao Cui and Hui Liu,
    Journal of Computational Mathematics,
    Vol. 27, No. 1 (January 2009), pp. 89-96,
    <https://www.jstor.org/stable/43693493>.

    See, e.g., <https://math.stackexchange.com/q/207073/36678> for a formula in
    all dimensions.
    """
    if symbolic:
        return sympy.prod([sympy.gamma(kk + 1) for kk in k]) / sympy.gamma(
            sum(k) + len(k) + 1
        )
    # exp-log to account for large values in numerator and denominator
    # import scipy.special
    return math.exp(
        math.fsum([scipy.special.gammaln(kk + 1) for kk in k])
        - scipy.special.gammaln(sum([kk + 1 for kk in k]) + 1)
    )

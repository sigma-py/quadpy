# -*- coding: utf-8 -*-
#
from ..nsphere.helpers import integrate_monomial_over_unit_nsphere


def volume_unit_ball(n, symbolic=False):
    return integrate_monomial_over_unit_nball(n * [0], symbolic=symbolic)


def integrate_monomial_over_unit_nball(exponents, symbolic=False):
    '''
    Gerald B. Folland,
    How to Integrate a Polynomial over a Sphere,
    The American Mathematical Monthly,
    Vol. 108, No. 5 (May, 2001), pp. 446-448,
    <https://doi.org/10.2307/2695802>.
    '''
    radius = 1
    n = len(exponents)
    alpha = n + sum(exponents)
    return (
        radius**alpha
        * integrate_monomial_over_unit_nsphere(exponents, symbolic=symbolic)
        / alpha
        )

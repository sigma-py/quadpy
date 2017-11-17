# -*- coding: utf-8 -*-
#
from ..nsphere.helpers import integrate_monomial_over_unit_nsphere


def volume_unit_ball(n):
    return integrate_monomial_over_unit_nball(n * [0])


def integrate_monomial_over_unit_nball(exp):
    '''
    Gerald B. Folland,
    How to Integrate a Polynomial over a Sphere,
    The American Mathematical Monthly,
    Vol. 108, No. 5 (May, 2001), pp. 446-448,
    <https://doi.org/10.2307/2695802>.
    '''
    radius = 1
    n = len(exp)
    alpha = n + sum(exp)
    return radius**alpha * integrate_monomial_over_unit_nsphere(exp) / alpha

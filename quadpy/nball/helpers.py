# -*- coding: utf-8 -*-
#
from ..nsphere.helpers import integrate_monomial_over_unit_nsphere


def integrate_monomial_over_unit_nball(exp):
    '''
    Gerald B. Folland,
    How to Integrate a Polynomial over a Sphere,
    The American Mathematical Monthly,
    Vol. 108, No. 5 (May, 2001), pp. 446-448,
    <https://dx.doi.org/10.2307/2695802>.
    '''
    radius = 1.0
    n = len(exp)
    alpha = n + sum(exp)
    return radius**alpha / alpha * integrate_monomial_over_unit_nsphere(exp)

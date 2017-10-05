# -*- coding: utf-8 -*-
#
from sympy import gamma, prod, Rational as fr


def integrate_monomial_over_unit_nsphere(alpha):
    '''
    Gerald B. Folland,
    How to Integrate a Polynomial over a Sphere,
    The American Mathematical Monthly,
    Vol. 108, No. 5 (May, 2001), pp. 446-448,
    <https://dx.doi.org/10.2307/2695802>.
    '''
    if any(a % 2 == 1 for a in alpha):
        return 0
    # Use lgamma since other with ordinary gamma, numerator and denominator
    # might overflow.
    # return 2 * math.exp(
    #     math.fsum([math.lgamma(0.5*(a+1)) for a in alpha])
    #     - math.lgamma(math.fsum([0.5*(a+1) for a in alpha]))
    #     )
    return 2 * (
        prod([gamma(fr(a+1, 2)) for a in alpha])
        / gamma(sum([fr(a+1, 2) for a in alpha]))
        )

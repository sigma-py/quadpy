from ..ncube._mcnamee_stenger import (
    _mcnamee_stenger_3,
    _mcnamee_stenger_5,
    _mcnamee_stenger_7,
    _mcnamee_stenger_9,
)
from ._helpers import Enr2Scheme, integrate_monomial_over_enr2


def integrator(n, k):
    """Returns the integral of the polynomial given by the coefficients k over Enr2.
    """
    alpha = k + [0] * (n - len(k))
    return integrate_monomial_over_enr2(alpha, symbolic=True)


def mcnamee_stenger_3(n):
    return Enr2Scheme(*_mcnamee_stenger_3(n, integrator))


def mcnamee_stenger_5(n):
    return Enr2Scheme(*_mcnamee_stenger_5(n, integrator))


def mcnamee_stenger_7a(n):
    return Enr2Scheme(*_mcnamee_stenger_7(n, integrator, False))


def mcnamee_stenger_7b(n):
    return Enr2Scheme(*_mcnamee_stenger_7(n, integrator, True))


def mcnamee_stenger_9a(n):
    return Enr2Scheme(*_mcnamee_stenger_9(n, integrator, False))


def mcnamee_stenger_9b(n):
    return Enr2Scheme(*_mcnamee_stenger_9(n, integrator, True))

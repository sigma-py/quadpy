from ..cn._mcnamee_stenger import (
    _mcnamee_stenger_3,
    _mcnamee_stenger_5,
    _mcnamee_stenger_7,
    _mcnamee_stenger_9,
)
from ._helpers import SnScheme, integrate_monomial_over_unit_nball


def integrator(n, k, symbolic):
    """Returns the integral of the polynomial given by the coefficients k over NBall.
    """
    alpha = k + [0] * (n - len(k))
    return integrate_monomial_over_unit_nball(alpha, symbolic)


def mcnamee_stenger_3(n, symbolic=False):
    return SnScheme(*_mcnamee_stenger_3(n, integrator, symbolic=symbolic))


def mcnamee_stenger_5(n, symbolic=False):
    return SnScheme(*_mcnamee_stenger_5(n, integrator, symbolic=symbolic))


def mcnamee_stenger_7a(n, symbolic=False):
    return SnScheme(*_mcnamee_stenger_7(n, integrator, False, symbolic=symbolic))


def mcnamee_stenger_7b(n, symbolic=False):
    return SnScheme(*_mcnamee_stenger_7(n, integrator, True, symbolic=symbolic))


def mcnamee_stenger_9a(n, symbolic=False):
    return SnScheme(*_mcnamee_stenger_9(n, integrator, False, symbolic=symbolic))


def mcnamee_stenger_9b(n, symbolic=False):
    return SnScheme(*_mcnamee_stenger_9(n, integrator, True, symbolic=symbolic))

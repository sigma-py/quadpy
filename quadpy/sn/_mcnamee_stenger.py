import ndim

from ..cn._mcnamee_stenger import (
    _mcnamee_stenger_3,
    _mcnamee_stenger_5,
    _mcnamee_stenger_7,
    _mcnamee_stenger_9,
)
from ._helpers import SnScheme


def integrator(n, k, symbolic):
    """Returns the integral of the polynomial given by the coefficients k over NBall."""
    alpha = k + [0] * (n - len(k))
    return ndim.nball.integrate_monomial(alpha, symbolic=symbolic)


def mcnamee_stenger_3(n, symbolic=False):
    return SnScheme(*_mcnamee_stenger_3(n, integrator, symbolic=symbolic))


def mcnamee_stenger_5(n, symbolic=False):
    return SnScheme(*_mcnamee_stenger_5(n, integrator, symbolic=symbolic))


def mcnamee_stenger_7a(n, symbolic=False):
    return SnScheme(
        *_mcnamee_stenger_7(n, integrator, False, symbolic=symbolic), 1.453e-13
    )


def mcnamee_stenger_7b(n, symbolic=False):
    return SnScheme(
        *_mcnamee_stenger_7(n, integrator, True, symbolic=symbolic), 1.455e-14
    )


def mcnamee_stenger_9a(n, symbolic=False):
    return SnScheme(*_mcnamee_stenger_9(n, integrator, False, symbolic=symbolic))


def mcnamee_stenger_9b(n, symbolic=False):
    return SnScheme(*_mcnamee_stenger_9(n, integrator, True, symbolic=symbolic))

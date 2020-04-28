import warnings

from ..ncube._mcnamee_stenger import (
    _mcnamee_stenger_3,
    _mcnamee_stenger_5,
    _mcnamee_stenger_7,
    _mcnamee_stenger_9,
)
from ._helpers import EnrScheme, integrate_monomial_over_enr


def integrator(n, k):
    """Returns the integral of the polynomial given by the coefficients k over Enr.
    """
    alpha = k + [0] * (n - len(k))
    return integrate_monomial_over_enr(alpha, symbolic=True)


def mcnamee_stenger_3(n):
    return EnrScheme(*_mcnamee_stenger_3(n, integrator))


def mcnamee_stenger_5(n):
    return EnrScheme(*_mcnamee_stenger_5(n, integrator))


def mcnamee_stenger_7a(n):
    return EnrScheme(*_mcnamee_stenger_7(n, integrator, False))


def mcnamee_stenger_7b(n):
    return EnrScheme(*_mcnamee_stenger_7(n, integrator, True))


def mcnamee_stenger_9a(n):
    scheme = EnrScheme(*_mcnamee_stenger_9(n, integrator, False))
    warnings.warn("{} is very badly conditioned.".format(scheme.name))
    return scheme


def mcnamee_stenger_9b(n):
    scheme = EnrScheme(*_mcnamee_stenger_9(n, integrator, True))
    warnings.warn("{} is very badly conditioned.".format(scheme.name))
    return scheme

import warnings

from ..cn._mcnamee_stenger import (
    _mcnamee_stenger_3,
    _mcnamee_stenger_5,
    _mcnamee_stenger_7,
    _mcnamee_stenger_9,
)
from ._helpers import EnrScheme, integrate_monomial_over_enr


def integrator(n, k, symbolic):
    """Returns the integral of the polynomial given by the coefficients k over Enr.
    """
    alpha = k + [0] * (n - len(k))
    return integrate_monomial_over_enr(alpha, symbolic)


def mcnamee_stenger_3(n, symbolic=False):
    return EnrScheme(*_mcnamee_stenger_3(n, integrator, symbolic=symbolic))


def mcnamee_stenger_5(n, symbolic=False):
    return EnrScheme(*_mcnamee_stenger_5(n, integrator, symbolic=symbolic))


def mcnamee_stenger_7a(n, symbolic=False):
    return EnrScheme(*_mcnamee_stenger_7(n, integrator, False, symbolic=symbolic))


def mcnamee_stenger_7b(n, symbolic=False):
    return EnrScheme(*_mcnamee_stenger_7(n, integrator, True, symbolic=symbolic))


def mcnamee_stenger_9a(n, symbolic=False):
    scheme = EnrScheme(*_mcnamee_stenger_9(n, integrator, False, symbolic=symbolic))
    warnings.warn("{} is very badly conditioned.".format(scheme.name))
    return scheme


def mcnamee_stenger_9b(n, symbolic=False):
    scheme = EnrScheme(*_mcnamee_stenger_9(n, integrator, True, symbolic=symbolic))
    warnings.warn("{} is very badly conditioned.".format(scheme.name))
    return scheme

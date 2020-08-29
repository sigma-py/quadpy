import warnings

import ndim

from ..cn._mcnamee_stenger import (
    _mcnamee_stenger_3,
    _mcnamee_stenger_5,
    _mcnamee_stenger_7,
    _mcnamee_stenger_9,
)
from ._helpers import EnrScheme


def integrator(n, k, symbolic):
    """Returns the integral of the polynomial given by the coefficients k over Enr."""
    alpha = k + [0] * (n - len(k))
    return ndim.enr.integrate_monomial(alpha, symbolic)


def mcnamee_stenger_3(n, symbolic=False):
    return EnrScheme(*_mcnamee_stenger_3(n, integrator, symbolic=symbolic))


def mcnamee_stenger_5(n, symbolic=False):
    return EnrScheme(*_mcnamee_stenger_5(n, integrator, symbolic=symbolic))


def mcnamee_stenger_7a(n, symbolic=False):
    return EnrScheme(
        *_mcnamee_stenger_7(n, integrator, False, symbolic=symbolic), 7.530e-14
    )


def mcnamee_stenger_7b(n, symbolic=False):
    return EnrScheme(*_mcnamee_stenger_7(n, integrator, True, symbolic=symbolic))


def mcnamee_stenger_9a(n, symbolic=False):
    scheme = EnrScheme(*_mcnamee_stenger_9(n, integrator, False, symbolic=symbolic))
    warnings.warn(f"{scheme.name} is very badly conditioned.")
    return scheme


def mcnamee_stenger_9b(n, symbolic=False):
    scheme = EnrScheme(*_mcnamee_stenger_9(n, integrator, True, symbolic=symbolic))
    warnings.warn(f"{scheme.name} is very badly conditioned.")
    return scheme

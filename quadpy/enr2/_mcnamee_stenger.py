from sympy import sqrt, pi, factorial

from ._helpers import Enr2Scheme
from ..ncube._mcnamee_stenger import (
    _mcnamee_stenger_3,
    _mcnamee_stenger_5,
    _mcnamee_stenger_7,
    _mcnamee_stenger_9,
)


def integrator(n, k):
    """Returns the integral of the polynomial given by the coefficients k over Enr2.
    """
    assert len(k) <= n

    if any(kk % 2 == 1 for kk in k):
        return 0

    def prod(factors):
        from functools import reduce
        import operator
        return reduce(operator.mul, factors, 1)

    # return numpy.prod([math.gamma((kk + 1) / 2.0) for kk in k])
    k2 = [kk // 2 for kk in k]
    out = prod([factorial(2 * kk) / 4 ** kk / factorial(kk) * sqrt(pi) for kk in k2])
    out *= sqrt(pi) ** n
    return out


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

from typing import Callable

import numpy as np

from .c1 import integrate_adaptive


# compatibility for scipy.quad
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.quad.html
def quad(
    f: Callable,
    a,
    b,
    args=(),
    epsabs: float = 1.49e-08,
    epsrel: float = 1.49e-08,
    limit: int = 50,
):
    # See <https://www.gnu.org/software/gsl/doc/html/integration.html> for the
    # variable transformations
    if a == -np.inf and b == np.inf:
        # x = (1 - t) / t
        # dx / dt = -1 / t**2
        a = 0.0
        b = 1.0

        def g(t):
            return (f((1 - t) / t, *args) + f(-(1 - t) / t)) / t ** 2

    elif b == np.inf:
        a_orig = a
        a = 0.0
        b = 1.0

        def g(t):
            return f(a_orig + (1 - t) / t, *args) / t ** 2

    elif a == -np.inf:
        b_orig = b
        a = 0.0
        b = 1.0

        def g(t):
            return f(b_orig - (1 - t) / t, *args) / t ** 2

    else:

        def g(x):
            return f(x, *args)

    swap = a > b
    if swap:
        a, b = b, a

    val, err = integrate_adaptive(
        g,
        [a, b],
        eps_abs=epsabs,
        eps_rel=epsrel,
        criteria_connection=np.any,
        max_num_subintervals=limit,
    )
    if swap:
        val *= -1

    return val, err

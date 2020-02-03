import numpy

from .line_segment import integrate_adaptive


# compatibility for scipy.quad
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.quad.html
def quad(f, a, b, args=(), epsabs=1.49e-08, epsrel=1.49e-08):
    assert a <= b

    # See <https://www.gnu.org/software/gsl/doc/html/integration.html> for the
    # variable transformations
    if a == -numpy.inf and b == numpy.inf:
        # x = (1 - t) / t
        # dx / dt = -1 / t**2
        a = 0.0
        b = 1.0

        def g(t):
            return (f((1 - t) / t, *args) + f(-(1 - t) / t)) / t ** 2

    elif b == numpy.inf:
        a_orig = a
        a = 0.0
        b = 1.0

        def g(t):
            return f(a_orig + (1 - t) / t, *args) / t ** 2

    elif a == -numpy.inf:
        b_orig = b
        a = 0.0
        b = 1.0

        def g(t):
            return f(b_orig - (1 - t) / t, *args) / t ** 2

    else:

        def g(x):
            return f(x, *args)

    return integrate_adaptive(g, [a, b], eps_abs=epsabs, eps_rel=epsrel)

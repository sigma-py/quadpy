from .line_segment import integrate_adaptive


# compatibility for scipy.quad
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.quad.html
def quad(f, a, b, args=(), epsabs=1.49e-08):
    def g(x):
        return f(x, *args)
    return integrate_adaptive(g, [a, b], eps=epsabs)

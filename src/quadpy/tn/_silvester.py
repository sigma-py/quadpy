import numpy as np
import sympy

from ..helpers import article, get_all_exponents, prod
from ._helpers import TnScheme

source = article(
    authors=["P. Silvester"],
    title="Symmetric quadrature formulae for simplexes",
    journal="Math. Comp.",
    volume="24",
    pages="95-100",
    year="1970",
    url="https://doi.org/10.1090/S0025-5718-1970-0258283-6",
)


def integrate_bary(k, symbolic):
    # This is very similar to integrate_monomial_over_unit_simplex, except that n =
    # len(k) - 1, not len(k). This is because we're integrating over a polynomial in
    # barycentric coordinates here. Also the integration for [0, ..., 0] is set equal to
    # 1. This ensures that the weights add up to 1, not the volume of the unit simplex.
    frac = sympy.Rational if symbolic else lambda a, b: a / b

    n = len(k) - 1
    assert all(kk >= 0 for kk in k)

    if all(kk == 0 for kk in k):
        # return frac(1, math.factorial(n))
        return 1

    # find first nonzero
    idx = next(i for i, j in enumerate(k) if j > 0)
    alpha = frac(k[idx], sum(k) + n)
    k2 = k.copy()
    k2[idx] -= 1
    return integrate_bary(k2, symbolic) * alpha


def _get_data(dim, n, point_fun, symbolic):
    # points
    idxs = np.array(get_all_exponents(dim + 1, n)[-1])
    points = point_fun(idxs)

    # weights
    weights = np.empty(len(points))
    kk = 0
    for idx in idxs:
        # Define the polynomial which to integrate over the simplex.
        t = sympy.DeferredVector("t")
        g = prod(
            sympy.poly((t[i] - point_fun(k)) / (point_fun(m) - point_fun(k)))
            for i, m in enumerate(idx)
            for k in range(m)
        )
        if isinstance(g, int):
            exp = [0] * (dim + 1)
            coeffs = g
        else:
            # make sure the exponents are all of the correct length
            exp = [list(m) + [0] * (dim - len(m) + 1) for m in g.monoms()]
            coeffs = g.coeffs()
        weights[kk] = sum(c * integrate_bary(m, symbolic) for c, m in zip(coeffs, exp))
        kk += 1

    return weights, points


def silvester(dim, variant, n, symbolic=False):
    frac = np.vectorize(sympy.Rational) if symbolic else lambda a, b: a / b

    if variant == "closed":
        degree = n
        tol = 1.034e-14

        def points1d(k):
            return frac(k, n)

    else:
        assert variant == "open"
        degree = 1 if n == 0 else n
        tol = 1.040e-13

        def points1d(k):
            return frac(k + 1, n + 1 + dim)

    weights, points = _get_data(dim, n, points1d, symbolic)
    points = np.ascontiguousarray(points.T)

    return TnScheme(
        f"Silvester ({variant}, dim={n})", dim, weights, points, degree, source, tol
    )

import numpy
import sympy

from ..helpers import article, get_all_exponents, prod
from ._helpers import NSimplexScheme

citation = article(
    authors=["P. Silvester"],
    title="Symmetric quadrature formulae for simplexes",
    journal="Math. Comp.",
    volume="24",
    pages="95-100",
    year="1970",
    url="https://doi.org/10.1090/S0025-5718-1970-0258283-6",
)


# frac(prod([math.factorial(l) for l in m]), math.factorial(sum(m) + dim))
def integrate_bary_over_unit_simplex(k, symbolic):
    # This is very similar to integrate_monomial_over_unit_simplex, except that n =
    # len(k) - 1, not len(k). This is because we're integrating over a polynomial in
    # barycentric coordinates here. Also the integration for [0, ..., 0] is set equal to
    # 1. This ensure that the weights add up to 1, not the volume of the unit simplex.
    frac = sympy.Rational if symbolic else lambda a, b: a / b

    n = len(k) - 1
    assert all(kk >= 0 for kk in k)

    if all(kk == 0 for kk in k):
        # return frac(1, math.factorial(n))
        return 1

    # find first nonzero
    idx = next((i for i, j in enumerate(k) if j > 0), None)
    alpha = frac(k[idx], sum(k) + n)
    k2 = k.copy()
    k2[idx] -= 1
    return integrate_bary_over_unit_simplex(k2, symbolic) * alpha


def _newton_cotes(dim, n, point_fun, symbolic):
    degree = n

    # points
    idxs = numpy.array(get_all_exponents(dim + 1, n)[-1])
    points = point_fun(idxs)

    # weights
    if n == 0:
        weights = numpy.array([1])
        return points, weights, degree

    weights = numpy.empty(len(points))
    kk = 0
    for idx in idxs:
        # Define the polynomial which to integrate over the simplex.
        t = sympy.DeferredVector("t")
        g = prod(
            sympy.poly((t[i] - point_fun(k)) / (point_fun(m) - point_fun(k)))
            for i, m in enumerate(idx)
            for k in range(m)
        )
        # make sure the exponents are all of the correct length
        exp = [list(m) + [0] * (dim - len(m) + 1) for m in g.monoms()]
        weights[kk] = sum(
            c * integrate_bary_over_unit_simplex(m, symbolic)
            for m, c in zip(exp, g.coeffs())
        )
        kk += 1

    return weights, points, degree


def silvester(dim, variant, n, symbolic=False):
    frac = numpy.vectorize(sympy.Rational) if symbolic else lambda a, b: a / b

    if variant == "closed":
        weights, points, degree = _newton_cotes(dim, n, lambda k: frac(k, n), symbolic)
    else:
        assert variant == "open"
        weights, points, degree = _newton_cotes(
            dim, n, lambda k: frac(k + 1, n + 1 + dim), symbolic
        )
        if n == 0:
            degree = 1

    return NSimplexScheme(
        f"Silvester ({variant}, {n})", dim, weights, points, degree, citation
    )

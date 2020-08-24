import sympy

from ..helpers import article
from ._helpers import S3Scheme, register

_source = article(
    authors=["Preston C. Hammer", "Arthur H. Stroud"],
    title="Numerical Evaluation of Multiple Integrals II",
    journal="Math. Comp.",
    number="12",
    year="1958",
    pages="272-280",
    url="https://doi.org/10.1090/S0025-5718-1958-0102176-6",
)

frac = sympy.Rational
sqrt = sympy.sqrt
pi = sympy.pi


def hammer_stroud_11_3():
    d = {"symm_r00": [[frac(1, 6)], [sqrt(frac(3, 5))]]}
    return S3Scheme("Hammer-Stroud 11-3", d, 3, _source)


def hammer_stroud_12_3():
    alpha = sqrt(frac(3, 7))
    d = {
        "zero3": [[frac(1, 15)]],
        "symm_r00": [[frac(7, 90)], [alpha]],
        "symm_rr0": [[frac(7, 180)], [alpha]],
    }
    return S3Scheme("Hammer-Stroud 12-3", d, 5, _source)


def hammer_stroud_14_3(variant_a=True):
    t = 1 if variant_a else -1

    sqrt14 = sqrt(14)

    # ERR The article incorrectly gives 0.50824... instead of 0.050824...
    a1 = frac(1, 125) * (9 + t * 2 * sqrt14)
    c1 = (71 - t * 12 * sqrt14) / 1000

    nu = sqrt((7 - t * sqrt14) / 7)
    eta1 = sqrt(5 / (21 - t * 2 * sqrt14))

    d = {"symm_r00": [[a1], [nu]], "symm_rrr": [[c1], [eta1]]}
    name = "Hammer-Stroud 14-3" + ("a" if variant_a else "b")
    return S3Scheme(name, d, 5, _source)


def _hammer_stroud_15_3(variant_a):
    t = 1 if variant_a else -1

    sqrt30 = sqrt(30)
    nu2 = (45 - t * sqrt30) / 57
    xi2 = (18 + t * sqrt30) / 42
    eta2 = 7 / (27 + t * 2 * sqrt30)

    # The extract expressions are from Stroud's book.
    a1 = 1 / nu2 ** 3 / 63
    b1 = 1 / xi2 ** 3 / 630
    c1 = 1 / eta2 ** 3 / 2520
    a0 = 1 - 6 * a1 - 12 * b1 - 8 * c1

    d = {
        "zero3": [[a0]],
        "symm_r00": [[a1], [sqrt(nu2)]],
        "symm_rr0": [[b1], [sqrt(xi2)]],
        "symm_rrr": [[c1], [sqrt(eta2)]],
    }
    name = "Hammer-Stroud 15-3" + ("a" if variant_a else "b")
    return S3Scheme(name, d, 7, _source)


def hammer_stroud_15_3a():
    return _hammer_stroud_15_3(True)


def hammer_stroud_15_3b():
    return _hammer_stroud_15_3(False)


register(
    [
        hammer_stroud_11_3,
        hammer_stroud_12_3,
        hammer_stroud_14_3,
        hammer_stroud_15_3a,
        hammer_stroud_15_3b,
    ]
)

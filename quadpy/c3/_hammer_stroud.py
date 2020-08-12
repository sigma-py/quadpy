from sympy import Rational as frac
from sympy import sqrt

from ..helpers import article, fsd, pm, untangle, z
from ._helpers import C3Scheme, expand_symmetries

_source = article(
    authors=["Preston C. Hammer", "Arthur H. Stroud"],
    title="Numerical Evaluation of Multiple Integrals II",
    journal="Math. Comp.",
    number="12",
    year="1958",
    pages="272-280",
    url="https://doi.org/10.1090/S0025-5718-1958-0102176-6",
)


def hammer_stroud_1_3():
    d = {"symm_r00": [[frac(1, 6)], [1]]}
    points, weights = expand_symmetries(d)
    return C3Scheme("Hammer-Stroud 1-3", weights, points, 3, _source)


def hammer_stroud_2_3():
    alpha = sqrt(frac(3, 5))
    d = {
        "zero": [[frac(56, 27)]],
        "symm_r00": [[-frac(20, 81)], [alpha]],
        "symm_rr0": [[+frac(50, 81)], [alpha]],
    }
    points, weights = expand_symmetries(d)
    weights /= 8
    return C3Scheme("Hammer-Stroud 2-3", weights, points, 5, _source)


def hammer_stroud_4_3():
    d = {
        "symm_r00": [[frac(320, 361)], [sqrt(frac(19, 30))]],
        "symm_rrr": [[frac(121, 361)], [sqrt(frac(19, 33))]],
    }
    points, weights = expand_symmetries(d)
    weights /= 8
    return C3Scheme("Hammer-Stroud 4-3", weights, points, 5, _source)


def _hammer_stroud_5_3(variant_a):
    i = 1 if variant_a else -1

    r2 = (33 - i * sqrt(165)) / 28
    s2 = (30 + i * sqrt(165)) / 35
    t2 = (195 - i * 4 * sqrt(165)) / 337

    r = sqrt(r2)
    s = sqrt(s2)
    t = sqrt(t2)

    B1 = 176 / r2 ** 3 / 945
    B2 = 8 / s2 ** 3 / 135
    B3 = 8 / t2 ** 3 / 216
    B0 = 8 - 6 * B1 - 12 * B2 - 8 * B3

    d = {
        "zero": [[B0]],
        "symm_r00": [[B1], [r]],
        "symm_rr0": [[B2], [s]],
        "symm_rrr": [[B3], [t]],
    }
    points, weights = expand_symmetries(d)
    weights /= 8
    variant = "a" if variant_a else "b"
    return C3Scheme(f"Hammer-Stroud 5-3{variant}", weights, points, 7, _source)


def hammer_stroud_5_3a():
    return _hammer_stroud_5_3(True)


def hammer_stroud_5_3b():
    return _hammer_stroud_5_3(False)


def hammer_stroud_6_3():
    alpha = sqrt(frac(6, 7))
    d = {
        "symm_r00": [[frac(1078, 3645)], [alpha]],
        "symm_rr0": [[frac(343, 3645)], [alpha]],
        "symm_rrr": [
            [0.2247031747656014, 0.4123338622714356],
            [0.7341125287521153, 0.4067031864267161],
        ],
    }
    points, weights = expand_symmetries(d)
    weights /= 8
    return C3Scheme("Hammer-Stroud 6-3", weights, points, 7, _source)

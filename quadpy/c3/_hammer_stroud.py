from sympy import Rational as frac
from sympy import sqrt

from ..helpers import article, fsd, pm, untangle, z
from ._helpers import C3Scheme

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
    data = [(frac(1, 6), fsd(3, (1, 1)))]
    points, weights = untangle(data)
    return C3Scheme("Hammer-Stroud 1-3", weights, points, 3, _source)


def hammer_stroud_2_3():
    alpha = sqrt(frac(3, 5))
    data = [
        (+frac(56, 27), z(3)),
        (-frac(20, 81), fsd(3, (alpha, 1))),
        (+frac(50, 81), fsd(3, (alpha, 2))),
    ]
    points, weights = untangle(data)
    weights /= 8
    return C3Scheme("Hammer-Stroud 2-3", weights, points, 5, _source)


def hammer_stroud_4_3():
    data = [
        (frac(320, 361), fsd(3, (sqrt(frac(19, 30)), 1))),
        (frac(121, 361), pm(3 * [sqrt(frac(19, 33))])),
    ]
    points, weights = untangle(data)
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

    data = [
        (B0, z(3)),
        (B1, fsd(3, (r, 1))),
        (B2, fsd(3, (s, 2))),
        (B3, pm([t, t, t])),
    ]
    points, weights = untangle(data)
    weights /= 8
    variant = "a" if variant_a else "b"
    return C3Scheme(f"Hammer-Stroud 5-3{variant}", weights, points, 7, _source)


def hammer_stroud_5_3a():
    return _hammer_stroud_5_3(True)


def hammer_stroud_5_3b():
    return _hammer_stroud_5_3(False)


def hammer_stroud_6_3():
    alpha = sqrt(frac(6, 7))
    data = [
        (frac(1078, 3645), fsd(3, (alpha, 1))),
        (frac(343, 3645), fsd(3, (alpha, 2))),
        (0.2247031747656014, pm(3 * [0.7341125287521153])),
        (0.4123338622714356, pm(3 * [0.4067031864267161])),
    ]
    points, weights = untangle(data)
    weights /= 8
    return C3Scheme("Hammer-Stroud 6-3", weights, points, 7, _source)

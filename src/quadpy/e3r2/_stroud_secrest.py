from sympy import Rational as frac
from sympy import sqrt

from ..helpers import article
from ._helpers import E3r2Scheme, register

source = article(
    authors=["A.H. Stroud", "D. Secrest"],
    title="Approximate integration formulas for certain spherically symmetric regions",
    journal="Math. Comp.",
    volume="17",
    year="1963",
    pages="105-135",
    url="https://doi.org/10.1090/S0025-5718-1963-0161473-0",
)


def stroud_secrest_07():
    # ERR
    # article:
    # nu, xi = sqrt((15 + plus_minus * 3*sqrt(5)))
    # A = 3/5
    # B = 1/30

    # book:
    nu, xi = (sqrt((5 - p_m * sqrt(5)) / 4) for p_m in [+1, -1])
    A = frac(2, 5)
    B = frac(1, 20)

    d = {"zero3": [[A]], "symm_rs0_roll": [[B], [nu], [xi]]}
    return E3r2Scheme("Stroud-Secrest VII", d, 5, source)


def stroud_secrest_08a():
    r = sqrt(frac(5, 4))
    s = sqrt(frac(5, 2))
    d = {
        "symm_r00": [[frac(4, 25)], [r]],
        "symm_rrr": [[frac(1, 200)], [s]],
    }
    return E3r2Scheme("Stroud-Secrest VIIIa", d, 5, source)


def stroud_secrest_08b():
    r = sqrt(frac(5, 2))
    s = sqrt(frac(5, 6))
    d = {
        "zero3": [[frac(2, 5)]],
        "symm_r00": [[frac(1, 25)], [r]],
        "symm_rrr": [[frac(9, 200)], [s]],
    }
    return E3r2Scheme("Stroud-Secrest VIIIb", d, 5, source)


def stroud_secrest_09():
    r, s = (sqrt((15 - p_m * 5 * sqrt(5)) / 12) for p_m in [+1, -1])
    t = sqrt(frac(5, 6))
    d = {
        "zero3": [[frac(2, 5)]],
        "symm_rs0_roll": [[frac(3, 100)], [r], [s]],
        "symm_rrr": [[frac(3, 100)], [t]],
    }
    return E3r2Scheme("Stroud-Secrest IX", d, 5, source)


def _stroud_secrest_10(positive):
    plus_minus = 1 if positive else -1

    sqrt15 = sqrt(15)

    r = sqrt((15 + plus_minus * sqrt15) / 4)
    s = sqrt((6 - plus_minus * sqrt15) / 2)
    t = sqrt((9 + plus_minus * 2 * sqrt15) / 2)
    A = (720 + plus_minus * 8 * sqrt15) / 2205
    B = (270 - plus_minus * 46 * sqrt15) / 15435
    C = (162 + plus_minus * 41 * sqrt15) / 6174
    D = (783 - plus_minus * 202 * sqrt15) / 24696

    d = {
        "zero3": [[A]],
        "symm_r00": [[B], [r]],
        "symm_rr0": [[C], [s]],
        "symm_rrr": [[D], [t]],
    }
    variant = "a" if positive else "b"
    return E3r2Scheme(f"Stroud-Secrest X{variant}", d, 7, source)


def stroud_secrest_10a():
    return _stroud_secrest_10(True)


def stroud_secrest_10b():
    return _stroud_secrest_10(False)


def _stroud_secrest_11(positive):
    p_m = 1 if positive else -1

    sqrt2 = sqrt(2)
    sqrt5 = sqrt(5)
    sqrt10 = sqrt(10)

    r = sqrt((25 + p_m * 15 * sqrt2 + 5 * sqrt5 + p_m * 3 * sqrt10) / 4)
    s = sqrt((25 + p_m * 15 * sqrt2 - 5 * sqrt5 - p_m * 3 * sqrt10) / 4)
    t = sqrt((3 - p_m * sqrt2) / 2)
    u = sqrt((9 - p_m * 3 * sqrt2 - 3 * sqrt5 + p_m * sqrt10) / 4)
    v = sqrt((9 - p_m * 3 * sqrt2 + 3 * sqrt5 - p_m * sqrt10) / 4)

    A = (80 + p_m * 8 * sqrt2) / 245
    B = (395 - p_m * 279 * sqrt2) / 13720
    C = (45 + p_m * 29 * sqrt2) / 2744

    d = {
        "zero3": [[A]],
        "symm_rs0_roll": [[B, C], [r, u], [s, v]],
        "symm_rrr": [[C], [t]],
    }
    return E3r2Scheme("Stroud-Secrest XI", d, 7, source)


def stroud_secrest_11a():
    return _stroud_secrest_11(True)


def stroud_secrest_11b():
    return _stroud_secrest_11(False)


register(
    [
        stroud_secrest_07,
        stroud_secrest_09,
        stroud_secrest_08a,
        stroud_secrest_08b,
        stroud_secrest_10a,
        stroud_secrest_10b,
        stroud_secrest_11a,
        stroud_secrest_11b,
    ]
)

from sympy import Rational as frac

from ..helpers import article
from ._helpers import C3Scheme, register

_source = article(
    authors=["D. Mustard", "J.N. Lyness", "J.M. Blatt"],
    title="Numerical quadrature in n dimensions",
    journal="The Computer Journal",
    year="1963",
    volume="6",
    number="1",
    pages="75-87",
    url="https://doi.org/10.1093/comjnl/6.1.75",
)


def mustard_lyness_blatt_1():
    d = {"zero3": [[frac(1, 2)]], "symm_rr0": [[frac(1, 24)], [1]]}
    return C3Scheme("Mustard-Lyness-Blatt 1", d, 3, _source)


def mustard_lyness_blatt_2():
    d = {
        "zero3": [[frac(2, 9)]],
        "symm_r00": [[frac(1, 9)], [1]],
        "symm_rrr": [[frac(1, 72)], [1]],
    }
    return C3Scheme("Mustard-Lyness-Blatt 2", d, 3, _source)


def mustard_lyness_blatt_3():
    d = {
        "symm_rr0": [[frac(1, 6)], [1]],
        "symm_rrr": [[-frac(1, 8)], [1]],
    }
    return C3Scheme("Mustard-Lyness-Blatt 3", d, 3, _source)


def mustard_lyness_blatt_4():
    d = {
        "zero3": [[-frac(2, 45)]],
        "symm_r00": [[+frac(2, 45)], [1]],
        "symm_rrr": [[frac(4, 45), frac(1, 120)], [frac(1, 2), 1]],
    }
    return C3Scheme("Mustard-Lyness-Blatt 4", d, 5, _source)


def mustard_lyness_blatt_5():
    d = {
        "zero3": [[-frac(19, 15)]],
        "symm_r00": [[frac(16, 45), -frac(1, 30)], [frac(1, 2), 1]],
        "symm_rr0": [[frac(1, 36)], [1]],
    }
    return C3Scheme("Mustard-Lyness-Blatt 5", d, 5, _source)


def mustard_lyness_blatt_6():
    d = {
        "zero3": [[-frac(4, 3)]],
        "symm_r00": [[frac(16, 45)], [frac(1, 2)]],
        "symm_rr0": [[frac(1, 90)], [1]],
        "symm_rrr": [[frac(1, 120)], [1]],
    }
    return C3Scheme("Mustard-Lyness-Blatt 6", d, 5, _source)


def mustard_lyness_blatt_7():
    d = {
        "zero3": [[frac(2, 45)]],
        "symm_rr0": [[frac(1, 45)], [1]],
        "symm_rrr": [[frac(4, 45), frac(-1, 360)], [frac(1, 2), 1]],
    }
    return C3Scheme("Mustard-Lyness-Blatt 7", d, 5, _source)


register(
    [
        mustard_lyness_blatt_1,
        mustard_lyness_blatt_2,
        mustard_lyness_blatt_3,
        mustard_lyness_blatt_4,
        mustard_lyness_blatt_5,
        mustard_lyness_blatt_6,
        mustard_lyness_blatt_7,
    ]
)

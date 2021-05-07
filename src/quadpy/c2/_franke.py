import math
import warnings

from sympy import Rational as frac
from sympy import sqrt

from ..helpers import article
from ._helpers import C2Scheme, register
from ._rabinowitz_richter import rabinowitz_richter_1
from ._tyler import tyler_2

source = article(
    authors=["Richard Franke"],
    title="Obtaining cubatures for rectangles and other planar regions by using orthogonal polynomials",
    journal="Math. Comp.",
    volume="25",
    year="1971",
    pages="803-817",
    url="https://doi.org/10.1090/S0025-5718-1971-0300440-5",
)


def franke_1(lmbda=2):
    assert -frac(9, 5) <= lmbda <= frac(9, 4)

    a = sqrt(frac(9 + 5 * lmbda, 15))
    b = sqrt(frac(9 - 4 * lmbda, 15))
    c = sqrt(frac(3, 5))

    d = {
        "zero2": [[frac(4 * (4 + 5 * lmbda), 9 * (9 + 5 * lmbda))]],
        "sxy": [[frac(25, 36 * (9 - 4 * lmbda))], [b], [c]],
        "c2_a0": [[frac(40, 36 * (9 + 5 * lmbda))], [a]],
        "c2_0a": [[frac(10 * (1 - lmbda), 9 * (9 - 4 * lmbda))], [c]],
    }
    return C2Scheme(f"Franke(1, {lmbda})", d, 5, source)


def franke_2a():
    a = math.sqrt((15 + 2 * sqrt(30)) / 35)
    b = math.sqrt((15 - 2 * sqrt(30)) / 35)

    # closed forms not appearing in the article:
    c = 5 * (18 + math.sqrt(30)) / 324 / 4
    d = 2 * (18 + math.sqrt(30)) / 81 / 4

    d = {
        "sxy": [
            [0.01094603802180727, c, 0.07601767326255625],
            [0.105784012371275e1, math.sqrt(3 / 5), 0.469253522127911],
            [a, b, a],
        ],
        "c2_0a": [[d], [b]],
    }
    return C2Scheme("Franke 2a", d, 7, source, 1.232e-14)


def franke_2b():
    a = math.sqrt((15 + 2 * sqrt(30)) / 35)
    b = math.sqrt((15 - 2 * sqrt(30)) / 35)

    d = {
        "sxy": [
            [0.0483131729357575, 0.0422624803047505, 0.120773808410886],
            [0.774596669241483, 0.915060523380880, 0.396191039748320],
            [a, b, b],
        ],
        "c2_0a": [[0.077301076697212], [a]],
    }
    return C2Scheme("Franke 2b", d, 7, source)


def franke_3a():
    a = math.sqrt(5.0 / 9.0 + 2.0 / 63.0 * math.sqrt(70))
    b = math.sqrt(5.0 / 9.0 - 2.0 / 63.0 * math.sqrt(70))

    d = {
        "sxy": [
            [
                0.705065140564012e-1,
                0.721121511007611e-1,
                0.971492736037507e-1,
                0.368549048677049,
            ],
            [
                0.845927799771709,
                0.628901636732253,
                0.959681421214621,
                0.436030596273468,
            ],
            [a, a, b, b],
        ],
        "c2_a0": [[0.316049382716049], [0.774596669241483]],
        "c2_0a": [[0.188616439798053, 0.258606964371341e-1], [a, b]],
        "zero2": [[0.505679012345679]],
    }

    for value in d.values():
        value[0] = [val / 4 for val in value[0]]

    return C2Scheme("Franke 3a", d, 9, source)


def franke_3b():
    a = math.sqrt(5.0 / 9.0 + 2.0 / 63.0 * math.sqrt(70))
    b = math.sqrt(5.0 / 9.0 - 2.0 / 63.0 * math.sqrt(70))

    d = {
        "sxy": [
            [
                0.499290623065150e-1,
                0.158445182284802,
                0.183383788151247,
                0.881476523665422e-1,
            ],
            [
                0.945813739519925,
                0.465346624836203,
                0.804253925742002,
                0.681385892163677,
            ],
            [a, a, b, b],
        ],
        "c2_a0": [
            [0.114456375561331, 0.454432513327558],
            [0.963018409085396, 0.428610143223121],
        ],
        "c2_0a": [[0.571052809297435e-1, 0.414194459963155], [a, b]],
    }
    for value in d.values():
        value[0] = [val / 4 for val in value[0]]
    return C2Scheme("Franke 3b", d, 9, source)


def franke_3c():
    a = math.sqrt(5.0 / 9.0 + 2.0 / 63.0 * math.sqrt(70))
    b = math.sqrt(5.0 / 9.0 - 2.0 / 63.0 * math.sqrt(70))

    d = {
        "sxy": [
            [0.494522019130682e-1, 0.163914731881061, 0.265904816944092],
            [0.949307350001342, 0.458177548931134, 0.774596669241483],
            [a, a, b],
        ],
        "c2_a0": [
            [0.113041839046410, 0.479922229600720],
            [0.967776908976724, 0.417754671502987],
        ],
        "c2_0a": [[0.471199025241204e-1, 0.425447707110548], [a, b]],
        "zero2": [[-0.481503595164821e-1]],
    }
    for value in d.values():
        value[0] = [val / 4 for val in value[0]]
    return C2Scheme("Franke 3c", d, 9, source)


def franke_5():
    # DUP Noted as a duplicate in the original article
    return tyler_2()


def franke_6():
    a = sqrt(frac(3, 2))
    b = sqrt(frac(3, 7) * (1 + sqrt(frac(10, 31))))
    c = sqrt(frac(3, 7) * (1 - sqrt(frac(10, 31))))

    d = {
        "zero2": [[frac(392, 405)]],
        "d4_aa": [[frac(16, 2025)], [a]],
        "d4_ab": [[frac(1519, 4050)], [b], [c]],
    }
    for value in d.values():
        value[0] = [val / 4 for val in value[0]]
    return C2Scheme("Franke 6", d, 7, source)


def franke_7():
    # DUP Noted as a duplicate in the original article
    return rabinowitz_richter_1()


def franke_8():
    # TODO find error in franke_8
    warnings.warn("Franke(8) only has degree 1.")

    d = {
        "d4_a0": [[0.454163960686749], [0.488926856974369]],
        "d4_aa": [
            [0.214200360926862, 0.427312318657758e-1],
            [0.690880550486344, 0.939565258096838],
        ],
        "d4_ab": [[0.144452223260307], [0.918620441056722], [0.344872025364404]],
    }
    for value in d.values():
        value[0] = [val / 4 for val in value[0]]
    return C2Scheme("Franke 8", d, 1, source)


register(
    [
        franke_1,
        franke_2a,
        franke_2b,
        franke_3a,
        franke_3b,
        franke_3c,
        franke_5,
        franke_6,
        franke_7,
        franke_8,
    ]
)

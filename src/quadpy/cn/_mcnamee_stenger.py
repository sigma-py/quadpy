import math

import numpy as np
import sympy

from ..helpers import article, comb, fsd, prod, untangle, z
from ._helpers import CnScheme

_source = article(
    authors=["J. McNamee", "F. Stenger"],
    title="Construction of Fully Symmetric Numerical Integration Formulas",
    journal="Numerische Mathematik",
    year="1967",
    volume="10",
    pages="327-344",
    url="https://doi.org/10.1007/BF02162032",
)


def _mcnamee_stenger_3(n, integrator, symbolic):
    sqrt = sympy.sqrt if symbolic else math.sqrt

    I0 = integrator(n, [], symbolic)
    I2 = integrator(n, [2], symbolic)

    u2 = I2 / I0
    # ERR The article says I0 / (2 * u**2)
    A1 = I2 / (2 * u2)
    A0 = (1 - n) * I0
    u = sqrt(u2)

    data = [
        (A0 / I0, z(n)),
        (A1 / I0, fsd(n, (u, 1))),
    ]

    points, weights = untangle(data)
    points = np.ascontiguousarray(points.T)
    return "McNamee-Stenger 3", n, weights, points, 3, _source


def _mcnamee_stenger_5(n, integrator, symbolic):
    assert n >= 2
    frac = sympy.Rational if symbolic else lambda a, b: a / b
    sqrt = sympy.sqrt if symbolic else math.sqrt

    I0 = integrator(n, [], symbolic)
    I2 = integrator(n, [2], symbolic)
    I4 = integrator(n, [4], symbolic)
    I22 = integrator(n, [2, 2], symbolic)

    u = sqrt(I4 / I2)
    A0 = I0 - n * (I2 / I4) ** 2 * (I4 - frac(n - 1, 2) * I22)
    A1 = frac(1, 2) * (I2 / I4) ** 2 * (I4 - (n - 1) * I22)
    A11 = frac(1, 4) * (I2 / I4) ** 2 * I22

    data = [
        (A0, z(n)),
        (A1, fsd(n, (u, 1))),
        (A11, fsd(n, (u, 2))),
    ]

    points, weights = untangle(data)
    points = np.ascontiguousarray(points.T)
    weights /= I0
    return "McNamee-Stenger 5", n, weights, points, 5, _source


def _mcnamee_stenger_7(n, integrator, switch_uv, symbolic):
    assert n >= 3
    frac = sympy.Rational if symbolic else lambda a, b: a / b
    sqrt = sympy.sqrt if symbolic else math.sqrt

    I0 = integrator(n, [], symbolic)
    I2 = integrator(n, [2], symbolic)
    I4 = integrator(n, [4], symbolic)
    I6 = integrator(n, [6], symbolic)
    I22 = integrator(n, [2, 2], symbolic)
    I24 = integrator(n, [2, 4], symbolic)
    I222 = integrator(n, [2, 2, 2], symbolic)

    # Choose u, v as solutions of a u**4 - b u**2 + c = 0.
    a = I2 ** 2 - I0 * I4
    b = I2 * I4 - I0 * I6
    c = I4 ** 2 - I2 * I6
    #
    u2 = (b + sqrt(b ** 2 - 4 * a * c)) / (2 * a)
    v2 = (b - sqrt(b ** 2 - 4 * a * c)) / (2 * a)

    if switch_uv:
        u2, v2 = v2, u2

    u4 = u2 ** 2
    u6 = u2 ** 3
    v4 = v2 ** 2
    v6 = v2 ** 3

    A111 = I222 / (2 * u2) ** 3

    # mat = [[u4, v4], [u6, v6]]
    det = u4 * v6 - v4 * u6
    inv = [[v6 / det, -v4 / det], [-u6 / det, u4 / det]]
    # vector
    vec = [
        I22 - 2 ** 3 * (n - 2) * u4 * A111,
        I24 - 2 ** 3 * (n - 2) * u6 * A111,
    ]
    A11 = frac(1, 2 ** 2) * (inv[0][0] * vec[0] + inv[0][1] * vec[1])
    A22 = frac(1, 2 ** 2) * (inv[1][0] * vec[0] + inv[1][1] * vec[1])

    # mat = [[u2, v2], [u4, v4]]
    det = u2 * v4 - v2 * u4
    inv = [[v4 / det, -v2 / det], [-u4 / det, u2 / det]]
    vec = [
        I2 - 2 ** 3 * comb(n - 1, 2) * u2 * A111,
        I4 - 2 ** 3 * comb(n - 1, 2) * u4 * A111,
    ]
    A1 = -2 * (n - 1) * A11 + frac(1, 2) * (inv[0][0] * vec[0] + inv[0][1] * vec[1])
    A2 = -2 * (n - 1) * A22 + frac(1, 2) * (inv[1][0] * vec[0] + inv[1][1] * vec[1])

    A0 = (
        I0
        - 2 * n * (A1 + A2)
        - 2 ** 2 * comb(n, 2) * (A11 + A22)
        - 2 ** 3 * comb(n, 3) * A111
    )

    u = sqrt(u2)
    v = sqrt(v2)

    data = [
        (A0, z(n)),
        (A1, fsd(n, (u, 1))),
        (A2, fsd(n, (v, 1))),
        (A11, fsd(n, (u, 2))),
        (A22, fsd(n, (v, 2))),
        (A111, fsd(n, (u, 3))),
    ]

    points, weights = untangle(data)
    points = np.ascontiguousarray(points.T)

    weights /= I0
    variant = "a" if not switch_uv else "b"
    return f"McNamee-Stenger 7{variant}", n, weights, points, 7, _source


def _mcnamee_stenger_9(n, integrator, switch_uv, symbolic):
    assert n >= 4
    sqrt = sympy.sqrt if symbolic else math.sqrt
    frac = sympy.Rational if symbolic else lambda a, b: a / b

    I0 = integrator(n, [], symbolic)
    I2 = integrator(n, [2], symbolic)
    I4 = integrator(n, [4], symbolic)
    I6 = integrator(n, [6], symbolic)
    I8 = integrator(n, [8], symbolic)
    I24 = integrator(n, [2, 4], symbolic)
    I26 = integrator(n, [2, 6], symbolic)
    I44 = integrator(n, [4, 4], symbolic)
    I222 = integrator(n, [2, 2, 2], symbolic)
    I224 = integrator(n, [2, 2, 4], symbolic)
    I2222 = integrator(n, [2, 2, 2, 2], symbolic)

    # Choose u, v as solutions of a u**4 - b u**2 + c = 0.
    a = I4 ** 2 - I2 * I6
    b = I4 * I6 - I2 * I8
    c = I6 ** 2 - I4 * I8
    #
    u2 = (b + sqrt(b ** 2 - 4 * a * c)) / (2 * a)
    v2 = (b - sqrt(b ** 2 - 4 * a * c)) / (2 * a)

    if switch_uv:
        u2, v2 = v2, u2

    u4 = u2 ** 2
    u6 = u2 ** 3
    u8 = u2 ** 4
    v4 = v2 ** 2
    v6 = v2 ** 3
    v8 = v2 ** 4

    A1111 = I2222 / (2 * u2) ** 4

    # mat = [[u6, v6], [u8, v8]]
    det = u6 * v8 - v6 * u8
    inv = [[v8 / det, -v6 / det], [-u8 / det, u6 / det]]
    # vector
    vec = [
        I222 - 2 ** 4 * (n - 3) * u6 * A1111,
        I224 - 2 ** 4 * (n - 3) * u8 * A1111,
    ]
    A111 = frac(1, 2 ** 3) * (inv[0][0] * vec[0] + inv[0][1] * vec[1])
    A222 = frac(1, 2 ** 3) * (inv[1][0] * vec[0] + inv[1][1] * vec[1])

    A12 = (I26 - I44) / (2 ** 2 * u2 * v2 * (u2 - v2) ** 2)

    vec = [
        I24 - 2 ** 2 * (u4 * v2 + u2 * v4) * A12 - 2 ** 4 * comb(n - 2, 2) * u6 * A1111,
        I26 - 2 ** 2 * (u6 * v2 + u2 * v6) * A12 - 2 ** 4 * comb(n - 2, 2) * u8 * A1111,
    ]
    A11 = -2 * (n - 2) * A111 + frac(1, 4) * (inv[0][0] * vec[0] + inv[0][1] * vec[1])
    A22 = -2 * (n - 2) * A222 + frac(1, 4) * (inv[1][0] * vec[0] + inv[1][1] * vec[1])

    # mat = [[u2, v2], [u4, v4]]
    det = u2 * v4 - v2 * u4
    inv = [[v4 / det, -v2 / det], [-u4 / det, u2 / det]]
    #
    vec = [
        I2 - 2 ** 4 * comb(n - 1, 3) * u2 * A1111,
        I4 - 2 ** 4 * comb(n - 1, 3) * u4 * A1111,
    ]
    A1 = (
        -2 * (n - 1) * (A11 + A12)
        - 2 ** 2 * comb(n - 1, 2) * A111
        + frac(1, 2) * (inv[0][0] * vec[0] + inv[0][1] * vec[1])
    )
    A2 = (
        -2 * (n - 1) * (A22 + A12)
        - 2 ** 2 * comb(n - 1, 2) * A222
        + frac(1, 2) * (inv[1][0] * vec[0] + inv[1][1] * vec[1])
    )

    A0 = (
        I0
        - 2 * n * (A1 + A2)
        - 2 ** 2 * comb(n, 2) * (A11 + 2 * A12 + A22)
        - 2 ** 3 * comb(n, 3) * (A111 + A222)
        - 2 ** 4 * comb(n, 4) * A1111
    )

    u = sqrt(u2)
    v = sqrt(v2)

    data = [
        (A0, z(n)),
        (A1, fsd(n, (u, 1))),
        (A2, fsd(n, (v, 1))),
        (A11, fsd(n, (u, 2))),
        (A12, fsd(n, (u, 1), (v, 1))),
        (A22, fsd(n, (v, 2))),
        (A111, fsd(n, (u, 3))),
        (A222, fsd(n, (v, 3))),
        (A1111, fsd(n, (u, 4))),
    ]

    points, weights = untangle(data)
    points = np.ascontiguousarray(points.T)

    weights /= I0
    variant = "a" if not switch_uv else "b"
    return f"McNamee-Stenger 9{variant}", n, weights, points, 9, _source


# Here starts the cube-specific section.
def integrator(n, k, symbolic):
    """Returns the integral of the polynomial given by the coefficients k over the
    n-dimensional cube.
    """
    assert len(k) <= n

    if any(kk % 2 == 1 for kk in k):
        return 0

    frac = sympy.Rational if symbolic else lambda a, b: a / b
    return frac(2 ** n, prod([kk + 1 for kk in k]))


def mcnamee_stenger_3(n, symbolic=False):
    return CnScheme(*_mcnamee_stenger_3(n, integrator, symbolic=symbolic), 6.374e-14)


def mcnamee_stenger_5(n, symbolic=False):
    return CnScheme(*_mcnamee_stenger_5(n, integrator, symbolic=symbolic), 4.187e-14)


def mcnamee_stenger_7a(n, symbolic=False):
    return CnScheme(
        *_mcnamee_stenger_7(n, integrator, False, symbolic=symbolic), 4.748e-11
    )


def mcnamee_stenger_7b(n, symbolic=False):
    return CnScheme(
        *_mcnamee_stenger_7(n, integrator, True, symbolic=symbolic), 1.669e-12
    )


def mcnamee_stenger_9a(n, symbolic=False):
    return CnScheme(
        *_mcnamee_stenger_9(n, integrator, False, symbolic=symbolic), 3.961e-12
    )


def mcnamee_stenger_9b(n, symbolic=False):
    return CnScheme(
        *_mcnamee_stenger_9(n, integrator, True, symbolic=symbolic), 1.512e-12
    )

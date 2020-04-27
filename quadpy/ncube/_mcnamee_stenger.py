from sympy import Rational as frac
from sympy import binomial, sqrt

from ..helpers import article, fsd, untangle, z
from ._helpers import NCubeScheme

_citation = article(
    authors=["J. McNamee", "F. Stenger"],
    title="Construction of Fully Symmetric Numerical Integration Formulas",
    journal="Numerische Mathematik",
    year="1967",
    volume="10",
    pages="327-344",
    url="https://doi.org/10.1007/BF02162032",
)


def mcnamee_stenger_3(n):
    I0 = 2 ** n
    I2 = frac(2 ** n, 3)

    u2 = frac(I2, I0)
    # ERR The article says I0 / (2 * u**2)
    A1 = frac(I2, 2 * u2)
    A0 = (1 - n) * I0
    u = sqrt(u2)

    data = [
        (A0, z(n)),
        (A1, fsd(n, (u, 1))),
    ]

    points, weights = untangle(data)
    # weights *= 2 ** n
    return NCubeScheme("McNamee-Stenger 3", n, weights, points, 3, _citation)


def mcnamee_stenger_5(n):
    assert n >= 2

    I0 = 2 ** n
    I2 = frac(2 ** n, 3)
    I4 = frac(2 ** n, 5)
    I22 = frac(2 ** n, 9)

    u = sqrt(frac(I4, I2))
    A0 = I0 - n * (I2 / I4) ** 2 * (I4 - frac(n - 1, 2) * I22)
    A1 = frac(1, 2) * frac(I2, I4) ** 2 * (I4 - (n - 1) * I22)
    A11 = frac(1, 4) * frac(I2, I4) ** 2 * I22

    data = [
        (A0, z(n)),
        (A1, fsd(n, (u, 1))),
        (A11, fsd(n, (u, 2))),
    ]

    points, weights = untangle(data)
    # weights *= 2 ** n
    return NCubeScheme("McNamee-Stenger 5", n, weights, points, 5, _citation)


def _mcnamee_stenger_7(n, switch_uv):
    assert n >= 3

    I0 = 2 ** n
    I2 = frac(2 ** n, 3)
    I4 = frac(2 ** n, 5)
    I6 = frac(2 ** n, 7)
    I22 = frac(2 ** n, 9)
    I24 = frac(2 ** n, 15)
    I222 = frac(2 ** n, 27)

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
        I2 - 2 ** 3 * binomial(n - 1, 2) * u2 * A111,
        I4 - 2 ** 3 * binomial(n - 1, 2) * u4 * A111,
    ]
    A1 = -2 * (n - 1) * A11 + frac(1, 2) * (inv[0][0] * vec[0] + inv[0][1] * vec[1])
    A2 = -2 * (n - 1) * A22 + frac(1, 2) * (inv[1][0] * vec[0] + inv[1][1] * vec[1])

    A0 = (
        I0
        - 2 * n * (A1 + A2)
        - 2 ** 2 * binomial(n, 2) * (A11 + A22)
        - 2 ** 3 * binomial(n, 3) * A111
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

    # weights *= 2 ** n
    return NCubeScheme("McNamee-Stenger 7", n, weights, points, 7, _citation)


def mcnamee_stenger_7a(n):
    return _mcnamee_stenger_7(n, switch_uv=False)


def mcnamee_stenger_7b(n):
    return _mcnamee_stenger_7(n, switch_uv=True)


def _mcnamee_stenger_9(n, switch_uv):
    assert n >= 4

    I0 = 2 ** n
    I2 = frac(2 ** n, 3)
    I4 = frac(2 ** n, 5)
    I6 = frac(2 ** n, 7)
    I8 = frac(2 ** n, 9)
    I24 = frac(2 ** n, 15)
    I26 = frac(2 ** n, 21)
    I44 = frac(2 ** n, 25)
    I222 = frac(2 ** n, 27)
    I224 = frac(2 ** n, 45)
    I2222 = frac(2 ** n, 81)

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
        I24
        - 2 ** 2 * (u4 * v2 + u2 * v4) * A12
        - 2 ** 4 * binomial(n - 2, 2) * u6 * A1111,
        I26
        - 2 ** 2 * (u6 * v2 + u2 * v6) * A12
        - 2 ** 4 * binomial(n - 2, 2) * u8 * A1111,
    ]
    A11 = -2 * (n - 2) * A111 + frac(1, 4) * (inv[0][0] * vec[0] + inv[0][1] * vec[1])
    A22 = -2 * (n - 2) * A222 + frac(1, 4) * (inv[1][0] * vec[0] + inv[1][1] * vec[1])

    # mat = [[u2, v2], [u4, v4]]
    det = u2 * v4 - v2 * u4
    inv = [[v4 / det, -v2 / det], [-u4 / det, u2 / det]]
    #
    vec = [
        I2 - 2 ** 4 * binomial(n - 1, 3) * u2 * A1111,
        I4 - 2 ** 4 * binomial(n - 1, 3) * u4 * A1111,
    ]
    A1 = (
        -2 * (n - 1) * (A11 + A12)
        - 2 ** 2 * binomial(n - 1, 2) * A111
        + frac(1, 2) * (inv[0][0] * vec[0] + inv[0][1] * vec[1])
    )
    A2 = (
        -2 * (n - 1) * (A22 + A12)
        - 2 ** 2 * binomial(n - 1, 2) * A222
        + frac(1, 2) * (inv[1][0] * vec[0] + inv[1][1] * vec[1])
    )

    A0 = (
        I0
        - 2 * n * (A1 + A2)
        - 2 ** 2 * binomial(n, 2) * (A11 + 2 * A12 + A22)
        - 2 ** 3 * binomial(n, 3) * (A111 + A222)
        - 2 ** 4 * binomial(n, 4) * A1111
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

    # weights *= 2 ** n
    return NCubeScheme("McNamee-Stenger 9", n, weights, points, 9, _citation)


def mcnamee_stenger_9a(n):
    return _mcnamee_stenger_9(n, switch_uv=False)


def mcnamee_stenger_9b(n):
    return _mcnamee_stenger_9(n, switch_uv=True)
from sympy import Rational as frac
from sympy import sqrt

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

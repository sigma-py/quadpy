import numpy as np

from ..helpers import article
from ._helpers import S2Scheme, register

_source = article(
    authors=["Frank G. Lether"],
    title="A Generalized Product Rule for the Circle",
    journal="SIAM Journal on Numerical Analysis",
    volume="8",
    number="2",
    month="jun",
    year="1971",
    pages="249-253",
    url="https://www.jstor.org/stable/2949473",
)


def lether(n: int):
    assert n >= 1
    p, w = np.polynomial.legendre.leggauss(n)

    mu = np.arange(1, n + 1)
    points = np.array(
        [
            np.tile(np.cos(mu * np.pi / (n + 1)), n),
            np.outer(p, np.sin(mu * np.pi / (n + 1))).flatten(),
        ]
    )
    weights = np.outer(w, np.sin(mu * np.pi / (n + 1)) ** 2).flatten() / (n + 1)
    d = {"plain": [weights, points[0], points[1]]}
    return S2Scheme(f"Lether({n})", d, 2 * n - 1, _source)


register([lether])

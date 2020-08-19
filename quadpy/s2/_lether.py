import numpy

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


def lether(n):
    assert n >= 1
    p, w = numpy.polynomial.legendre.leggauss(n)

    mu = numpy.arange(1, n + 1)
    points = numpy.array(
        [
            numpy.tile(numpy.cos(mu * numpy.pi / (n + 1)), n),
            numpy.outer(p, numpy.sin(mu * numpy.pi / (n + 1))).flatten(),
        ]
    )
    weights = numpy.outer(w, numpy.sin(mu * numpy.pi / (n + 1)) ** 2).flatten() / (
        n + 1
    )
    d = {"plain": [weights, points[0], points[1]]}
    return S2Scheme(f"Lether({n})", d, 2 * n - 1, _source)


register([lether])

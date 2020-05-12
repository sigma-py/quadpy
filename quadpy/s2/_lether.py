import numpy

from ..helpers import article
from ._helpers import DiskScheme

_citation = article(
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
    points = numpy.column_stack(
        [
            numpy.tile(numpy.cos(mu * numpy.pi / (n + 1)), n),
            numpy.outer(p, numpy.sin(mu * numpy.pi / (n + 1))).flatten(),
        ]
    )

    weights = (
        numpy.pi
        / (n + 1)
        * numpy.outer(w, numpy.sin(mu * numpy.pi / (n + 1)) ** 2).flatten()
    )
    return DiskScheme(f"Lether({n})", weights, points, 2 * n - 1, _citation)

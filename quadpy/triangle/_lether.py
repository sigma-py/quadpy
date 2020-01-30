# The results have later been reproduced in:
#
# H.T. Rathod, K.V. Nagaraja, B. Venkatesudu,
# Symmetric Gauss Legendre quadrature formulas for composite numerical integration over
# a triangular surface,
# Applied Mathematics and Computation 188 (2007) 865–876
# <https://doi.org/10.1016/j.amc.2006.10.041>.
#
# Reported to Elsevier on June 24, 2019.
#

import numpy

from ..helpers import article
from ..line_segment import gauss_legendre
from ._helpers import TriangleScheme

citation = article(
    authors=["Frank G. Lether"],
    title="Computation of double integrals over a triangle",
    journal="Journal of Computational and Applied Mathematics",
    volume="2",
    number="3",
    month="sep",
    year="1976",
    pages="219–224",
    url="https://doi.org/10.1016/0771-050X(76)90008-5",
)


def lether(n):
    gl = gauss_legendre(n)

    w = numpy.outer((1 + gl.points) * gl.weights, gl.weights) / 4
    x = numpy.outer(1 - gl.points, numpy.ones(n)) / 2
    y = numpy.outer(1 + gl.points, 1 - gl.points) / 4

    points = numpy.array([x.flatten(), y.flatten()]).T
    weights = w.flatten()

    points = numpy.array([points[:, 0], points[:, 1], 1 - numpy.sum(points, axis=1)]).T

    degree = 2 * (n - 1)
    return TriangleScheme(f"Lether({n})", weights, points, degree, citation)

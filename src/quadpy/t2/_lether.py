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

import numpy as np

from ..c1 import gauss_legendre
from ..helpers import article
from ._helpers import T2Scheme, register

source = article(
    authors=["Frank G. Lether"],
    title="Computation of double integrals over a triangle",
    journal="Journal of Computational and Applied Mathematics",
    volume="2",
    number="3",
    month="sep",
    year="1976",
    pages="219–224",
    url="https://doi.org/10.1016/0771-050X%2876%2990008-5",
)


def lether(n):
    gl = gauss_legendre(n)

    w = np.outer((1 + gl.points) * gl.weights, gl.weights) / 4
    x = np.outer(1 - gl.points, np.ones(n)) / 2
    y = np.outer(1 + gl.points, 1 - gl.points) / 4

    points = np.array([x.flatten(), y.flatten()])
    weights = w.flatten()

    points = np.array([points[0], points[1], 1 - points[0] - points[1]])

    degree = 2 * (n - 1)
    return T2Scheme(
        f"Lether({n})", {"plain": np.vstack([weights, points])}, degree, source
    )


register([lether])

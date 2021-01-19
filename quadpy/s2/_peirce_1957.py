import numpy as np

from ..helpers import article
from ._helpers import S2Scheme, register

_source = article(
    authors=["W.H. Peirce"],
    title="Numerical integration over the planar annulus",
    journal="J. Soc. Indust. Appl. Math.",
    volume="5",
    number="2",
    month="jun",
    year="1957",
    url="https://www.jstor.org/stable/2098722",
)


# TODO symbolic
def peirce_1957(m):
    k = 4 * m + 3
    theta = 2 * np.pi * np.arange(1, k + 2) / (k + 1)
    p, w = np.polynomial.legendre.leggauss(m + 1)
    # scale points to [r0, r1] (where r0 = 0, r1 = 1 for now)
    p = np.sqrt(0.5 * (p + 1.0))
    p_theta = np.dstack(np.meshgrid(p, theta)).reshape(-1, 2).T
    points = np.array(
        [p_theta[0] * np.cos(p_theta[1]), p_theta[0] * np.sin(p_theta[1])]
    )
    # When integrating between 0 and 1, the weights are exactly the Gauss-Legendre
    # weights, scaled according to the disk area.
    weights = np.tile(0.5 / (k + 1) * w, k + 1)
    d = {"plain": [weights, points[0], points[1]]}
    return S2Scheme("Peirce 1957", d, k, _source)


register([peirce_1957])

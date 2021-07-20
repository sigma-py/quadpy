import numpy as np

from ..helpers import techreport
from ._helpers import C2Scheme, register

source = techreport(
    authors=["Shayne Waldron"],
    title="Symmetries of linear functionals",
    institution="University of Wisconsin-Madison, Center for Mathematical Sciences",
    note="CMS Technical Summary Report",
    month="oct",
    year="1994",
    url="https://www.semanticscholar.org/paper/Symmetries-of-Linear-Functionals-Waldron/086b084bac739c4616a06cc503dfcbda9f3cc7e8",
)


def waldron(r=1, alpha=0):
    assert r ** 2 >= 1 / 3

    R = r / np.sqrt(3 * r ** 2 - 1)

    beta = alpha + np.pi / 2

    sin_alpha = np.sin(alpha)
    cos_alpha = np.cos(alpha)
    sin_beta = np.sin(beta)
    cos_beta = np.cos(beta)

    d = {
        "plain": [
            [1 / 6 / r ** 2, 1 / 6 / r ** 2, 1 / 6 / R ** 2, 1 / 6 / R ** 2],
            [+r * cos_alpha, -r * cos_alpha, +R * cos_beta, -R * cos_beta],
            [+r * sin_alpha, -r * sin_alpha, +R * sin_beta, -R * sin_beta],
        ]
    }
    return C2Scheme("Waldron", d, 3, source)


register([waldron])

from sympy import Rational as frac
from sympy import sqrt

from ..helpers import article, untangle
from ._helpers import C3Scheme, fs_r00, fs_rr0, fs_rrs

source = article(
    authors=["Michael Sadowsky"],
    title="A Formula for Approximate Computation of a Triple Integral",
    journal="The American Mathematical Monthly",
    volume="47",
    number="8",
    month="oct",
    year="1940",
    pages="539-543",
    url="https://doi.org/10.2307/2303834",
)


def sadowsky():
    data = [
        (frac(91, 450), fs_r00(1)),
        (frac(-20, 225), fs_rr0(1)),
        (frac(8, 225), fs_rrs(sqrt(frac(5, 8)), 1)),
    ]
    points, weights = untangle(data)
    return C3Scheme("Sadowsky", weights, points, 5, source)

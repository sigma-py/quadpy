from sympy import Rational as frac
from sympy import sqrt

from ..helpers import article
from ._helpers import C3Scheme, register

source = article(
    authors=["Michael Sadowsky"],
    title="A Formula for Approximate Computation of a Triple Integral",
    journal="The American Mathematical Monthly",
    volume="47",
    number="8",
    month="oct",
    year="1940",
    pages="539-543",
    url="https://doi.org/10.1080/00029890.1940.11991016",
)


def sadowsky():
    d = {
        "symm_r00": [[frac(91, 450)], [1]],
        "symm_rr0": [[frac(-20, 225)], [1]],
        "symm_rrs": [[frac(8, 225)], [sqrt(frac(5, 8))], [1]],
    }
    return C3Scheme("Sadowsky", d, 5, source)


register([sadowsky])

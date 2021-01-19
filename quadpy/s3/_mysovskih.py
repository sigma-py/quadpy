import numpy as np
import sympy

from ..helpers import article
from ._helpers import S3Scheme, register

_source = article(
    authors=["I.P. Mysovskih"],
    title="On the construction of cubature formulas for the simplest regions",
    journal="Z. Vychisl. Mat. i. Mat. Fiz.",
    number="4",
    pages="3-14",
    year="1964",
)

frac = sympy.Rational
sqrt = np.vectorize(sympy.sqrt)
pi = sympy.pi
pm_ = np.array([+1, -1])


def mysovskih():
    sqrt17770 = sqrt(17770)
    r, s = sqrt((1715 - pm_ * 7 * sqrt17770) / 2817)
    t = sqrt(frac(7, 18))
    u = sqrt(frac(7, 27))

    B1, B2 = (2965 * sqrt17770 + pm_ * 227816) / 72030 / sqrt17770
    B3 = frac(324, 12005)
    B4 = frac(2187, 96040)

    d = {
        "symm_r00": [[B1, B2], [r, s]],
        "symm_rr0": [[B3], [t]],
        "symm_rrr": [[B4], [u]],
    }
    return S3Scheme("Mysovskih", d, 7, _source)


register([mysovskih])

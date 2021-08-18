# ENH closed forms of all schemes
from sympy import Rational as frac
from sympy import sqrt

from ..helpers import article
from ._helpers import C2Scheme, register

source = article(
    authors=["Preston C. Hammer", "Arthur H. Stroud"],
    title="Numerical Evaluation of Multiple Integrals II",
    journal="Math. Comp.",
    volume="12",
    year="1958",
    pages="272-280",
    url="https://doi.org/10.1090/S0025-5718-1958-0102176-6",
)


def hammer_stroud_1_2():
    d = {"d4_a0": [[frac(1, 4)], [sqrt(frac(2, 3))]]}
    return C2Scheme("Hammer-Stroud 1-2", d, 3, source)


def hammer_stroud_2_2():
    alpha = sqrt(frac(3, 5))
    d = {
        "zero2": [[frac(16, 81)]],
        "d4_a0": [[frac(10, 81)], [alpha]],
        "d4_aa": [[frac(25, 4 * 81)], [alpha]],
    }
    return C2Scheme("Hammer-Stroud 2-2", d, 5, source)


def hammer_stroud_3_2():
    xi1, xi2 = (sqrt(frac(3, 287) * (38 - i * sqrt(583))) for i in [+1, -1])
    d = {
        "d4_a0": [[frac(98, 4 * 405)], [sqrt(frac(6, 7))]],
        "d4_aa": [[0.5205929166673945 / 4, 0.2374317746906302 / 4], [xi1, xi2]],
    }
    return C2Scheme("Hammer-Stroud 3-2", d, 7, source, 4.441e-16)


register([hammer_stroud_1_2, hammer_stroud_2_2, hammer_stroud_3_2])

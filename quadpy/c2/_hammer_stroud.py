# ENH closed forms of all schemes
from sympy import Rational as frac
from sympy import sqrt

from ..helpers import article
from ._helpers import C2Scheme, expand_symmetries

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
    d = {
        "symm_r0": [[frac(1, 4)], [sqrt(frac(2, 3))]]
    }
    points, weights = expand_symmetries(d)
    return C2Scheme("Hammer-Stroud 1-2", weights, points, 3, source)


def hammer_stroud_2_2():
    alpha = sqrt(frac(3, 5))
    d = {
        "zero": [[frac(64, 81)]],
        "symm_r0": [[frac(40, 81)], [alpha]],
        "symm_s": [[frac(25, 81)], [alpha]]
    }
    points, weights = expand_symmetries(d)
    weights /= 4
    return C2Scheme("Hammer-Stroud 2-2", weights, points, 5, source)


def hammer_stroud_3_2():
    xi1, xi2 = [sqrt(frac(3, 287) * (38 - i * sqrt(583))) for i in [+1, -1]]
    d = {
        "symm_r0": [[frac(98, 405)], [sqrt(frac(6, 7))]],
        "symm_s": [[0.5205929166673945, 0.2374317746906302], [xi1, xi2]]
    }
    points, weights = expand_symmetries(d)
    weights /= 4
    return C2Scheme("Hammer-Stroud 3-2", weights, points, 7, source, 4.441e-16)

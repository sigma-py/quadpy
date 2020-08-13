from sympy import Rational as frac
from sympy import sqrt

from ..helpers import article
from ._helpers import C2Scheme, expand_symmetries

source = article(
    authors=["J.C. Maxwell"],
    title="On Approximate Multiple Integration between Limits by Summation",
    journal="Cambridge Library Collection - Physical Sciences",
    pages="604-611",
    # publisher="Cambridge University Press",
    url="https://doi.org/10.1017/CBO9780511710377.061",
    note="In W. Niven (Ed.), The Scientific Papers of James Clerk Maxwell. First published in 1890.",
)


def maxwell():
    r = sqrt(frac(12, 35))
    s, t = [sqrt((93 + i * 3 * sqrt(186)) / 155) for i in [+1, -1]]

    d = {
        "zero": [[frac(1, 81)]],
        "symm_r0": [[frac(49, 324)], [r]],
        # ERR typo in Stroud: 648 vs 649
        "symm_s_t": [[frac(31, 648)], [s], [t]],
    }
    points, weights = expand_symmetries(d)
    return C2Scheme("Maxwell", weights, points, 7, source)

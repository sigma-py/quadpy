from sympy import Rational as frac
from sympy import sqrt

from ..helpers import article
from ._helpers import C2Scheme, register

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
    s, t = (sqrt((93 + i * 3 * sqrt(186)) / 155) for i in [+1, -1])

    d = {
        "zero2": [[frac(1, 81)]],
        "d4_a0": [[frac(49, 324)], [r]],
        # ERR typo in Stroud: 648 vs 649
        "d4_ab": [[frac(31, 648)], [s], [t]],
    }
    return C2Scheme("Maxwell", d, 7, source)


register([maxwell])

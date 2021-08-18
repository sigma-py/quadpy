import warnings

from sympy import Rational as frac
from sympy import sqrt

from .. import cn
from ..helpers import book
from ._albrecht_collatz import albrecht_collatz_1 as stroud_c2_3_4
from ._albrecht_collatz import albrecht_collatz_2 as stroud_c2_5_1
from ._albrecht_collatz import albrecht_collatz_3 as stroud_c2_5_2
from ._albrecht_collatz import albrecht_collatz_4 as stroud_c2_5_6
from ._burnside import burnside as stroud_c2_5_3
from ._helpers import C2Scheme, register
from ._irwin import irwin_1 as stroud_c2_3_5
from ._irwin import irwin_2 as stroud_c2_5_7
from ._maxwell import maxwell as stroud_c2_7_3
from ._meister import meister as stroud_c2_7_6
from ._miller import miller as stroud_c2_1_2
from ._phillips import phillips as stroud_c2_7_2
from ._rabinowitz_richter import rabinowitz_richter_1 as stroud_c2_9_1
from ._rabinowitz_richter import rabinowitz_richter_2 as stroud_c2_11_1
from ._rabinowitz_richter import rabinowitz_richter_3 as stroud_c2_11_2
from ._rabinowitz_richter import rabinowitz_richter_4 as stroud_c2_13_1
from ._rabinowitz_richter import rabinowitz_richter_5 as stroud_c2_15_1
from ._rabinowitz_richter import rabinowitz_richter_6 as stroud_c2_15_2
from ._tyler import tyler_1 as stroud_c2_5_5
from ._tyler import tyler_2 as stroud_c2_7_1
from ._tyler import tyler_3 as stroud_c2_7_5

source = book(
    authors=["Arthur Stroud"],
    title="Approximate Calculation of Multiple Integrals",
    publisher="Prentice Hall",
    year="1971",
)


def stroud_c2_1_1():
    # product trapezoidal
    d = {"d4_aa": [[frac(1, 4)], [1]]}
    return C2Scheme("Stroud C2 1-1", d, 1, source)


def stroud_c2_3_1():
    # ERR misprint in Stroud: sqrt(1/3) vs 1/3
    d = {"d4_aa": [[frac(1, 4)], [sqrt(frac(1, 3))]]}
    return C2Scheme("Stroud C2 3-1", d, 3, source)


def stroud_c2_3_2():
    return cn.ewing(2)


def stroud_c2_3_3():
    return cn.stroud_cn_3_6(2)


def stroud_c2_5_4():
    r = sqrt(frac(3, 5))
    d = {
        "zero2": [[frac(16, 81)]],
        "d4_a0": [[frac(10, 81)], [r]],
        "d4_aa": [[frac(25, 324)], [r]],
    }
    return C2Scheme("Stroud C2 5-4", d, 5, source)


def stroud_c2_7_4():
    # product Gauss 7
    r, s = (sqrt((15 - i * 2 * sqrt(30)) / 35) for i in [+1, -1])

    B1, B2 = ((59 + i * 6 * sqrt(30)) / 864 for i in [+1, -1])
    B3 = frac(49, 864)

    r = sqrt(frac(3, 5))

    d = {"d4_aa": [[B1, B2], [r, s]], "d4_ab": [[B3], [r], [s]]}

    # TODO fix
    warnings.warn("Formula only has degree 1!")
    return C2Scheme("Stroud C2 7-4", d, 1, source)


register(
    [
        stroud_c2_1_1,
        stroud_c2_1_2,
        stroud_c2_3_1,
        stroud_c2_3_2,
        stroud_c2_3_3,
        stroud_c2_3_4,
        stroud_c2_3_5,
        stroud_c2_5_1,
        stroud_c2_5_2,
        stroud_c2_5_3,
        stroud_c2_5_4,
        stroud_c2_5_5,
        stroud_c2_5_6,
        stroud_c2_5_7,
        stroud_c2_7_1,
        stroud_c2_7_2,
        stroud_c2_7_3,
        stroud_c2_7_4,
        stroud_c2_7_5,
        stroud_c2_7_6,
        stroud_c2_9_1,
        stroud_c2_11_1,
        stroud_c2_11_2,
        stroud_c2_13_1,
        stroud_c2_15_1,
        stroud_c2_15_2,
    ]
)

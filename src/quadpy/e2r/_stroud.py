from sympy import Rational as frac
from sympy import sqrt

from ..helpers import book
from ._helpers import E2rScheme, register

# ERR misprint in Stroud copied from original article; rabinowitz_richter_4 as stroud_13_1,
from ._rabinowitz_richter import rabinowitz_richter_1 as stroud_9_1
from ._rabinowitz_richter import rabinowitz_richter_2 as stroud_11_1
from ._rabinowitz_richter import rabinowitz_richter_3 as stroud_11_2
from ._rabinowitz_richter import rabinowitz_richter_5 as stroud_15_1
from ._stroud_secrest import stroud_secrest_5 as stroud_5_1
from ._stroud_secrest import stroud_secrest_6 as stroud_7_1

_source = book(
    authors=["Arthur Stroud"],
    title="Approximate Calculation of Multiple Integrals",
    publisher="Prentice Hall",
    year="1971",
)


def stroud_4_1():
    d = {
        "zero2": [[frac(7, 10)]],
        "d5.0": [[frac(3, 50)], [[2 * sqrt(5)]]],
    }
    return E2rScheme("Stroud 4-1", d, 4, _source)


register(
    [
        stroud_4_1,
        stroud_5_1,
        stroud_7_1,
        stroud_9_1,
        stroud_11_1,
        stroud_11_2,
        stroud_15_1,
    ]
)

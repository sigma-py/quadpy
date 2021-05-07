import numpy as np
import sympy

from ..helpers import book, untangle, z
from ._albrecht import albrecht_4 as stroud_s2_9_1
from ._albrecht import albrecht_5 as stroud_s2_11_2
from ._albrecht import albrecht_6 as stroud_s2_13_2
from ._albrecht import albrecht_7 as stroud_s2_15_2
from ._albrecht import albrecht_8 as stroud_s2_17_1
from ._albrecht_collatz import albrecht_collatz as stroud_s2_3_2
from ._hammer_stroud import hammer_stroud_11_2 as stroud_s2_3_1
from ._hammer_stroud import hammer_stroud_12_2 as stroud_s2_5_2
from ._hammer_stroud import hammer_stroud_18 as stroud_s2_7_2
from ._helpers import S2Scheme, register
from ._mysovskih import mysovskih_1 as stroud_s2_4_1
from ._mysovskih import mysovskih_2 as stroud_s2_11_1
from ._mysovskih import mysovskih_3 as stroud_s2_15_1
from ._peirce_1956 import peirce_1956_1 as stroud_s2_7_1
from ._peirce_1956 import peirce_1956_2 as stroud_s2_9_5
from ._peirce_1956 import peirce_1956_3 as stroud_s2_11_4
from ._rabinowitz_richter import rabinowitz_richter_1 as stroud_s2_9_2
from ._rabinowitz_richter import rabinowitz_richter_2 as stroud_s2_9_4
from ._rabinowitz_richter import rabinowitz_richter_4 as stroud_s2_11_3
from ._rabinowitz_richter import rabinowitz_richter_5 as stroud_s2_13_1
from ._radon import radon

_source = book(
    authors=["Arthur Stroud"],
    title="Approximate Calculation of Multiple Integrals",
    publisher="Prentice Hall",
    year="1971",
)


def stroud_s2_5_1():
    return radon(0)


def stroud_s2_9_3():
    # spherical product gauss 9
    sqrt = np.vectorize(sympy.sqrt)
    pm_ = np.array([+1, -1])
    cos = np.vectorize(sympy.cos)
    sin = np.vectorize(sympy.sin)
    frac = sympy.Rational
    pi = sympy.pi

    r1, r2 = sqrt((6 - pm_ * sqrt(6)) / 10)

    a = 2 * (np.arange(10) + 1) * pi / 10
    x = np.array([cos(a), sin(a)]).T

    B0 = frac(1, 9)
    B1, B2 = (16 + pm_ * sqrt(6)) / 360

    data = [(B0, z(2)), (B1, r1 * x), (B2, r2 * x)]
    points, weights = untangle(data)
    points = np.ascontiguousarray(points.T)
    d = {"plain": [weights, points[0], points[1]]}
    return S2Scheme("Stroud S2 9-3", d, 9, _source)


register(
    [
        stroud_s2_3_1,
        stroud_s2_3_2,
        stroud_s2_4_1,
        stroud_s2_5_1,
        stroud_s2_5_2,
        stroud_s2_7_1,
        stroud_s2_7_2,
        stroud_s2_9_1,
        stroud_s2_9_2,
        stroud_s2_9_3,
        stroud_s2_9_4,
        stroud_s2_9_5,
        stroud_s2_11_1,
        stroud_s2_11_2,
        stroud_s2_11_3,
        stroud_s2_11_4,
        stroud_s2_13_1,
        stroud_s2_13_2,
        stroud_s2_15_1,
        stroud_s2_15_2,
        stroud_s2_17_1,
    ]
)

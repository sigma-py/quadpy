# -*- coding: utf-8 -*-
#
import numpy
import sympy

from ._helpers import _s, NCubeScheme
from ..helpers import untangle, article


_citation = article(
    authors=["Henry C. Thacher"],
    title="An efficient composite formula for multidimensional quadrature",
    journal="Communications of the ACM",
    volume="7",
    number="1",
    month="jan",
    year="1964",
    pages="23-25",
    url="https://doi.org/10.1145/363872.363897",
)


def thacher(n, symbolic=False):
    sqrt = sympy.sqrt if symbolic else numpy.sqrt

    r = sqrt(3) / 6
    data = [(1.0, [n * [2 * r]]), (+r, _s(n, -1, r)), (-r, _s(n, +1, r))]

    points, weights = untangle(data)
    reference_volume = 2 ** n
    weights *= reference_volume
    return NCubeScheme("Thacher", n, weights, points, 2, _citation)

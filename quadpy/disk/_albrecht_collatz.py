# -*- coding: utf-8 -*-
#
from __future__ import division

import numpy
import sympy

from ..helpers import untangle, pm, article
from ._helpers import DiskScheme


_citation = article(
    authors=["J. Albrecht", "L. Collatz"],
    title="Zur numerischen Auswertung mehrdimensionaler Integrale",
    journal="ZAMM",
    volume="38",
    number="1-2",
    year="1958",
    pages="1–15",
    url="https://doi.org/10.1002/zamm.19580380102",
)


def albrecht_collatz(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    pi = sympy.pi if symbolic else numpy.pi

    data = [
        # ERR Wrongly stated in Stroud as sqrt(1/2) instead of 1/2
        (frac(1, 4), pm(2, frac(1, 2)))
    ]
    points, weights = untangle(data)
    weights *= pi
    return DiskScheme("Albrecht-Collatz", weights, points, 3, _citation)

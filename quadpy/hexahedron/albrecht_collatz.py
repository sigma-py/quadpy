# -*- coding: utf-8 -*-
#
from __future__ import division

import sympy

from .helpers import fs_r00, fs_rr0, z, HexahedronScheme
from ..helpers import untangle, article

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

    data = [(frac(1, 4), z()), (frac(1, 12), fs_r00(1)), (frac(1, 48), fs_rr0(1))]

    points, weights = untangle(data)
    weights *= 8
    return HexahedronScheme("Albrecht-Collatz", weights, points, 3, _citation)

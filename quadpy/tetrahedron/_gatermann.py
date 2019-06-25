# -*- coding: utf-8 -*-
#
from ._helpers import untangle2, TetrahedronScheme
from ..helpers import article

citation = article(
    authors=["Karin Gatermann"],
    title="Linear Representations of Finite Groups and The Ideal Theoretical Construction of G-Invariant Cubature Formulas",
    journal="Numerical Integration",
    pages="25-35",
    note="Part of the NATO ASI Series book series (ASIC, volume 357)",
)


def gatermann():
    degree = 5
    data = {
        "s31": [
            [
                9.73033316198362119165356216965707e-06,
                0.656936552995394536166881327385593,
            ],
            [
                8.99031481668747219698547129902142e-03,
                0.0801424420792727848879183805550907,
            ],
        ],
        "s22": [
            [
                2.17777476778781405656596945369837e-02,
                0.404475329343454044779549906725159,
            ]
        ],
    }
    bary, weights = untangle2(data)
    weights *= 6
    return TetrahedronScheme("Gatermann", weights, bary, degree, citation)

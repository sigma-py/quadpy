# -*- coding: utf-8 -*-
#
from ._helpers import QuadrilateralScheme, concat, pm, pm2
from ..helpers import article

citation = article(
    authors=["Ann Haegemans", "Robert Piessens"],
    title="Construction of Cubature Formulas of Degree Seven and Nine Symmetric Planar Regions Using Orthogonal Polynomials",
    journal="SIAM Journal on Numerical Analysis",
    volume="14",
    number="3",
    month="jun",
    year="1977",
    pages="492-508",
    url="https://www.jstor.org/stable/2156699",
)


def haegemans_piessens():
    weights, points = concat(
        pm2(
            [0.213057211620949126, 0.9171178223127705862, 0.547931206828092323],
            [0.17400948894689560610, 0.61126876646532841440, 0.93884325665885830459],
        ),
        pm(
            [0.63585388344327977182, +0.52942280204265532589, 0.0],
            [0.59001271542103076297, 0.0, +0.62704137378039531763],
        ),
    )
    return QuadrilateralScheme("Haegemans-Piessens", weights, points, 7, citation)

# -*- coding: utf-8 -*-
#
from ..helpers import article
from ._helpers import TetrahedronScheme, untangle2

citation = article(
    authors=["D.M. Williams", "L. Shunn", "A. Jameson"],
    title="Symmetric quadrature rules for simplexes based on sphere close packed lattice arrangements",
    journal="Journal of Computational and Applied Mathematics",
    volume="266",
    year="2014",
    pages="18–38",
    url="https://doi.org/10.1016/j.cam.2014.01.007",
)


def williams_shunn_jameson():
    degree = 9
    data = {
        "s31": [
            [0.002144935144316, 0.026878474414817],
            [0.020826641690769, 0.187140675803470],
            [0.023000681669286, 0.322111431830857],
        ],
        "s22": [
            [0.007210136064455, 0.473575835127937],
            [0.030798919159712, 0.352045262027356],
        ],
        "s211": [
            [0.004357844813864, 0.020953442220056, 0.225783205866940],
            [0.008593530677833, 0.096989733123466, 0.158462939666092],
            [0.004863063904912, 0.097608162890442, 0.011844417749498],
        ],
        "s1111": [
            [0.015595140078259, 0.541184412800237, 0.133558160703568, 0.296501020543124]
        ],
    }
    points, weights = untangle2(data)
    return TetrahedronScheme(
        "Williams-Shunn-Jameson", weights, points, degree, citation
    )

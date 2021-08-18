from sympy import Rational as frac
from sympy import sqrt

from ..helpers import article
from ._classical import seven_point as lyness_jespersen_04
from ._helpers import T2Scheme, register
from ._strang_fix_cowper import strang_fix_cowper_09 as lyness_jespersen_10

source = article(
    authors=["J.N. Lyness", "D. Jespersen"],
    title="Moderate Degree Symmetric Quadrature Rules for the Triangle",
    journal="IMA Journal of Applied Mathematics",
    year="1975",
    month="feb",
    volume="15",
    number="1",
    pages="19-32",
    url="https://doi.org/10.1093/imamat/15.1.19",
)


def lyness_jespersen_01():
    d = {"d3_aa": [[frac(1, 3)], [frac(1, 2)]]}
    return T2Scheme("Lyness-Jespersen 1", d, 2, source)


def lyness_jespersen_02():
    d = {
        "centroid": [[frac(3, 4)]],
        "vertex": [[frac(1, 12)]],
    }
    return T2Scheme("Lyness-Jespersen 2", d, 2, source)


def lyness_jespersen_03():
    d = {
        "centroid": [[-frac(9, 16)]],
        "d3_aa": [[frac(25, 48)], [frac(1, 5)]],
    }
    return T2Scheme("Lyness-Jespersen 3", d, 3, source)


def lyness_jespersen_05():
    d = {
        "d3_aa": [
            [3.298552309659655e-01 / 3, 6.701447690340345e-01 / 3],
            [9.157621350977073e-02, 4.459484909159649e-01],
        ]
    }
    return T2Scheme("Lyness-Jespersen 5", d, 4, source)


def lyness_jespersen_06():
    a0, a1 = ((3 + i * sqrt(3)) / 6 for i in [+1, -1])
    d = {
        "centroid": [frac(9, 20)],
        "vertex": [[-frac(1, 60)]],
        "d3_ab": [[frac(1, 10)], [a0], [a1]],
    }
    return T2Scheme("Lyness-Jespersen 6", d, 4, source)


def lyness_jespersen_07():
    sqrt13 = sqrt(13)
    d = {
        "vertex": [[(11 - sqrt13) / 360]],
        "d3_aa": [
            [(10 - 2 * sqrt13) / 45, (29 + 17 * sqrt13) / 360],
            [frac(1, 2), (7 - sqrt13) / 18],
        ],
    }
    return T2Scheme("Lyness-Jespersen 7", d, 4, source)


def lyness_jespersen_08():
    sqrt15 = sqrt(15)
    b1, b2 = ((155 - i * sqrt15) / 1200 for i in [+1, -1])
    r1, r2 = ((6 - i * sqrt15) / 21 for i in [+1, -1])

    d = {"centroid": [[frac(9, 40)]], "d3_aa": [[b1, b2], [r1, r2]]}
    return T2Scheme("Lyness-Jespersen 8", d, 5, source)


def lyness_jespersen_09():
    d = {
        "centroid": [frac(81, 320)],
        "vertex": [[frac(1, 90)]],
        "d3_aa": [[frac(16, 225), frac(2401, 14400)], [frac(1, 2), frac(1, 7)]],
    }
    return T2Scheme("Lyness-Jespersen 9", d, 5, source)


def lyness_jespersen_11():
    c, d = ((3 + i * sqrt(6)) / 6 for i in [+1, -1])
    d = {
        "centroid": [[-frac(81, 140)]],
        "vertex": [[-frac(5, 252)]],
        "d3_aa": [[frac(17, 315), frac(128, 315)], [frac(1, 2), frac(1, 4)]],
        "d3_ab": [[frac(9, 210)], [c], [d]],
    }
    return T2Scheme("Lyness-Jespersen 11", d, 6, source)


def lyness_jespersen_12():
    d = {
        "centroid": [[1.527089667883523e-01]],
        "d3_aa": [
            [2.944076042366762e-01 / 3, 3.887052878418766e-01 / 3],
            [4.738308139536513e-01, 1.721176696308175e-01],
        ],
        "d3_ab": [[1.641781411330949e-01 / 6], [0], [8.653073540834571e-01]],
    }
    return T2Scheme("Lyness-Jespersen 12", d, 6, source, 2.628e-15)


def lyness_jespersen_13():
    d = {
        "centroid": [[-1.495700444677495e-01]],
        "d3_aa": [
            [5.268457722996328e-01 / 3, 1.600417068265167e-01 / 3],
            [2.603459660790466e-01, 6.513010290221623e-02],
        ],
        "d3_ab": [
            [4.626825653415500e-01 / 6],
            [6.384441885698096e-01],
            [4.869031542531756e-02],
        ],
    }
    return T2Scheme("Lyness-Jespersen 13", d, 7, source, 4.119e-14)


def lyness_jespersen_14():
    d = {
        "centroid": [[1.763126156005252e-01]],
        "vertex": [[1.210901532763310e-02 / 3]],
        "d3_aa": [
            [3.499561757697094e-01 / 3, 3.195119754425220e-01 / 3],
            [1.549360602237604e-01, 4.691507461438120e-01],
        ],
        "d3_ab": [[1.421102178595603e-01 / 6], [0], [8.392991722729236e-01]],
    }
    return T2Scheme("Lyness-Jespersen 14", d, 7, source, 2.805e-13)


def lyness_jespersen_15():
    d = {
        "centroid": [[1.443156076777862e-01]],
        "d3_aa": [
            [
                2.852749028018549e-01 / 3,
                9.737549286959440e-02 / 3,
                3.096521116041552e-01 / 3,
            ],
            [4.592925882927229e-01, 5.054722831703103e-02, 1.705693077517601e-01],
        ],
        "d3_ab": [
            [1.633818850466092e-01 / 6],
            [8.394777409957211e-03],
            [7.284923929554041e-01],
        ],
    }
    return T2Scheme("Lyness-Jespersen 15", d, 8, source, 8.854e-15)


def lyness_jespersen_16():
    d = {
        "vertex": [[+1.207273935292775e-02 / 3]],
        "d3_aa": [
            [
                -8.491579879151455e-01 / 3,
                +1.042367468891334e00 / 3,
                +1.947229791412260e-01 / 3,
                +4.511852767201322e-01 / 3,
            ],
            [0.5, 4.956813941755582e-01, 9.032775751426533e-02, 2.341547497073052e-01],
        ],
        "d3_ab": [[+1.488095238055238e-01 / 6], [0], [7.236067977499750e-01]],
    }
    return T2Scheme("Lyness-Jespersen 16", d, 8, source, 7.282e-12)


def lyness_jespersen_17():
    d = {
        "centroid": [[-2.834183851113958e-01]],
        "d3_aa": [
            [
                2.097208857979572e-01 / 3,
                5.127273801480265e-02 / 3,
                6.564896469913508e-01 / 3,
            ],
            [4.766654393821525e-01, 3.377184405448033e-02, 2.703478891654040e-01],
        ],
        "d3_ab": [
            [3.659351143072855e-01 / 6],
            [5.146433548666149e-02],
            [7.458294907672514e-01],
        ],
    }
    return T2Scheme("Lyness-Jespersen 17", d, 8, source)


def lyness_jespersen_18():
    d = {
        "centroid": [[9.713579628279610e-02]],
        "d3_aa": [
            [
                9.400410068141950e-02 / 3,
                2.334826230143263e-01 / 3,
                2.389432167816271e-01 / 3,
                7.673302697609430e-02 / 3,
            ],
            [
                4.896825191987370e-01,
                4.370895914929355e-01,
                1.882035356190322e-01,
                4.472951339445297e-02,
            ],
        ],
        "d3_ab": [
            [2.597012362637364e-01 / 6],
            [3.683841205473626e-02],
            [7.411985987844980e-01],
        ],
    }
    return T2Scheme("Lyness-Jespersen 18", d, 9, source, 8.604e-15)


def lyness_jespersen_19():
    d = {
        "centroid": [[1.133624844599192e-01]],
        "vertex": [[1.062573789846330e-03 / 3]],
        "d3_aa": [
            [
                4.803411513859279e-02 / 3,
                2.524243006337300e-01 / 3,
                7.819254371487040e-02 / 3,
                2.472227459993048e-01 / 3,
            ],
            [0.5, 4.497793381870162e-01, 4.694744319909033e-02, 1.918719127374489e-01],
        ],
        "d3_ab": [
            [2.597012362637364e-01 / 6],
            [3.683841205473626e-02],
            [7.411985987844980e-01],
        ],
    }
    return T2Scheme("Lyness-Jespersen 19", d, 9, source, 3.733e-15)


def lyness_jespersen_20():
    d = {
        "d3_aa": [
            [
                4.097919300803106e-02 / 3,
                1.085536215102866e-01 / 3,
                2.781018986881812e-03 / 3,
                1.779689321422668e-01 / 3,
                2.314486047444677e-01 / 3,
            ],
            [
                3.236494811127173e-02,
                1.193509122825931e-01,
                5.346110482707572e-01,
                2.033099004312816e-01,
                3.989693029658558e-01,
            ],
        ],
        "d3_ab": [
            [3.140226717732234e-01 / 6, 1.242459578348437e-01 / 6],
            [5.017813831049474e-02, 2.102201653616613e-02],
            [5.932012134282132e-01, 8.074890031597923e-01],
        ],
    }
    return T2Scheme("Lyness-Jespersen 20", d, 11, source, 1.361e-13)


def lyness_jespersen_21():
    d = {
        "centroid": [[8.797730116222190e-02]],
        "d3_aa": [
            [
                2.623293466120857e-02 / 3,
                1.142447159818060e-01 / 3,
                5.656634416839376e-02 / 3,
                2.164790926342230e-01 / 3,
                2.079874161166116e-01 / 3,
            ],
            [
                2.598914092828833e-02,
                9.428750264792270e-02,
                4.946367750172147e-01,
                2.073433826145142e-01,
                4.389078057004907e-01,
            ],
        ],
        "d3_ab": [
            [4.417430269980344e-02 / 6, 2.463378925757316e-01 / 6],
            [0, 4.484167758913055e-02],
            [8.588702812826364e-01, 6.779376548825902e-01],
        ],
    }
    return T2Scheme("Lyness-Jespersen 21", d, 11, source, 3.543e-14)


triex_19 = lyness_jespersen_18
triex_28 = lyness_jespersen_21


register(
    [
        lyness_jespersen_01,
        lyness_jespersen_02,
        lyness_jespersen_03,
        lyness_jespersen_04,
        lyness_jespersen_05,
        lyness_jespersen_06,
        lyness_jespersen_07,
        lyness_jespersen_08,
        lyness_jespersen_09,
        lyness_jespersen_10,
        lyness_jespersen_11,
        lyness_jespersen_12,
        lyness_jespersen_13,
        lyness_jespersen_14,
        lyness_jespersen_15,
        lyness_jespersen_16,
        lyness_jespersen_17,
        lyness_jespersen_18,
        lyness_jespersen_19,
        lyness_jespersen_20,
        lyness_jespersen_21,
    ]
)

from sympy import Rational as frac

from ..helpers import article
from ._helpers import T2Scheme, concat, s1, s2, s3

source = article(
    authors=["D.M. Williams", "L. Shunn", "A. Jameson"],
    title="Symmetric quadrature rules for simplexes based on sphere close packed lattice arrangements",
    journal="Journal of Computational and Applied Mathematics",
    volume="266",
    year="2014",
    pages="18–38",
    url="https://doi.org/10.1016/j.cam.2014.01.007",
)


def williams_shunn_jameson_1():
    weights, points = s3(1)
    return T2Scheme("Williams-Shunn-Jameson 1", weights, points, 1, source)


def williams_shunn_jameson_2():
    weights, points = s2([frac(1, 3), frac(1, 6)])
    return T2Scheme("Williams-Shunn-Jameson 2", weights, points, 2, source)


def williams_shunn_jameson_3():
    weights, points = s2(
        [0.109951743655333, 0.091576213509780], [0.223381589678000, 0.445948490915964]
    )
    return T2Scheme("Williams-Shunn-Jameson 3", weights, points, 4, source)


def williams_shunn_jameson_4():
    weights, points = concat(
        s3(0.201542988584730),
        s2([0.041955512996649, 0.055564052669793]),
        s1([0.112098412070887, 0.295533711735893, 0.634210747745723]),
    )
    return T2Scheme("Williams-Shunn-Jameson 4", weights, points, 5, source)


def williams_shunn_jameson_5():
    weights, points = concat(
        s2(
            [0.017915455012303, 0.035870877695734],
            [0.127712195881265, 0.241729395767967],
            [0.076206062385535, 0.474308787777079],
        ),
        s1([0.055749810027115, 0.201503881881800, 0.751183631106484]),
    )
    return T2Scheme("Williams-Shunn-Jameson 5", weights, points, 7, source)


def williams_shunn_jameson_6():
    weights, points = concat(
        s2(
            [0.010359374696538, 0.028112952182664],
            [0.075394884326738, 0.177139098469317],
            [0.097547802373242, 0.405508595867433],
        ),
        s1(
            [0.028969269372473, 0.148565812270887, 0.817900980028499],
            [0.046046366595935, 0.357196298615681, 0.604978911775132],
        ),
    )
    return T2Scheme("Williams-Shunn-Jameson 6", weights, points, 8, source)


def williams_shunn_jameson_7():
    weights, points = concat(
        s3(0.083608212215637),
        s2(
            [0.005272170280495, 0.019977187122193],
            [0.044552936679504, 0.131721767529998],
            [0.033815712804198, 0.485135346793461],
        ),
        s1(
            [0.015710461340183, 0.107951981846011, 0.867911210117951],
            [0.028205136280616, 0.270840772921567, 0.700872570380723],
            [0.066995957127830, 0.316549598844617, 0.536654684206138],
        ),
    )
    return T2Scheme("Williams-Shunn-Jameson 7", weights, points, 10, source)


def williams_shunn_jameson_8():
    weights, points = concat(
        s2(
            [0.005639123786910, 0.021171422779465],
            [0.027148968192278, 0.100584397395888],
            [0.063100912533359, 0.271038307711932],
            [0.051752795679899, 0.440191258403832],
        ),
        s1(
            [0.009866753574646, 0.101763679498021, 0.879979641427232],
            [0.022008204800147, 0.394033271669987, 0.582562022863673],
            [0.016644570076736, 0.226245530909229, 0.751530614542782],
            [0.044326238118914, 0.635737183263105, 0.249079227621332],
        ),
    )
    return T2Scheme("Williams-Shunn-Jameson 8", weights, points, 12, source)

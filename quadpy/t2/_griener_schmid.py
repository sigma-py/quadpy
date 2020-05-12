import numpy

from ..helpers import article
from ._helpers import T2Scheme

citation = article(
    authors=["Bernhard Griener", "Hans Joachim Schmid"],
    title="An interactive tool to visualize common zeros of two-dimensional polynomials",
    journal="Journal of Computational and Applied Mathematics",
    volume="112",
    year="1999",
    pages="83-94",
    url="https://doi.org/10.1016/S0377-0427(99)00215-0",
)

# c2 = article(
#     authors="G.G. Rasputin",
#     title="Construction of cubature formulas containing prespecied knots",
#     journal="Metody Vychisl.",
#     volume="13",
#     year="1983",
#     pages="122–129",
#     note="in Russian."
# )


def griener_schmid_1():
    # According to the article, this scheme appeared earlier in `c2`.
    points = numpy.array(
        [
            [1.28867990757340236072, -0.00879199714420631034],
            [-0.00879199714420631034, 1.28867990757340236072],
            [0.78822920956686404839, 0.06468880605601090940],
            [0.06468880605601090940, 0.78822920956686404839],
            [0.40172877323475981682, 0.07682571684908210504],
            [0.07682571684908210504, 0.40172877323475981682],
            [0.64933214716985020413, 0.29502511936941515372],
            [0.29502511936941515372, 0.64933214716985020413],
            [0.08533316117031069729, 0.08533316117031069729],
            [0.35369054666996447962, 0.35369054666996447962],
        ]
    )
    weights = numpy.array(
        [
            0.00051542875928455448,
            0.00051542875928455448,
            0.10479531313284680990,
            0.10479531313284680990,
            0.14855110930104331913,
            0.14855110930104331913,
            0.09497723917756742321,
            0.09497723917756742321,
            0.09099352359044946853,
            0.21132829566806631803,
        ]
    )
    points = numpy.array([points[:, 0], points[:, 1], 1 - numpy.sum(points, axis=1)]).T
    return T2Scheme("Griener-Schmid 1", weights, points, 6, citation)


def griener_schmid_2():
    points = numpy.array(
        [
            [0.62261842170067743793, 0.29681206443213099158],
            [0.29681206443213099158, 0.62261842170067743793],
            [0.86560137104657306697, 0.05848274947077690073],
            [0.05848274947077690073, 0.86560137104657306697],
            [0.30783381768828903704, 0.30783381768828903704],
            [0.56391552905457826463, 0.06887470996336075494],
            [0.06887470996336075494, 0.56391552905457826463],
            [-0.17255363572833024944, -0.17255363572833024944],
            [0.19794325221245567307, 0.04812894161425963489],
            [0.04812894161425963481, 0.19794325221245567307],
        ]
    )
    weights = numpy.array(
        [
            0.12493284637910931301,
            0.12493284637910931301,
            0.05736898913621279825,
            0.05736898913621279825,
            0.21037615383745954659,
            0.12547963681681714946,
            0.12547963681681714946,
            0.00042783665681296109,
            0.08681653242072448544,
            0.08681653242072448544,
        ]
    )

    points = numpy.array([points[:, 0], points[:, 1], 1 - numpy.sum(points, axis=1)]).T
    return T2Scheme("Griener-Schmid 2", weights, points, 6, citation)

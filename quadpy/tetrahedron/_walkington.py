# -*- coding: utf-8 -*-
#
import numpy

from ..helpers import techreport
from ._helpers import TetrahedronScheme

citation = techreport(
    authors=["Noel J. Walkington"],
    title="Quadrature on simplices of arbitrary dimension",
    institution="CMU",
    year="2000",
    url="https://www.math.cmu.edu/~nw0z/publications/00-CNA-023/023abs/",
)


def walkington_p5():
    degree = 5
    weights = 6 * numpy.concatenate(
        [
            numpy.full(4, 0.018781320953002641800),
            numpy.full(4, 0.012248840519393658257),
            numpy.full(6, 0.0070910034628469110730),
        ]
    )
    points = numpy.concatenate(
        [
            _xi1(0.31088591926330060980),
            _xi1(0.092735250310891226402),
            _xi11(0.045503704125649649492),
        ]
    )
    return TetrahedronScheme("Walkington p5", weights, points, degree, citation)


def _xi1(a):
    b = 1 - 3 * a
    return numpy.array([[b, a, a, a], [a, b, a, a], [a, a, b, a], [a, a, a, b]])


def _xi11(a):
    b = (1 - 2 * a) / 2
    return numpy.array(
        [
            [b, b, a, a],
            [b, a, b, a],
            [b, a, a, b],
            [a, b, a, b],
            [a, a, b, b],
            [a, b, b, a],
        ]
    )

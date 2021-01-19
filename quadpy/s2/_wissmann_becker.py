import math

import numpy as np

from ..helpers import article, untangle
from ._helpers import S2Scheme, register

_source = article(
    authors=["Johannes W. Wissmann", "Thomas Becker"],
    title="Partially Symmetric Cubature Formulas for Even Degrees of Exactness",
    journal="SIAM J. Numer. Anal.",
    volume="23",
    number="3",
    year="1986",
    pages="676–685",
    url="https://doi.org/10.1137/0723043",
)


def wissmann_becker_6_1():
    data = [
        (0.192106330409167, _z(0.914874261544932)),
        (0.642883172453235, _z(-0.206355543412215)),
        (0.030327815305342, _m(0.978839047064283, +0.662952453661044)),
        (0.502358408812314, _m(0.458236515708684, +0.464889889783278)),
        (0.340540094500688, _m(0.802017004072343, -0.174207378839556)),
        (0.280075256745352, _m(0.373683864304499, -0.770749853148807)),
    ]
    points, weights = untangle(data)
    weights /= math.pi
    points = np.ascontiguousarray(points.T)
    d = {"plain": [weights, points[0], points[1]]}
    return S2Scheme("Wissmann-Becker 6-1", d, 6, _source, 1.108e-14)


def wissmann_becker_6_2():
    data = [
        (0.552687596464871, _z(0.0)),
        (0.076090439302545, _z(1.032546063890842)),
        (0.441690572122440, _z(-0.726359846042944)),
        (0.076090439302545, _m(0.982009662438297, +0.319074281217230)),
        (0.076090439302545, _m(0.606915348667678, -0.835347313162651)),
        (0.441690572122440, _m(0.426943605361474, +0.587637459480312)),
        (0.441690572122440, _m(0.690809264754287, -0.224457536458840)),
    ]
    points, weights = untangle(data)
    weights /= math.pi
    points = np.ascontiguousarray(points.T)
    d = {"plain": [weights, points[0], points[1]]}
    return S2Scheme("Wissmann-Becker 6-2", d, 6, _source)


def wissmann_becker_8_1():
    data = [
        (0.334521439965580, _z(0.000000000000000)),
        (0.087911387704639, _z(0.953321175521807)),
        (0.167628039106469, _z(-0.882458098900126)),
        (0.305874815913735, _z(0.582334188120499)),
        (0.087911387704639, _m(0.906662316102171, +0.294592444333741)),
        (0.087911387704639, _m(0.560348127669843, -0.771253032094644)),
        (0.167628039106469, _m(0.518695856299547, +0.713923598834010)),
        (0.167628039106469, _m(0.839267525316398, -0.272694549383947)),
        (0.305874815913735, _m(0.553832724273449, +0.179951160534772)),
        (0.305874815913735, _m(0.342287447682940, -0.471118254595022)),
    ]
    points, weights = untangle(data)
    weights /= math.pi
    points = np.ascontiguousarray(points.T)
    d = {"plain": [weights, points[0], points[1]]}
    return S2Scheme("Wissmann-Becker 8-1", d, 8, _source)


def _z(a):
    return np.array([[0.0, a]])


def _m(a, b):
    return np.array([[+a, +b], [-a, +b]])


register([wissmann_becker_6_1, wissmann_becker_6_2, wissmann_becker_8_1])

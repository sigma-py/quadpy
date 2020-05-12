import numpy
from sympy import Rational as frac
from sympy import sqrt

from ..helpers import article
from ._helpers import QuadrilateralScheme

citation = article(
    authors=["H.J. Schmid"],
    title="On cubature formulae with a minimal number of knots",
    journal="Numerische Mathematik",
    month="sep",
    year="1978",
    volume="31",
    number="3",
    pages="281–297",
    url="https://eudml.org/doc/132580",
)


def schmid_2():
    points = numpy.array(
        [
            [-sqrt(frac(1, 3)), +sqrt(frac(2, 3))],
            [-sqrt(frac(1, 3)), -sqrt(frac(2, 3))],
            [+sqrt(frac(1, 3)), 0],
        ]
    )
    weights = numpy.array([frac(1, 4), frac(1, 4), frac(1, 2)])
    weights *= 4
    return QuadrilateralScheme("Schmid 2", weights, points, 2, citation)


def schmid_4():
    points = numpy.array(
        [
            [0, (sqrt(3) + sqrt(15)) / 6],
            [0, (sqrt(3) - sqrt(15)) / 6],
            [+sqrt(15) / 5, (+sqrt(87) - 2 * sqrt(3)) / 15],
            [-sqrt(15) / 5, (+sqrt(87) - 2 * sqrt(3)) / 15],
            [+sqrt(15) / 5, (-sqrt(87) - 2 * sqrt(3)) / 15],
            [-sqrt(15) / 5, (-sqrt(87) - 2 * sqrt(3)) / 15],
        ]
    )
    weights = numpy.array(
        [
            frac(2, 9) - 2 * sqrt(5) / 45,
            frac(2, 9) + 2 * sqrt(5) / 45,
            frac(5, 36) + 5 * sqrt(29) / 18 / 29,
            frac(5, 36) + 5 * sqrt(29) / 18 / 29,
            frac(5, 36) - 5 * sqrt(29) / 18 / 29,
            frac(5, 36) - 5 * sqrt(29) / 18 / 29,
        ]
    )
    weights *= 4
    return QuadrilateralScheme("Schmid 4", weights, points, 4, citation)


def schmid_6():
    # TODO better-quality points/weights for Schmidt
    points = numpy.array(
        [
            [+0.000000000000, +0.774596669241],
            [+0.563604836881, -0.795508520349],
            [+0.838331011044, +0.845091361153],
            [+0.651030930900, +0.166755021097],
            [-0.484792881050, -0.927694708202],
            [-0.914603935097, -0.520771886130],
            [-0.135220856964, -0.279191827433],
            [-0.731697727745, +0.417391901524],
            [-0.887824220291, +1.075479856096],
            [+1.101172842321, -0.485302501018],
        ]
    )
    weights = numpy.array(
        [
            0.140845070423,
            0.113931725656,
            0.049023075184,
            0.168918151204,
            0.063463914536,
            0.066611011696,
            0.214897708035,
            0.145149421990,
            0.014704280797,
            0.022455640481,
        ]
    )
    weights *= 4
    return QuadrilateralScheme("Schmid 6", weights, points, 6, citation)

# -*- coding: utf-8 -*-
#
from __future__ import division

import numpy
import sympy


class Schmid(object):
    """
    H.J. Schmid,
    On cubature formulae with a minimal number of knots,
    Numerische Mathematik, September 1978, Volume 31, Issue 3, pp 281–297,
    <https://eudml.org/doc/132580>.

    Abstract:
    In this paper an approach is outlined to the two-dimensional analogon of the
    Gaussian quadrature problem. The main results are necessary and sufficient
    conditions for the existence of cubature formulae which are exact for all
    polynomials of degree ≦m and which have a minimal number of 1/2k(k+1)
    knots,k=[m/2]+1. Ifm is odd, similar results are due to I.P. Mysovskikh ([5, 6])
    which will be derived in a new way as a special case of the general characterization
    given here. Furthermore, it will be shown how this characterization can be used to
    construct minimal formulae of even degree.
    """

    def __init__(self, index, symbolic=False):
        frac = sympy.Rational if symbolic else lambda x, y: x / y
        sqrt = sympy.sqrt if symbolic else numpy.sqrt

        self.name = "Schmid({})".format(index)

        if index == 2:
            self.degree = 2
            self.points = numpy.array(
                [
                    [-sqrt(frac(1, 3)), +sqrt(frac(2, 3))],
                    [-sqrt(frac(1, 3)), -sqrt(frac(2, 3))],
                    [+sqrt(frac(1, 3)), 0],
                ]
            )
            self.weights = numpy.array([frac(1, 4), frac(1, 4), frac(1, 2)])
        elif index == 4:
            self.degree = 4
            self.points = numpy.array(
                [
                    [0, (sqrt(3) + sqrt(15)) / 6],
                    [0, (sqrt(3) - sqrt(15)) / 6],
                    [+sqrt(15) / 5, (+sqrt(87) - 2 * sqrt(3)) / 15],
                    [-sqrt(15) / 5, (+sqrt(87) - 2 * sqrt(3)) / 15],
                    [+sqrt(15) / 5, (-sqrt(87) - 2 * sqrt(3)) / 15],
                    [-sqrt(15) / 5, (-sqrt(87) - 2 * sqrt(3)) / 15],
                ]
            )
            self.weights = numpy.array(
                [
                    frac(2, 9) - 2 * sqrt(5) / 45,
                    frac(2, 9) + 2 * sqrt(5) / 45,
                    frac(5, 36) + 5 * sqrt(29) / 18 / 29,
                    frac(5, 36) + 5 * sqrt(29) / 18 / 29,
                    frac(5, 36) - 5 * sqrt(29) / 18 / 29,
                    frac(5, 36) - 5 * sqrt(29) / 18 / 29,
                ]
            )
        else:
            # TODO get more digits
            assert index == 6
            self.degree = 6
            self.points = numpy.array(
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
            self.weights = numpy.array(
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

        self.weights *= 4
        return

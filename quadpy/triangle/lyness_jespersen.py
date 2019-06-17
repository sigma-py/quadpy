# -*- coding: utf-8 -*-
#
"""
J.N. Lyness, D. Jespersen,
Moderate Degree Symmetric Quadrature Rules for the Triangle,
J. Inst. Maths Applies (1975) 15, 19-32,
<https://doi.org/10.1093/imamat/15.1.19>.

Abstract:
A variant formulation of the moment fitting equations for the construction of D3
(triangularly symmetric) quadrature rules for the triangle is derived. These equations
are solved to produce weights and abscissas for quadrature rules of polynomial degree up
to 11 for the triangle, some of which require fewer function evaluations than any
presently available rule of the same polynomial degree. Cytolic rules of degrees up to 9
are also derived.
"""
from __future__ import division

import numpy
import sympy

from .helpers import TriangleScheme, s2, concat, s3, s1


def lyness_jespersen_1(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y

    weights, bary = s2([frac(1, 3), frac(1, 2)])
    return TriangleScheme("Lyness-Jespersen 1", 2, weights, bary)


def lyness_jespersen_2(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y

    weights, bary = concat(s3(frac(3, 4)), s2([frac(1, 12), 0]))
    return TriangleScheme("Lyness-Jespersen 2", 2, weights, bary)


def lyness_jespersen_3(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y

    weights, bary = concat(s3(-frac(9, 16)), s2([frac(25, 48), frac(1, 5)]))
    return TriangleScheme("Lyness-Jespersen 3", 3, weights, bary)


def lyness_jespersen_4(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y

    weights, bary = concat(
        s3(frac(9, 20)), s2([frac(1, 20), 0], [frac(2, 15), frac(1, 2)])
    )
    return TriangleScheme("Lyness-Jespersen 4", 3, weights, bary)


def lyness_jespersen_5():
    weights, bary = s2(
        [3.298552309659655e-01 / 3, 9.157621350977073e-02],
        [6.701447690340345e-01 / 3, 4.459484909159649e-01],
    )
    return TriangleScheme("Lyness-Jespersen 5", 4, weights, bary)


def lyness_jespersen_6(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt

    a0, a1 = [(3 + i * sqrt(3)) / 6 for i in [+1, -1]]

    weights, bary = concat(
        s3(+frac(9, 20)), s2([-frac(1, 60), 0]), s1([+frac(1, 10), a0, a1])
    )
    return TriangleScheme("Lyness-Jespersen 6", 4, weights, bary)


def lyness_jespersen_7(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt

    sqrt13 = sqrt(13)

    weights, bary = s2(
        [(11 - sqrt13) / 360, 0],
        [(10 - 2 * sqrt13) / 45, frac(1, 2)],
        [(29 + 17 * sqrt13) / 360, (7 - sqrt13) / 18],
    )
    return TriangleScheme("Lyness-Jespersen 7", 4, weights, bary)


def lyness_jespersen_8(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt

    sqrt15 = sqrt(15)
    b1, b2 = [(155 - i * sqrt15) / 1200 for i in [+1, -1]]
    r1, r2 = [(6 - i * sqrt15) / 21 for i in [+1, -1]]

    weights, bary = concat(s3(frac(9, 40)), s2([b1, r1], [b2, r2]))
    return TriangleScheme("Lyness-Jespersen 8", 5, weights, bary)


def lyness_jespersen_9(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y

    weights, bary = concat(
        s3(frac(81, 320)),
        s2(
            [frac(1, 90), 0],
            [frac(16, 225), frac(1, 2)],
            [frac(2401, 14400), frac(1, 7)],
        ),
    )
    return TriangleScheme("Lyness-Jespersen 9", 5, weights, bary)


def lyness_jespersen_10():
    weights, bary = concat(
        s2(
            [3.503588271790222e-01 / 3, 2.492867451709329e-01],
            [1.525347191106164e-01 / 3, 6.308901449150177e-02],
        ),
        s1([4.971064537103375e-01 / 6, 6.365024991213939e-01, 5.314504984483216e-02]),
    )
    return TriangleScheme("Lyness-Jespersen 10", 6, weights, bary)


def lyness_jespersen_11(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt

    c, d = [(3 + i * sqrt(6)) / 6 for i in [+1, -1]]

    weights, bary = concat(
        s3(-frac(81, 140)),
        s2(
            [-frac(5, 252), 0],
            [frac(17, 315), frac(1, 2)],
            [frac(128, 315), frac(1, 4)],
        ),
        s1([frac(9, 210), c, d]),
    )
    return TriangleScheme("Lyness-Jespersen 11", 6, weights, bary)


def lyness_jespersen_12():
    weights, bary = concat(
        s3(1.527089667883523e-01),
        s2(
            [2.944076042366762e-01 / 3, 4.738308139536513e-01],
            [3.887052878418766e-01 / 3, 1.721176696308175e-01],
        ),
        s1([1.641781411330949e-01 / 6, 0, 8.653073540834571e-01]),
    )
    return TriangleScheme("Lyness-Jespersen 12", 6, weights, bary)


def lyness_jespersen_13():
    weights, bary = concat(
        s3(-1.495700444677495e-01),
        s2(
            [+5.268457722996328e-01 / 3, 2.603459660790466e-01],
            [+1.600417068265167e-01 / 3, 6.513010290221623e-02],
        ),
        s1([+4.626825653415500e-01 / 6, 6.384441885698096e-01, 4.869031542531756e-02]),
    )
    return TriangleScheme("Lyness-Jespersen 13", 7, weights, bary)


def lyness_jespersen_14():
    weights, bary = concat(
        s3(1.763126156005252e-01),
        s2(
            [1.210901532763310e-02 / 3, 0],
            [3.499561757697094e-01 / 3, 1.549360602237604e-01],
            [3.195119754425220e-01 / 3, 4.691507461438120e-01],
        ),
        s1([1.421102178595603e-01 / 6, 0, 8.392991722729236e-01]),
    )
    return TriangleScheme("Lyness-Jespersen 14", 7, weights, bary)


def lyness_jespersen_15():
    weights, bary = concat(
        s3(1.443156076777862e-01),
        s2(
            [2.852749028018549e-01 / 3, 4.592925882927229e-01],
            [9.737549286959440e-02 / 3, 5.054722831703103e-02],
            [3.096521116041552e-01 / 3, 1.705693077517601e-01],
        ),
        s1([1.633818850466092e-01 / 6, 8.394777409957211e-03, 7.284923929554041e-01]),
    )
    return TriangleScheme("Lyness-Jespersen 15", 8, weights, bary)


def lyness_jespersen_16():
    weights, bary = concat(
        s2(
            [+1.207273935292775e-02 / 3, 0],
            [-8.491579879151455e-01 / 3, 0.5],
            [+1.042367468891334e00 / 3, 4.956813941755582e-01],
            [+1.947229791412260e-01 / 3, 9.032775751426533e-02],
            [+4.511852767201322e-01 / 3, 2.341547497073052e-01],
        ),
        s1([+1.488095238055238e-01 / 6, 0, 7.236067977499750e-01]),
    )
    return TriangleScheme("Lyness-Jespersen 16", 8, weights, bary)


def lyness_jespersen_17():
    weights, bary = concat(
        s3(-2.834183851113958e-01),
        s2(
            [2.097208857979572e-01 / 3, 4.766654393821525e-01],
            [5.127273801480265e-02 / 3, 3.377184405448033e-02],
            [6.564896469913508e-01 / 3, 2.703478891654040e-01],
        ),
        s1([3.659351143072855e-01 / 6, 5.146433548666149e-02, 7.458294907672514e-01]),
    )
    return TriangleScheme("Lyness-Jespersen 17", 8, weights, bary)


def lyness_jespersen_18():
    weights, bary = concat(
        s3(9.713579628279610e-02),
        s2(
            [9.400410068141950e-02 / 3, 4.896825191987370e-01],
            [2.334826230143263e-01 / 3, 4.370895914929355e-01],
            [2.389432167816271e-01 / 3, 1.882035356190322e-01],
            [7.673302697609430e-02 / 3, 4.472951339445297e-02],
        ),
        s1([2.597012362637364e-01 / 6, 3.683841205473626e-02, 7.411985987844980e-01]),
    )
    return TriangleScheme("Lyness-Jespersen 18", 9, weights, bary)


def lyness_jespersen_19():
    weights, bary = concat(
        s3(1.133624844599192e-01),
        s2(
            [1.062573789846330e-03 / 3, 0],
            [4.803411513859279e-02 / 3, 0.5],
            [2.524243006337300e-01 / 3, 4.497793381870162e-01],
            [7.819254371487040e-02 / 3, 4.694744319909033e-02],
            [2.472227459993048e-01 / 3, 1.918719127374489e-01],
        ),
        s1([2.597012362637364e-01 / 6, 3.683841205473626e-02, 7.411985987844980e-01]),
    )
    return TriangleScheme("Lyness-Jespersen 19", 9, weights, bary)


def lyness_jespersen_20():
    weights, bary = concat(
        s2(
            [4.097919300803106e-02 / 3, 3.236494811127173e-02],
            [1.085536215102866e-01 / 3, 1.193509122825931e-01],
            [2.781018986881812e-03 / 3, 5.346110482707572e-01],
            [1.779689321422668e-01 / 3, 2.033099004312816e-01],
            [2.314486047444677e-01 / 3, 3.989693029658558e-01],
        ),
        s1(
            [3.140226717732234e-01 / 6, 5.017813831049474e-02, 5.932012134282132e-01],
            [1.242459578348437e-01 / 6, 2.102201653616613e-02, 8.074890031597923e-01],
        ),
    )
    return TriangleScheme("Lyness-Jespersen 20", 11, weights, bary)


def lyness_jespersen_21():
    weights, bary = concat(
        s3(8.797730116222190e-02),
        s2(
            [2.623293466120857e-02 / 3, 2.598914092828833e-02],
            [1.142447159818060e-01 / 3, 9.428750264792270e-02],
            [5.656634416839376e-02 / 3, 4.946367750172147e-01],
            [2.164790926342230e-01 / 3, 2.073433826145142e-01],
            [2.079874161166116e-01 / 3, 4.389078057004907e-01],
        ),
        s1(
            [4.417430269980344e-02 / 6, 0, 8.588702812826364e-01],
            [2.463378925757316e-01 / 6, 4.484167758913055e-02, 6.779376548825902e-01],
        ),
    )
    return TriangleScheme("Lyness-Jespersen 21", 11, weights, bary)


LynessJespersen = {
    1: lyness_jespersen_1,
    2: lyness_jespersen_2,
    3: lyness_jespersen_3,
    4: lyness_jespersen_4,
    5: lyness_jespersen_5,
    6: lyness_jespersen_6,
    7: lyness_jespersen_7,
    8: lyness_jespersen_8,
    9: lyness_jespersen_9,
    10: lyness_jespersen_10,
    11: lyness_jespersen_11,
    12: lyness_jespersen_12,
    13: lyness_jespersen_13,
    14: lyness_jespersen_14,
    15: lyness_jespersen_15,
    16: lyness_jespersen_16,
    17: lyness_jespersen_17,
    18: lyness_jespersen_18,
    19: lyness_jespersen_19,
    20: lyness_jespersen_20,
    21: lyness_jespersen_21,
}

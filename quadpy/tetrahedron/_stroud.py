
import numpy

from ..helpers import book
from ..line_segment import gauss_legendre
from ..nsimplex._stroud import stroud_tn_5_1
from ._helpers import TetrahedronScheme

citation = book(
    authors=["Arthur Stroud"],
    title="Approximate Calculation of Multiple Integrals",
    publisher="Prentice Hall",
    year="1971",
)


def stroud_t3_5_1():
    return stroud_tn_5_1(3)


def stroud_t3_7_1():
    degree = 7

    gl4 = gauss_legendre(4)
    r = (gl4.points + 1) / 2
    A = gl4.weights / 2

    # Generate Gauss formula for int_0^1 (1-s) * f(s) ds.
    # ```
    # k = numpy.arange(8)
    # moments = 1 / (k**2 + 3*k + 2)
    # alpha, beta = orthopy.line.chebyshev(moments)
    # s, B = orthopy.line.schemes.custom(alpha, beta, mode='numpy')
    # ```
    s = numpy.array(
        [
            5.710419611452533e-02,
            2.768430136381415e-01,
            5.835904323689318e-01,
            8.602401356562251e-01,
        ]
    )
    B = numpy.array(
        [
            1.355069134315012e-01,
            2.034645680102685e-01,
            1.298475476082247e-01,
            3.118097095000554e-02,
        ]
    )

    # Generate Gauss formula for int_0^1 (1-t)^2 * f(t) ds.
    # ```
    # k = numpy.arange(8)
    # moments = 2 / (k**3 + 6*k**2 + 11*k + 6)
    # alpha, beta = orthopy.line.chebyshev(moments)
    # t, C = orthopy.line.schemes.custom(alpha, beta, mode='numpy')
    # ```
    t = numpy.array(
        [
            4.850054944699245e-02,
            2.386007375518456e-01,
            5.170472951043522e-01,
            7.958514178967657e-01,
        ]
    )
    C = numpy.array(
        [
            1.108884156112685e-01,
            1.434587897992167e-01,
            6.863388717292915e-02,
            1.035224074991912e-02,
        ]
    )

    weights = numpy.array(
        [6 * A[i] * B[j] * C[k] for i in range(4) for j in range(4) for k in range(4)]
    )

    points = numpy.array(
        [
            [
                t[k],
                s[j] * (1 - t[k]),
                r[i] * (1 - s[j]) * (1 - t[k]),
                (1 - r[i]) * (1 - s[j]) * (1 - t[k]),
            ]
            for i in range(4)
            for j in range(4)
            for k in range(4)
        ]
    )
    return TetrahedronScheme("Stroud T3 7-1", weights, points, degree, citation)

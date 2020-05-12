import numpy

from ..helpers import book, untangle
from ..u3 import _stroud as sphere_stroud
from ._helpers import E3r2Scheme
from ._stroud_secrest import stroud_secrest_07 as stroud_e3r2_5_1
from ._stroud_secrest import stroud_secrest_08a as stroud_e3r2_5_2a
from ._stroud_secrest import stroud_secrest_08b as stroud_e3r2_5_2b
from ._stroud_secrest import stroud_secrest_09 as stroud_e3r2_5_3
from ._stroud_secrest import stroud_secrest_10a as stroud_e3r2_7_1a
from ._stroud_secrest import stroud_secrest_10b as stroud_e3r2_7_1b
from ._stroud_secrest import stroud_secrest_11a as stroud_e3r2_7_2a
from ._stroud_secrest import stroud_secrest_11b as stroud_e3r2_7_2b

citation = book(
    authors=["Arthur Stroud"],
    title="Approximate Calculation of Multiple Integrals",
    publisher="Prentice Hall",
    year="1971",
)


def stroud_e3r2_14_1(symbolic=False):
    # Get the moments corresponding to monomials and the weight function omega(x) = x^2
    # * exp(-x^2):
    #
    #    int_{-infty}^{infty} x^2 exp(-x^2) x^k dx \
    #
    #           / 0 for k odd,
    #        = {
    #           \ Gamma((k+3)/2) if k even
    #
    # In this particular case, we don't need to compute the recurrence coefficients
    # numerically, but they are given analytically.
    # ```
    # n = 8
    # alpha = numpy.zeros(n)
    # beta = numpy.empty(n)
    # beta[0] = numpy.sqrt(numpy.pi)/2
    # beta[1::2] = numpy.arange(n//2) + 1.5
    # beta[2::2] = numpy.arange(n//2-1) + 1.0
    # points, weights = \
    #     orthopy.line.schemes.custom(alpha, beta, mode='numpy')
    # r = points[-4:]
    # A = weights[-4:]
    # ```
    r = numpy.array(
        [
            7.235510187528402e-01,
            1.468553289216669e00,
            2.266580584531844e00,
            3.190993201781527e00,
        ]
    )
    A = numpy.array(
        [
            2.265043732793035e-01,
            1.908084800858996e-01,
            2.539731378612040e-02,
            4.032955750550135e-04,
        ]
    )

    spherical_scheme = sphere_stroud.stroud_u3_14_1()
    v = spherical_scheme.points
    B = spherical_scheme.weights

    # Normalize the weights to 1
    B /= numpy.sqrt(numpy.pi) / 4

    data = [
        (A[i] * B[j], r[i] * numpy.array([v[j]])) for i in range(4) for j in range(72)
    ]

    points, weights = untangle(data)
    weights *= numpy.sqrt(numpy.pi) ** 3
    return E3r2Scheme("Stroud E3r2 14-1", weights, points, 14, citation)


__all__ = [
    "stroud_e3r2_5_1",
    "stroud_e3r2_5_2a",
    "stroud_e3r2_5_2b",
    "stroud_e3r2_5_3",
    "stroud_e3r2_7_1a",
    "stroud_e3r2_7_1b",
    "stroud_e3r2_7_2a",
    "stroud_e3r2_7_2b",
    "stroud_e3r2_14_1",
]

import numpy
from sympy import Rational as frac
from sympy import pi, sqrt

from ..helpers import article, get_nsimplex_points, fsd
from ._helpers import NBallScheme

citation = article(
    authors=["Srebra B. Stoyanova"],
    title="Cubature formulae of the seventh degree of accuracy for the hypersphere",
    journal="Journal of Computational and Applied Mathematics",
    year="1997",
    pages="15-21",
    url="https://doi.org/10.1016/S0377-0427(97)00094-0",
)


# Compute the area of a hypersphere, <https://tauday.com/tau-manifesto>.
# No roots, no gamma functions. Such a nice recursion!
def surface_hypersphere(n):
    assert n > 0
    if n == 1:
        return 2
    elif n == 2:
        return 2 * pi
    return 2 * pi / (n - 2) * surface_hypersphere(n - 2)


def stoyanova(n, delta=1):
    delta = 0.98945

    pts_a = get_nsimplex_points(n)
    # simplex edge midpoints projected onto the sphere
    pts_b = numpy.array(
        [
            sqrt(frac(n, 2 * (n - 1))) * (pts_a[i1] + pts_a[i2])
            for i1 in range(n + 1)
            for i2 in range(i1)
        ]
    )
    # face midpoints projected onto the sphere
    pts_c = numpy.array(
        [
            sqrt(frac(n, 3 * (n - 2))) * (pts_a[i1] + pts_a[i2] + pts_a[i3])
            for i1 in range(n + 1)
            for i2 in range(i1)
            for i3 in range(i2)
        ]
    )
    # (1/4)-points on the edges connecting the simplex points
    pts_b14 = numpy.array(
        [
            sqrt(frac(n, 10 * n - 6)) * (pts_a[k] + 3 * pts_a[l])
            for k in range(n + 1)
            for l in range(k)
        ] + [
            sqrt(frac(n, 10 * n - 6)) * (3 * pts_a[k] + pts_a[l])
            for k in range(n + 1)
            for l in range(k)
        ]
    )

    n1 = 2 * (n + 1)
    n2 = n * (n + 1)
    n4 = 2 * n * (n + 1)

    print(n)
    print()
    print(len(pts_a), n1)
    print(len(pts_b), n2)

    assert len(pts_a) == n1
    assert len(pts_b) == n2

    exit(1)

    assert n >= 5

    s1 = (n - 2) ** 5 + 243
    s2 = (n - 1) ** 5 + 32
    s3 = (3 * n - 1) ** 6 + (n - 3) ** 6 + 4096 * (n - 1)
    s4 = n ** 5 + 1
    s5 = n ** 2 + 4 * n - 8
    s6 = (n + 1) * (n + 2) * (n + 4)
    s7 = n * (n - 1)
    s8 = n * (n - 2)
    s9 = (n - 1) * (n - 2)
    s10 = (n + 1) * (n + 6) * s6
    s11 = (5 * n - 3) ** 3
    s12 = (n - 1) * (n - 3)
    s13 = n ** 2 - 7 * n + 19
    s14 = n ** 2 - n + 1
    s15 = 2 * (n - 3) * (n + 1) ** 2
    s16 = n ** 3 - 9 * n ** 2 + 33 * n - 38
    s17 = 3 * (n - 2)
    q1 = (
        9
        * (n + 1)
        * (
            10 * n ** 5 * (n + 1) ** 2 * s17
            - n ** 3 * s2 * s5
            + 4 * (n - 1) * s4 * s5 * s7
            - 2 * s4 * s6 * s8
        )
        - n ** 2 * s3 * s8
        + 8 * s4 * s8 * s11
        + 108 * n * s2 * s7 ** 2
        - 432 * (n - 1) ** 3 * s4 * s7
    )
    q2 = (
        n * s1 * s7
        - 3 * s4 * s9 * s17
        + 12 * (n - 1) * s4 * s12
        - 3 * n ** 2 * (n - 3) * s2
    )
    b1 = (
        2
        * (n - 1)
        * (3 * (n + 1) * s5 * s7 * q2 - s12 * q1 - 36 * (n - 1) ** 2 * s7 * q2)
        / (s10 * s17 * q2)
    )
    c1 = s9 * q1 / (2 * s10 * q2)
    e1 = 4 * n * s11 / (9 * s10)
    a1 = n / (n + 6) - b1 - c1 - e1
    y1 = n / (n + 2) - e1 / delta ** 4
    y2 = n / (n + 4) - e1 / delta ** 2
    y3 = 3 * n ** 2 / ((n + 2) * (n + 4)) - e1 * (
        41 * n ** 3 - 101 * n ** 2 + 155 * n - 87
    ) / (2 * (5 * n - 3) ** 2 * delta ** 2)
    w0 = s17 * (y3 * (s15 - s14 * s17) + s13 * s14 * y2) / (s13 * s15 * c1)
    u0 = n * (s17 * y3 - s13 * y2) / (s15 * a1)
    p1 = n * (n - 4) * b1 / (4 * s12 * a1)
    p2 = s16 * s17 * b1 / (4 * s12 * s13 * c1)
    d1 = a1 * p1 ** 2 + c1 * p2 ** 2 + b1
    d2 = a1 * p1 * u0 + c1 * p2 * w0
    d3 = a1 * u0 ** 2 + c1 * w0 ** 2 - y1
    d0 = d2 ** 2 - d1 * d3
    v = (d2 + numpy.sqrt(d0)) / d1  # TODO pm
    u = u0 - p1 * v
    w = w0 - p2 * v
    lmbda2 = 1 / u
    beta2 = 1 / v
    gamma2 = 1 / w
    delta2 = delta ** 2

    lmbda6 = lmbda2 ** 3
    beta6 = beta2 ** 3
    gamma6 = gamma2 ** 3
    delta6 = delta2 ** 3

    n1 = 2 * (n + 1)
    n2 = n * (n + 1)
    n4 = 2 * n * (n + 1)
    if n == 5:
        n3 = cn1 ** 3
    else:
        assert n > 5
        n3 = 2 * cn1 ** 3

    a = a1 / n1 / lmbda6
    b = b1 / n2 / beta6
    c = c1 / n3 / gamma6
    e = e1 / n4 / delta6
    d = 1 - n1 * a - n2 * b - n3 * c - n4 * e

    print(n)
    print()
    print(numpy.sqrt(lmbda2))
    print(numpy.sqrt(beta2))
    print(numpy.sqrt(gamma2))
    print()
    exit(1)
    print(a)

    exit(1)

    return NBallScheme("Stoyanova", n, weights, points, 7, citation)

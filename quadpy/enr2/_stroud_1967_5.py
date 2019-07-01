# TODO sympyfy

import numpy

from ..helpers import article, rd, untangle
from ._helpers import Enr2Scheme

citation = article(
    authors=["A.H. Stroud"],
    title="Some fifth degree integration formulas for symmetric regions II",
    journal="Numerische Mathematik",
    volume="9",
    number="5",
    month="apr",
    year="1967",
    pages="460-468",
    url="https://doi.org/10.1007/BF02162160",
)


def stroud_1967_5_a(n):
    if n == 2:
        eta = 0.446103183094540
        lmbda = 0.136602540378444e1
        xi = -0.366025403784439
        mu = 0.198167882945871e1
        A = 0.328774019778636
        B = 0.833333333333333e-1
        C = 0.455931355469736e-2
    elif n == 3:
        eta = 0.476731294622796
        lmbda = 0.935429018879534
        xi = -0.731237647787132
        mu = 0.433155309477649
        gamma = 0.266922328697744e1
        A = 0.242000000000000
        B = 0.810000000000000e-1
        C = 0.500000000000000e-2
    elif n == 4:
        eta = 0.523945658287507
        lmbda = 0.119433782552719e1
        xi = -0.398112608509063
        mu = -0.318569372920112
        gamma = 0.185675837424096e1
        A = 0.155502116982037
        B = 0.777510584910183e-1
        C = 0.558227484231506e-2
    elif n == 5:
        eta = 0.214972564378798e1
        lmbda = 0.464252986016289e1
        xi = -0.623201054093728
        mu = -0.447108700673434
        gamma = 0.812171426076331
        A = 0.487749259189752e-3
        B = A
        C = 0.497073504444862e-1
    else:
        assert n == 6
        eta = 1.0
        lmbda = numpy.sqrt(2)
        xi = 0.0
        mu = -1.0
        gamma = 1.0
        A = 0.78125e-2
        B = 0.625e-1
        C = A
    # TODO double precision
    # else:
    #     assert n == 7
    #     eta = 0
    #     lmbda = 2.0095056
    #     xi = 0.27745483
    #     mu = -1.06221560
    #     gamma = 0.66983521
    #     A = 1 / 9
    #     B = 1 / 72
    #     C = B

    data = [
        (B, rd(n, [(+lmbda, 1), (+xi, n - 1)])),
        (B, rd(n, [(-lmbda, 1), (-xi, n - 1)])),
    ]

    if n == 2:
        data += [(C, numpy.array([[+mu, +mu]])), (C, numpy.array([[-mu, -mu]]))]
    else:
        data += [
            (C, rd(n, [(+mu, 2), (+gamma, n - 2)])),
            (C, rd(n, [(-mu, 2), (-gamma, n - 2)])),
        ]

    if n == 7:
        data += [(2 * A, numpy.full((1, n), 0.0))]
    else:
        data += [(A, numpy.full((1, n), +eta)), (A, numpy.full((1, n), -eta))]

    points, weights = untangle(data)
    weights *= numpy.sqrt(numpy.pi) ** n

    return Enr2Scheme("Stroud 1967-5 a", n, weights, points, 5, citation)


def stroud_1967_5_b(n):
    if n == 3:
        eta = 0.476731294622796
        lmbda = 0.128679320334269e1
        xi = -0.379873463323979
        mu = -0.192386729447751e1
        # ERR Stroud's book falsely lists 0.31_33_0068...
        gamma = 0.312200683022281
        A = 0.242
        B = 0.81e-1
        C = 0.5e-2
    elif n == 5:
        eta = 0.615369528365158
        lmbda = 0.132894698387445e1
        xi = -0.178394363877324
        mu = -0.745963266507289
        gamma = 0.135503972310817e1
        A = 0.726415024414905e-1
        B = A
        C = 0.641509853510569e-2
    elif n == 6:
        eta = 1.0
        lmbda = 0.942809041582063
        xi = -0.471404520791032
        mu = -0.166666666666667e1
        gamma = 1 / 3
        A = 0.78125e-2
        B = 0.62500e-1
        C = A
    else:
        assert n == 7, "n must be in [3, 5, 6, 7]"
        eta = 0
        lmbda = 0.959724318748357
        xi = -0.772326488820521
        mu = -0.141214270131942e1
        gamma = 0.319908106249452
        A = 1 / 9
        B = 1 / 72
        C = B

    data = [
        (B, rd(n, [(+lmbda, 1), (+xi, n - 1)])),
        (B, rd(n, [(-lmbda, 1), (-xi, n - 1)])),
    ]

    data += [
        (C, rd(n, [(+mu, 2), (+gamma, n - 2)])),
        (C, rd(n, [(-mu, 2), (-gamma, n - 2)])),
    ]

    if n == 7:
        data += [(2 * A, numpy.full((1, n), 0.0))]
    else:
        data += [(A, numpy.full((1, n), +eta)), (A, numpy.full((1, n), -eta))]

    points, weights = untangle(data)
    weights *= numpy.sqrt(numpy.pi) ** n

    return Enr2Scheme("Stroud 1967-5 b", n, weights, points, 5, citation)

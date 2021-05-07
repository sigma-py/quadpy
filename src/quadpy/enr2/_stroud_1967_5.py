# TODO sympyfy

import numpy as np

from ..helpers import article, rd, untangle
from ._helpers import Enr2Scheme

source = article(
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
        B = 1 / 12
        C = 0.455931355469736e-2
    elif n == 3:
        eta = 0.476731294622796
        lmbda = 0.935429018879534
        xi = -0.731237647787132
        mu = 0.433155309477649
        gamma = 0.266922328697744e1
        A = 121 / 500
        B = 81 / 1000
        C = 1 / 200
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
        eta = 1
        lmbda = np.sqrt(2)
        xi = 0
        mu = -1
        gamma = 1
        A = 1 / 128
        B = 1 / 16
        C = A
    # else:
    #     assert n == 7
    #     # ENH These double-precision values have been obtained via a minimization.
    #     # The original reference only has single-prevision.
    #     # TODO make them even more precise
    #     eta = 0
    #     lmbda = 2.009505637083749e+00
    #     xi = 2.774548295173737e-01
    #     mu = -1.062215595206724e+00
    #     gamma = 6.698352123613097e-01
    #     A = 1 / 9
    #     B = 1 / 72
    #     C = B

    data = [
        (B, rd(n, [(+lmbda, 1), (+xi, n - 1)])),
        (B, rd(n, [(-lmbda, 1), (-xi, n - 1)])),
    ]

    if n == 2:
        data += [(C, np.array([[+mu, +mu]])), (C, np.array([[-mu, -mu]]))]
    else:
        data += [
            (C, rd(n, [(+mu, 2), (+gamma, n - 2)])),
            (C, rd(n, [(-mu, 2), (-gamma, n - 2)])),
        ]

    if n == 7:
        data += [(2 * A, np.full((1, n), eta))]
    else:
        data += [(A, np.full((1, n), +eta)), (A, np.full((1, n), -eta))]

    points, weights = untangle(data)
    points = np.ascontiguousarray(points.T)
    return Enr2Scheme("Stroud 1967-5 a", n, weights, points, 5, source)


def stroud_1967_5_b(n):
    if n == 3:
        eta = 0.476731294622796
        lmbda = 0.128679320334269e1
        xi = -0.379873463323979
        mu = -0.192386729447751e1
        # ERR Stroud's book incorrectly lists 0.31_33_0068...
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
        eta = 1
        lmbda = 0.942809041582063
        xi = -0.471404520791032
        mu = -5 / 3
        gamma = 1 / 3
        A = 1 / 128
        B = 1 / 16
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
        data += [(2 * A, np.full((1, n), 0.0))]
    else:
        data += [(A, np.full((1, n), +eta)), (A, np.full((1, n), -eta))]

    points, weights = untangle(data)
    points = np.ascontiguousarray(points.T)
    return Enr2Scheme("Stroud 1967-5 b", n, weights, points, 5, source)

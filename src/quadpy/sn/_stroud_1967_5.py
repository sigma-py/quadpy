import numpy as np

from ..helpers import article, combine, untangle
from ._helpers import SnScheme

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
    if n == 4:
        eta = 0.118301270189222e01
        lmbda = 0.889256950212571
        xi = -0.296418983404190
        mu = -0.531672220881311
        gamma = 0.912205316291852e-01
        A = 0.498588678642297e-03
        B = 0.210826276726545e-01
        C = 0.691951501051233e-01
    elif n == 5:
        eta = 1.0 / 3.0
        lmbda = 0.586530598075514
        xi = -0.229965982852212
        mu = -0.689897948556636
        gamma = 0.126598632371090
        A = 0.535714285714286e-01
        B = A
        C = 0.178571428571429e-01
    elif n == 6:
        eta = 0.395842806007865
        lmbda = 0.495664045946913
        xi = -0.247832022973456
        mu = -0.659738010013108
        gamma = 0.131947602002622
        A = 0.159099570971686e-01
        B = 0.409067810742171e-01
        C = A
    else:
        assert n == 7
        eta = 0.0
        lmbda = 0.409227824522100
        xi = -0.329321121353897
        mu = -0.602139671035296
        gamma = 0.136409274840700
        A = 0.493827160493827e-01
        B = 0.169753086419753e-01
        C = B

    data = [
        (B, combine(((+lmbda,), 1), ((+xi,), n - 1))),
        (B, combine(((-lmbda,), 1), ((-xi,), n - 1))),
        (C, combine(((+mu,), 2), ((+gamma,), n - 2))),
        (C, combine(((-mu,), 2), ((-gamma,), n - 2))),
    ]
    if n == 7:
        data += [(A, np.full((1, n), 0.0))]
    else:
        data += [(A, np.full((1, n), +eta)), (A, np.full((1, n), -eta))]

    points, weights = untangle(data)
    points = np.ascontiguousarray(points.T)
    return SnScheme("Stroud 1967-5 a", n, weights, points, 5, source)


def stroud_1967_5_b(n):
    if n == 4:
        eta = 0.316987298107781
        lmbda = 0.660677533878188
        xi = -0.220225844626063
        mu = -0.715618729879241
        gamma = 0.122780763070861
        A = 0.967236335435799e-01
        B = 0.691951501051233e-01
        C = 0.210826276726545e-01
    elif n == 5:
        eta = 1.0 / 3.0
        lmbda = 0.719863931408848
        xi = -0.966326495188785e-01
        mu = -0.289897948556636
        gamma = 0.526598632371090
        A = 0.535714285714286e-01
        B = A
        C = 0.178571428571429e-01
    elif n == 6:
        eta = 0.395842806007865
        lmbda = 0.743496068920369
        xi = 0.0
        mu = -eta
        gamma = eta
        A = 0.159099570971686e-01
        B = 0.409067810742171e-01
        C = A
    else:
        assert n == 7
        eta = 0.0
        lmbda = 0.856856082693894
        xi = 0.118307136817898
        mu = -0.452930251644698
        gamma = 0.285618694231298
        A = 0.493827160493827e-01
        B = 0.169753086419753e-01
        C = B

    data = [
        (B, combine(((+lmbda,), 1), ((+xi,), n - 1))),
        (B, combine(((-lmbda,), 1), ((-xi,), n - 1))),
        (C, combine(((+mu,), 2), ((+gamma,), n - 2))),
        (C, combine(((-mu,), 2), ((-gamma,), n - 2))),
    ]
    if n == 7:
        data += [(A, np.full((1, n), 0.0))]
    else:
        data += [(A, np.full((1, n), +eta)), (A, np.full((1, n), -eta))]

    points, weights = untangle(data)
    points = np.ascontiguousarray(points.T)
    return SnScheme("Stroud 1967-5 b", n, weights, points, 5, source)

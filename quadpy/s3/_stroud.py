import math

import numpy as np
import sympy

from ..helpers import book, untangle
from ..u3 import _stroud as sphere_stroud
from ._ditkin import ditkin_1 as stroud_5_1
from ._ditkin import ditkin_2 as stroud_5_2
from ._ditkin import ditkin_3 as stroud_7_3
from ._hammer_stroud import hammer_stroud_11_3 as stroud_3_1
from ._hammer_stroud import hammer_stroud_15_3a as stroud_7_1a
from ._hammer_stroud import hammer_stroud_15_3b as stroud_7_1b
from ._helpers import S3Scheme, register
from ._mysovskih import mysovskih as stroud_7_2

_source = book(
    authors=["Arthur Stroud"],
    title="Approximate Calculation of Multiple Integrals",
    publisher="Prentice Hall",
    year="1971",
)

pi = sympy.pi
sqrt = np.vectorize(sympy.sqrt)


def stroud_7_4():
    # spherical product gauss
    # ENH Stroud only gives decimals, sophisticated guesswork gives the analytical
    # expressions.
    pm = np.array([+1, -1])

    # 0.9061798459, 0.5384691101
    alpha, beta = sqrt((35 + pm * 2 * sqrt(70)) / 63)
    rho = np.array([-alpha, -beta, beta, alpha])

    # 0.8611363116, 0.3399810436
    alpha, beta = sqrt((15 + pm * 2 * sqrt(30)) / 35)
    u = np.array([-alpha, -beta, beta, alpha])

    # 0.9238795325, 0.3826834324
    alpha, beta = sqrt((2 + pm * sqrt(2)) / 4)
    v = np.array([-alpha, -beta, beta, alpha])

    # 0.1945553342, 0.1387779991
    alpha, beta = (50 + pm * sqrt(70)) / 300
    A = np.array([alpha, beta, beta, alpha])

    # 0.3478548451, 0.6521451549
    alpha, beta = (18 - pm * sqrt(30)) / 36
    B = np.array([alpha, beta, beta, alpha])

    C = np.full(4, pi / 4)

    def outer3(a, b, c):
        """Given 3 1-dimensional vectors a, b, c, the output is of
        shape (len(a), len(b), len(c)) and contains the values

           out[i, j, k] = a[i] * b[j] * c[k]
        """
        return np.multiply.outer(np.multiply.outer(a, b), c)

    r = outer3(rho, sqrt(1 - u ** 2), sqrt(1 - v ** 2))
    s = outer3(rho, sqrt(1 - u ** 2), v)
    t = outer3(rho, u, 4 * [1])

    data = [
        ((A[i] * B[j] * C[k]), np.array([[r[i][j][k], s[i][j][k], t[i][j][k]]]))
        for i in range(4)
        for j in range(4)
        for k in range(4)
    ]
    points, weights = untangle(data)
    d = {"plain": [weights, points[:, 0], points[:, 1], points[:, 2]]}
    weights /= 4 / 3 * math.pi
    return S3Scheme("Stroud S3 7-4", d, 7, _source)


def stroud_14_1():
    # Get the moments corresponding to the Legendre polynomials and the weight function
    # omega(x) = x^2:
    #
    #                                    / 2/3   if k == 0,
    #    int_{-1}^{+1} |x^2| P_k(x) dx ={  8/45  if k == 2,
    #                                    \ 0     otherwise.
    #
    # In this case, the recurrence coefficients can be determined analytically.
    # ```
    # n = 8
    # alpha = np.full(n, fr(0))
    # k = np.arange(n)
    # beta = np.full(n, fr(0))
    # beta[0] = fr(2, 3)
    # # beta[1::2] = fr((k[1::2]+2)**2, ((2*k[1::2]+2)**2 - 1))
    # for k in range(1, n, 2):
    #     beta[k] = fr((k+2)**2, (2*k+2)**2 - 1)
    # # beta[2::2] = fr(k[2::2]**2, ((2*k[2::2]+2)**2 - 1))
    # for k in range(2, n, 2):
    #     beta[k] = fr(k**2, (2*k+2)**2 - 1)
    #
    # # symbolic computation of the points and weights takes 4orever. Keep an eye on
    # # <https://math.stackexchange.com/q/2450401/36678> for a better algorithm to be
    # # implemented in orthopy.
    # flt = np.vectorize(float)
    # alpha = flt(alpha)
    # beta = flt(beta)
    # points, weights = \
    #     orthopy.line.schemes.custom(alpha, beta, mode='numpy')
    #
    # r = points[-4:]
    # A = weights[-4:]
    # ```
    # TODO get symbolic expressions here
    r = np.array(
        [
            3.242534234038097e-01,
            6.133714327005908e-01,
            8.360311073266362e-01,
            9.681602395076261e-01,
        ]
    )
    A = np.array(
        [
            3.284025994586210e-02,
            9.804813271549834e-02,
            1.262636728646019e-01,
            7.618126780737085e-02,
        ]
    )

    spherical_scheme = sphere_stroud.stroud_u3_14_1()
    v = spherical_scheme.points.T
    B = spherical_scheme.weights

    data = [(A[i] * B[j], r[i] * np.array([v[j]])) for i in range(4) for j in range(72)]

    points, weights = untangle(data)
    d = {"plain": [weights, points[:, 0], points[:, 1], points[:, 2]]}
    weights *= 3
    # weights *= 4 * math.pi
    # weights /= 4 / 3 * math.pi
    return S3Scheme("Stroud S3 14-1", d, 14, _source)


register(
    [
        stroud_3_1,
        stroud_5_1,
        stroud_5_2,
        stroud_7_1a,
        stroud_7_1b,
        stroud_7_2,
        stroud_7_3,
        stroud_7_4,
        stroud_14_1,
    ]
)

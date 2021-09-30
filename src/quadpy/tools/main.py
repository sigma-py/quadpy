"""
[1] Gene H. Golub and John H. Welsch,
    Calculation of Gauss Quadrature Rules,
    Mathematics of Computation,
    Vol. 23, No. 106 (Apr., 1969), pp. 221-230+s1-s10,
    <https://doi.org/10.1090/S0025-5718-69-99647-1>,
    <https://www.semanticscholar.org/paper/Calculation-of-Gauss-Quadrature-Rules-%2A-By-Gene-Golub-Welsch/c715119d5464f614fd8ec590b732ccfea53e72c4?p2df>.

[2] W. Gautschi,
    Algorithm 726: ORTHPOL–a package of routines for generating orthogonal polynomials
    and Gauss-type quadrature rules,
    ACM Transactions on Mathematical Software (TOMS),
    Volume 20, Issue 1, March 1994,
    Pages 21-62,
    <https://doi.org/10.1145/174603.174605>,
    <https://www.cs.purdue.edu/archives/2002/wxg/codes/gauss.m>,

[3] W. Gautschi,
    How and how not to check Gaussian quadrature formulae,
    BIT Numerical Mathematics,
    June 1983, Volume 23, Issue 2, pp 209–216,
    <https://doi.org/10.1007/BF02218441>.

[4] D. Boley and G.H. Golub,
    A survey of matrix inverse eigenvalue problems,
    Inverse Problems, 1987, Volume 3, Number 4,
    <https://doi.org/10.1088/0266-5611/3/4/010>.
"""
from __future__ import annotations

import math

import numpy as np
import sympy
from mpmath import mp
from mpmath.matrices.eigen_symmetric import tridiag_eigen
from scipy.linalg import eigh_tridiagonal
from scipy.linalg.lapack import get_lapack_funcs


def coefficients_from_gauss(points, weights):
    """Given the points and weights of a Gaussian quadrature rule, this method
    reconstructs the recurrence coefficients alpha, beta as appearing in the tridiagonal
    Jacobi matrix tri(b, a, b). This is using "Method 2--orthogonal reduction" from
    (section 3.2 in [4]).  The complexity is O(n^3); a faster method is suggested in 3.3
    in [4].
    """
    n = len(points)
    assert n == len(weights)

    flt = np.vectorize(float)
    points = flt(points)
    weights = flt(weights)

    A = np.zeros((n + 1, n + 1))

    # In sytrd, the _last_ row/column of Q are e, so put the values there.
    a00 = 1.0
    A[n, n] = a00
    k = np.arange(n)
    A[k, k] = points
    A[n, :-1] = np.sqrt(weights)
    A[:-1, n] = np.sqrt(weights)

    # Implemented in
    # <https://github.com/scipy/scipy/issues/7775>
    sytrd, sytrd_lwork = get_lapack_funcs(("sytrd", "sytrd_lwork"))

    # query lwork (optional)
    lwork, info = sytrd_lwork(n + 1)
    assert info == 0

    _, d, e, _, info = sytrd(A, lwork=lwork)
    assert info == 0

    alpha = d[:-1][::-1]
    beta = e[::-1] ** 2
    int_1 = beta[0]
    beta[0] = math.nan

    return alpha, beta, int_1


def _sympy_tridiag(a, b):
    """Creates the tridiagonal sympy matrix tridiag(b, a, b)."""
    n = len(a)
    assert n == len(b)
    A = [[0 for _ in range(n)] for _ in range(n)]
    for i in range(n):
        A[i][i] = a[i]
    for i in range(n - 1):
        A[i][i + 1] = b[i + 1]
        A[i + 1][i] = b[i + 1]
    return sympy.Matrix(A)


def scheme_from_rc(alpha, beta, int_1, mode: str | None = None):
    alpha = np.asarray(alpha)
    beta = np.asarray(beta)

    if mode is None:
        # try and guess the mode
        if alpha.dtype in [np.float32, np.float64]:
            mode = "numpy"
        else:
            raise ValueError(
                'Please specify the `mode` ("sympy", "numpy", or "mpmath").'
            )

    fun = {
        "sympy": _scheme_from_rc_sympy,
        "numpy": _scheme_from_rc_numpy,
        "mpmath": _scheme_from_rc_mpmath,
    }[mode]

    return fun(alpha, beta, int_1)


# Compute the Gauss nodes and weights from the recurrence coefficients associated with a
# set of orthogonal polynomials. See [2] and
# <http://www.scientificpython.net/pyblog/radau-quadrature>.
def _scheme_from_rc_sympy(alpha, beta, int_1):
    beta[0] = int_1

    # Construct the triadiagonal matrix [sqrt(beta), alpha, sqrt(beta)]
    A = _sympy_tridiag(alpha, [sympy.sqrt(bta) for bta in beta])

    # Extract points and weights from eigenproblem
    x = []
    w = []
    for item in A.eigenvects():
        val, multiplicity, vec = item
        assert multiplicity == 1
        assert len(vec) == 1
        vec = vec[0]
        x.append(val)
        norm2 = sum(v ** 2 for v in vec)
        # simplifiction takes really long
        # w.append(sympy.simplify(beta[0] * vec[0]**2 / norm2))
        w.append(beta[0] * vec[0] ** 2 / norm2)
    # sort by x
    order = sorted(range(len(x)), key=lambda i: x[i])
    x = [x[i] for i in order]
    w = [w[i] for i in order]
    return x, w


def _scheme_from_rc_mpmath(alpha, beta, int_1):
    # Create vector
    # cut off the first value of beta (None)
    n = len(alpha)
    b = mp.zeros(n, 1)
    for i in range(n - 1):
        b[i] = mp.sqrt(beta[i + 1])

    z = mp.zeros(1, n)
    z[0, 0] = 1
    d = mp.matrix(alpha)
    tridiag_eigen(mp, d, b, z)

    # nx1 matrix -> list of mpf
    x = np.array([mp.mpf(sympy.N(xx, mp.dps)) for xx in d])
    w = np.array([mp.mpf(sympy.N(int_1, mp.dps)) * mp.power(ww, 2) for ww in z])
    return x, w


def _scheme_from_rc_numpy(alpha, beta, int_1):
    alpha = alpha.astype(np.float64)
    beta = beta.astype(np.float64)
    x, V = eigh_tridiagonal(alpha, np.sqrt(beta[1:]))
    w = int_1 * V[0, :] ** 2
    return x, w

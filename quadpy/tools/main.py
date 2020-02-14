"""
[1] Gene H. Golub and John H. Welsch,
    Calculation of Gauss Quadrature Rules,
    Mathematics of Computation,
    Vol. 23, No. 106 (Apr., 1969), pp. 221-230+s1-s10,
    <https://dx.doi.org/10.2307/2004418>,
    <https://pdfs.semanticscholar.org/c715/119d5464f614fd8ec590b732ccfea53e72c4.pdf>.

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

import numpy
import sympy
from mpmath import mp
from mpmath.matrices.eigen_symmetric import tridiag_eigen
from scipy.linalg import eigh_tridiagonal
from scipy.linalg.lapack import get_lapack_funcs


def golub_welsch(moments):
    """Given moments

    mu_k = int_a^b omega(x) x^k dx,  k = {0, 1,...,2N}

    (with omega being a nonnegative weight function), this method creates the recurrence
    coefficients of the corresponding orthogonal polynomials, see section 4
    ("Determining the Three Term Relationship from the Moments") in Golub-Welsch [1].
    Numerically unstable, see [2].
    """
    assert len(moments) % 2 == 1
    n = (len(moments) - 1) // 2

    M = numpy.array([[moments[i + j] for j in range(n + 1)] for i in range(n + 1)])
    R = numpy.linalg.cholesky(M).T

    # (upper) diagonal
    Rd = R.diagonal()
    q = R.diagonal(1) / Rd[:-1]

    alpha = q.copy()
    alpha[+1:] -= q[:-1]

    # TODO don't square here, but adapt _gauss to accept square-rooted values
    #      as input
    beta = numpy.hstack([Rd[0], Rd[1:-1] / Rd[:-2]]) ** 2
    return alpha, beta


def stieltjes(w, a, b, n):
    t = sympy.Symbol("t")

    alpha = n * [None]
    beta = n * [None]
    mu = n * [None]
    pi = n * [None]

    k = 0
    pi[k] = 1
    mu[k] = sympy.integrate(pi[k] ** 2 * w(t), (t, a, b))
    alpha[k] = sympy.integrate(t * pi[k] ** 2 * w(t), (t, a, b)) / mu[k]
    beta[k] = mu[0]  # not used, by convention mu[0]

    k = 1
    pi[k] = (t - alpha[k - 1]) * pi[k - 1]
    mu[k] = sympy.integrate(pi[k] ** 2 * w(t), (t, a, b))
    alpha[k] = sympy.integrate(t * pi[k] ** 2 * w(t), (t, a, b)) / mu[k]
    beta[k] = mu[k] / mu[k - 1]

    for k in range(2, n):
        pi[k] = (t - alpha[k - 1]) * pi[k - 1] - beta[k - 1] * pi[k - 2]
        mu[k] = sympy.integrate(pi[k] ** 2 * w(t), (t, a, b))
        alpha[k] = sympy.integrate(t * pi[k] ** 2 * w(t), (t, a, b)) / mu[k]
        beta[k] = mu[k] / mu[k - 1]

    return alpha, beta


def chebyshev(moments):
    """Given the first 2n moments `int t^k dt`, this method uses the Chebyshev algorithm
    (see, e.g., [2]) for computing the associated recurrence coefficients.

    WARNING: Ill-conditioned, see [2].
    """
    m = len(moments)
    assert m % 2 == 0
    if isinstance(moments[0], tuple(sympy.core.all_classes)):
        dtype = sympy.Rational
    else:
        dtype = moments.dtype
    zeros = numpy.zeros(m, dtype=dtype)
    return chebyshev_modified(moments, zeros, zeros)


def chebyshev_modified(nu, a, b):
    """Given the first 2n modified moments `nu_k = int p_k(t) dt`, where the p_k are
    orthogonal polynomials with recurrence coefficients a, b, this method implements the
    modified Chebyshev algorithm (see, e.g., [2]) for computing the associated
    recurrence coefficients.
    """
    m = len(nu)
    assert m % 2 == 0, "Need an even number of moments."

    n = m // 2

    alpha = numpy.empty(n, dtype=a.dtype)
    beta = numpy.empty(n, dtype=a.dtype)
    # Actually overkill. One could alternatively make sigma a list, and store the
    # shrinking rows there, only ever keeping the last two.
    sigma = numpy.empty((n, 2 * n), dtype=a.dtype)

    if n > 0:
        k = 0
        sigma[k, k : 2 * n - k] = nu
        alpha[0] = a[0] + nu[1] / nu[0]
        beta[0] = nu[0]

    if n > 1:
        k = 1
        L = numpy.arange(k, 2 * n - k)
        sigma[k, L] = (
            sigma[k - 1, L + 1]
            - (alpha[k - 1] - a[L]) * sigma[k - 1, L]
            + b[L] * sigma[k - 1, L - 1]
        )
        alpha[k] = (
            a[k] + sigma[k, k + 1] / sigma[k, k] - sigma[k - 1, k] / sigma[k - 1, k - 1]
        )
        beta[k] = sigma[k, k] / sigma[k - 1, k - 1]

    for k in range(2, n):
        L = numpy.arange(k, 2 * n - k)
        sigma[k, L] = (
            sigma[k - 1, L + 1]
            - (alpha[k - 1] - a[L]) * sigma[k - 1, L]
            - beta[k - 1] * sigma[k - 2, L]
            + b[L] * sigma[k - 1, L - 1]
        )
        alpha[k] = (
            a[k] + sigma[k, k + 1] / sigma[k, k] - sigma[k - 1, k] / sigma[k - 1, k - 1]
        )
        beta[k] = sigma[k, k] / sigma[k - 1, k - 1]

    return alpha, beta


def integrate(f, a, b):
    """Symbolically calculate the integrals

      int_a^b f_k(x) dx.

    Useful for computing the moments `w(x) * P_k(x)`, e.g.,

    moments = quadpy.tools.integrate(
            lambda x: [x**k for k in range(5)],
            -1, +1
            )
    """
    x = sympy.Symbol("x")
    return numpy.array([sympy.integrate(fun, (x, a, b)) for fun in f(x)])


def coefficients_from_gauss(points, weights):
    """Given the points and weights of a Gaussian quadrature rule, this method
    reconstructs the recurrence coefficients alpha, beta as appearing in the tridiagonal
    Jacobi matrix tri(b, a, b). This is using "Method 2--orthogonal reduction" from
    (section 3.2 in [4]).  The complexity is O(n^3); a faster method is suggested in 3.3
    in [4].
    """
    n = len(points)
    assert n == len(weights)

    flt = numpy.vectorize(float)
    points = flt(points)
    weights = flt(weights)

    A = numpy.zeros((n + 1, n + 1))

    # In sytrd, the _last_ row/column of Q are e, so put the values there.
    a00 = 1.0
    A[n, n] = a00
    k = numpy.arange(n)
    A[k, k] = points
    A[n, :-1] = numpy.sqrt(weights)
    A[:-1, n] = numpy.sqrt(weights)

    # Implemented in
    # <https://github.com/scipy/scipy/issues/7775>
    sytrd, sytrd_lwork = get_lapack_funcs(("sytrd", "sytrd_lwork"))

    # query lwork (optional)
    lwork, info = sytrd_lwork(n + 1)
    assert info == 0

    _, d, e, _, info = sytrd(A, lwork=lwork)
    assert info == 0

    return d[:-1][::-1], e[::-1] ** 2


def check_coefficients(moments, alpha, beta):
    """In his article [3], Walter Gautschi suggests a method for checking if a Gauss
    quadrature rule is sane. This method implements test #3 for the article.
    """
    n = len(alpha)
    assert len(beta) == n

    D = numpy.empty(n + 1)
    Dp = numpy.empty(n + 1)
    D[0] = 1.0
    D[1] = moments[0]
    Dp[0] = 0.0
    Dp[1] = moments[1]
    for k in range(2, n + 1):
        A = numpy.array([moments[i : i + k] for i in range(k)])
        D[k] = numpy.linalg.det(A)
        #
        A[:, -1] = moments[k : 2 * k]
        Dp[k] = numpy.linalg.det(A)

    errors_alpha = numpy.zeros(n)
    errors_beta = numpy.zeros(n)

    errors_alpha[0] = abs(alpha[0] - (Dp[1] / D[1] - Dp[0] / D[0]))
    errors_beta[0] = abs(beta[0] - D[1])
    for k in range(1, n):
        errors_alpha[k] = abs(alpha[k] - (Dp[k + 1] / D[k + 1] - Dp[k] / D[k]))
        errors_beta[k] = abs(beta[k] - D[k + 1] * D[k - 1] / D[k] ** 2)

    return errors_alpha, errors_beta


def _sympy_tridiag(a, b):
    """Creates the tridiagonal sympy matrix tridiag(b, a, b).
    """
    n = len(a)
    assert n == len(b)
    A = [[0 for _ in range(n)] for _ in range(n)]
    for i in range(n):
        A[i][i] = a[i]
    for i in range(n - 1):
        A[i][i + 1] = b[i + 1]
        A[i + 1][i] = b[i + 1]
    return sympy.Matrix(A)


def scheme_from_rc(alpha, beta, mode=None):
    alpha = numpy.asarray(alpha)
    beta = numpy.asarray(beta)

    if mode is None:
        # try and guess the mode
        if alpha.dtype in [numpy.float32, numpy.float64]:
            mode = "numpy"
        else:
            raise ValueError("Please specify the `mode` (\"sympy\", \"numpy\", or \"mpmath\").")

    if mode == "sympy":
        return _scheme_from_rc_sympy(alpha, beta)
    elif mode == "numpy":
        return _scheme_from_rc_numpy(alpha, beta)

    assert mode == "mpmath"
    return _scheme_from_rc_mpmath(alpha, beta)


# Compute the Gauss nodes and weights from the recurrence coefficients associated with a
# set of orthogonal polynomials. See [2] and
# <http://www.scientificpython.net/pyblog/radau-quadrature>.
def _scheme_from_rc_sympy(alpha, beta):
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
        norm2 = sum([v ** 2 for v in vec])
        # simplifiction takes really long
        # w.append(sympy.simplify(beta[0] * vec[0]**2 / norm2))
        w.append(beta[0] * vec[0] ** 2 / norm2)
    # sort by x
    order = sorted(range(len(x)), key=lambda i: x[i])
    x = [x[i] for i in order]
    w = [w[i] for i in order]
    return x, w


def _scheme_from_rc_mpmath(alpha, beta):
    # Create vector cut of the first value of beta
    n = len(alpha)
    b = mp.zeros(n, 1)
    for i in range(n - 1):
        b[i] = mp.sqrt(beta[i + 1])

    z = mp.zeros(1, n)
    z[0, 0] = 1
    d = mp.matrix(alpha)
    tridiag_eigen(mp, d, b, z)

    # nx1 matrix -> list of mpf
    x = numpy.array([mp.mpf(sympy.N(xx, mp.dps)) for xx in d])
    w = numpy.array([mp.mpf(sympy.N(beta[0], mp.dps)) * mp.power(ww, 2) for ww in z])
    return x, w


def _scheme_from_rc_numpy(alpha, beta):
    alpha = alpha.astype(numpy.float64)
    beta = beta.astype(numpy.float64)
    x, V = eigh_tridiagonal(alpha, numpy.sqrt(beta[1:]))
    w = beta[0] * V[0, :] ** 2
    return x, w

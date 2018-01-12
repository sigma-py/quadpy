# -*- coding: utf-8 -*-
#
from mpmath import mp
from mpmath.matrices.eigen_symmetric import tridiag_eigen
import numpy
import scipy
from scipy.linalg import eig_banded
import sympy

from . import recurrence_coefficients
from .. import e1r
from .. import e1r2


def custom(alpha, beta, mode='mpmath', decimal_places=32):
    '''Compute the Gauss nodes and weights from the recurrence coefficients
    associated with a set of orthogonal polynomials. See [2] and
    <http://www.scientificpython.net/pyblog/radau-quadrature>.
    '''

    if mode == 'sympy':
        x, w = _gauss_from_coefficients_sympy(alpha, beta)
    elif mode == 'mpmath':
        x, w = _gauss_from_coefficients_mpmath(alpha, beta, decimal_places)
    else:
        assert mode == 'numpy'
        x, w = _gauss_from_coefficients_numpy(alpha, beta)
    return x, w


def _sympy_tridiag(a, b):
    '''Creates the tridiagonal sympy matrix tridiag(b, a, b).
    '''
    n = len(a)
    assert n == len(b)
    A = [[0 for _ in range(n)] for _ in range(n)]
    for i in range(n):
        A[i][i] = a[i]
    for i in range(n-1):
        A[i][i+1] = b[i+1]
        A[i+1][i] = b[i+1]
    return sympy.Matrix(A)


def _gauss_from_coefficients_sympy(alpha, beta):
    assert isinstance(alpha[0], sympy.Rational)
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
        norm2 = sum([v**2 for v in vec])
        # simplifiction takes really long
        # w.append(sympy.simplify(beta[0] * vec[0]**2 / norm2))
        w.append(beta[0] * vec[0]**2 / norm2)
    # sort by x
    order = sorted(range(len(x)), key=lambda i: x[i])
    x = [x[i] for i in order]
    w = [w[i] for i in order]
    return x, w


def _gauss_from_coefficients_mpmath(alpha, beta, decimal_places):
    mp.dps = decimal_places

    # Create vector cut of the first value of beta
    n = len(alpha)
    b = mp.zeros(n, 1)
    for i in range(n-1):
        # work around <https://github.com/fredrik-johansson/mpmath/issues/382>
        x = beta[i+1]
        if isinstance(x, numpy.int64):
            x = int(x)
        b[i] = mp.sqrt(x)

    z = mp.zeros(1, n)
    z[0, 0] = 1
    d = mp.matrix(alpha)
    tridiag_eigen(mp, d, b, z)

    # nx1 matrix -> list of sympy floats
    x = numpy.array([sympy.Float(xx) for xx in d])
    w = numpy.array([beta[0] * mp.power(ww, 2) for ww in z])
    return x, w


def _gauss_from_coefficients_numpy(alpha, beta):
    assert isinstance(alpha, numpy.ndarray)
    assert isinstance(beta, numpy.ndarray)

    # eigh_tridiagonal is only available from scipy 1.0.0
    try:
        from scipy.linalg import eigh_tridiagonal
    except ImportError:
        # Use eig_banded
        x, V = eig_banded(numpy.vstack((numpy.sqrt(beta), alpha)), lower=False)
        w = beta[0]*scipy.real(scipy.power(V[0, :], 2))
    else:
        x, V = eigh_tridiagonal(alpha, numpy.sqrt(beta[1:]))
        w = beta[0] * V[0, :]**2

    return x, w


def legendre(n, decimal_places):
    _, _, alpha, beta = \
        recurrence_coefficients.legendre(n, 'monic', symbolic=True)
    return custom(alpha, beta, mode='mpmath', decimal_places=decimal_places)


def jacobi(n, a, b, decimal_places):
    _, _, alpha, beta = \
        recurrence_coefficients.jacobi(n, a, b, 'monic', symbolic=True)
    return custom(alpha, beta, mode='mpmath', decimal_places=decimal_places)


def chebyshev1(n, decimal_places):
    # There are explicit representations, too, but for the sake of consistency
    # go for the recurrence coefficients approach here.
    _, _, alpha, beta = \
        recurrence_coefficients.chebyshev1(n, 'monic', symbolic=True)
    beta[0] = sympy.N(beta[0], decimal_places)
    return custom(alpha, beta, mode='mpmath', decimal_places=decimal_places)


def chebyshev2(n, decimal_places):
    # There are explicit representations, too, but for the sake of consistency
    # go for the recurrence coefficients approach here.
    _, _, alpha, beta = \
        recurrence_coefficients.chebyshev2(n, 'monic', symbolic=True)
    beta[0] = sympy.N(beta[0], decimal_places)
    return custom(alpha, beta, mode='mpmath', decimal_places=decimal_places)


def laguerre(n, decimal_places):
    _, _, alpha, beta = \
        e1r.recurrence_coefficients(n, 0, 'monic', symbolic=True)
    return custom(alpha, beta, mode='mpmath', decimal_places=decimal_places)


def laguerre_generalized(n, a, decimal_places):
    _, _, alpha, beta = \
        e1r.recurrence_coefficients(n, a, 'monic', symbolic=True)
    return custom(alpha, beta, mode='mpmath', decimal_places=decimal_places)


def hermite(n, decimal_places):
    _, _, alpha, beta = \
        e1r2.recurrence_coefficients(n, 'monic', symbolic=True)

    # For some reason, the parameters have to be adapted here.
    beta[1:] /= 2

    return custom(alpha, beta, mode='mpmath', decimal_places=decimal_places)

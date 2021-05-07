import numpy as np
import sympy
from mpmath import mp

from ._helpers import C1Scheme


def chebyshev_gauss_1(n, mode="numpy"):
    """Chebyshev-Gauss quadrature for \\int_{-1}^1 f(x) / sqrt(1+x^2) dx."""
    degree = n if n % 2 == 1 else n + 1

    if mode == "numpy":
        points = np.cos((2 * np.arange(1, n + 1) - 1) / (2 * n) * np.pi)
        weights = np.full(n, np.pi / n)
    elif mode == "sympy":
        points = np.array(
            [
                sympy.cos(sympy.Rational(2 * k - 1, 2 * n) * sympy.pi)
                for k in range(1, n + 1)
            ]
        )
        weights = np.full(n, sympy.pi / n)
    else:
        assert mode == "mpmath"
        points = np.array(
            [mp.cos(mp.mpf(2 * k - 1) / (2 * n) * mp.pi) for k in range(1, n + 1)]
        )
        weights = np.full(n, mp.pi / n)
    return C1Scheme("Chebyshev-Gauss 1", degree, weights, points)


def chebyshev_gauss_2(n, mode="numpy", decimal_places=None):
    """Chebyshev-Gauss quadrature for \\int_{-1}^1 f(x) * sqrt(1+x^2) dx."""
    degree = n if n % 2 == 1 else n + 1

    # TODO make explicit for all modes
    if mode == "numpy":
        points = np.cos(np.pi * np.arange(1, n + 1) / (n + 1))
        weights = np.pi / (n + 1) * (np.sin(np.pi * np.arange(1, n + 1) / (n + 1))) ** 2
    elif mode == "sympy":
        points = np.array(
            [sympy.cos(sympy.Rational(k, n + 1) * sympy.pi) for k in range(1, n + 1)]
        )
        weights = np.array(
            [
                sympy.pi / (n + 1) * sympy.sin(sympy.pi * sympy.Rational(k, n + 1)) ** 2
                for k in range(1, n + 1)
            ]
        )
    else:
        assert mode == "mpmath"
        points = np.array(
            [mp.cos(mp.mpf(k) / (n + 1) * mp.pi) for k in range(1, n + 1)]
        )
        weights = np.array(
            [
                mp.pi / (n + 1) * mp.sin(mp.pi * mp.mpf(k) / (n + 1)) ** 2
                for k in range(1, n + 1)
            ]
        )
    return C1Scheme("Chebyshev-Gauss 2", degree, weights, points)

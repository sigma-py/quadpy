# -*- coding: utf-8 -*-
#
import math

from mpmath import mp
import numpy
import pytest
import quadpy

from helpers import check_degree


@pytest.mark.parametrize(
    "scheme,tol", [(quadpy.e1r.gauss_laguerre(n), 1.0e-14) for n in range(1, 10)]
)
def test_scheme(scheme, tol):
    assert scheme.points.dtype == numpy.float64, scheme.name
    assert scheme.weights.dtype == numpy.float64, scheme.name

    degree = check_degree(
        lambda poly: scheme.integrate(poly),
        lambda k: math.factorial(k[0]),
        1,
        scheme.degree + 1,
        tol=tol,
    )
    assert degree == scheme.degree, "Observed: {}   expected: {}".format(
        degree, scheme.degree
    )
    return


@pytest.mark.parametrize("scheme", [quadpy.e1r.gauss_laguerre(1)])
def test_show(scheme):
    scheme.show()
    return


def test_laguerre_mpmath():
    scheme = quadpy.e1r.gauss_laguerre(2, mode="mpmath", decimal_places=51)

    tol = 1.0e-50

    x1 = 2 - mp.sqrt(2)
    x2 = 2 + mp.sqrt(2)
    assert (abs(scheme.points - [x1, x2]) < tol).all()

    w1 = (2 + mp.sqrt(2)) / 4
    w2 = (2 - mp.sqrt(2)) / 4
    assert (abs(scheme.weights - [w1, w2]) < tol).all()
    return


def test_laguerre_generalized_mpmath():
    scheme = quadpy.e1r.gauss_laguerre(2, alpha=1, mode="mpmath", decimal_places=51)

    tol = 1.0e-50

    x1 = 3 - mp.sqrt(3)
    x2 = 3 + mp.sqrt(3)
    assert (abs(scheme.points - [x1, x2]) < tol).all()

    w1 = 2 / ((-1 + mp.sqrt(3)) ** 2 * (1 + 2 / (-1 + mp.sqrt(3)) ** 2))
    w2 = 2 / ((-1 - mp.sqrt(3)) ** 2 * (1 + 2 / (-1 - mp.sqrt(3)) ** 2))
    assert (abs(scheme.weights - [w1, w2]) < tol).all()
    return


if __name__ == "__main__":
    scheme_ = quadpy.e1r.gauss_laguerre(3)
    test_scheme(scheme_, 1.0e-14)
    test_show(scheme_)

# -*- coding: utf-8 -*-
#
from mpmath import mp
import numpy
import pytest
import quadpy

from helpers import check_degree, integrate_monomial_over_enr2


@pytest.mark.parametrize(
    "scheme,tol", [(quadpy.e1r2.gauss_hermite(n), 1.0e-14) for n in range(1, 10)]
)
def test_scheme(scheme, tol):
    assert scheme.points.dtype == numpy.float64, scheme.name
    assert scheme.weights.dtype == numpy.float64, scheme.name

    degree = check_degree(
        lambda poly: scheme.integrate(poly),
        integrate_monomial_over_enr2,
        1,
        scheme.degree + 1,
        tol=tol,
    )
    assert degree == scheme.degree, "Observed: {}   expected: {}".format(
        degree, scheme.degree
    )
    return


@pytest.mark.parametrize("scheme", [quadpy.e1r2.gauss_hermite(2)])
def test_show(scheme):
    scheme.show()
    return


def test_hermite_mpmath():
    scheme = quadpy.e1r2.gauss_hermite(4, mode="mpmath", decimal_places=51)

    tol = 1.0e-50

    x1 = mp.sqrt((3 - mp.sqrt(6)) / 2)
    x2 = mp.sqrt((3 + mp.sqrt(6)) / 2)
    assert (abs(scheme.points - [-x2, -x1, +x1, +x2]) < tol).all()

    w1 = mp.sqrt(mp.pi) / 4 / (3 - mp.sqrt(6))
    w2 = mp.sqrt(mp.pi) / 4 / (3 + mp.sqrt(6))

    assert (abs(scheme.weights - [w2, w1, w1, w2]) < tol).all()
    return


if __name__ == "__main__":
    scheme_ = quadpy.e1r2.gauss_hermite(10)
    test_scheme(scheme_, 1.0e-14)
    test_show(scheme_)

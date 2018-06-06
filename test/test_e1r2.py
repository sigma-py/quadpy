# -*- coding: utf-8 -*-
#
from mpmath import mp
import numpy
import pytest
import quadpy

from helpers import check_degree, integrate_monomial_over_enr2


@pytest.mark.parametrize(
    "scheme,tol", [(quadpy.e1r2.GaussHermite(n), 1.0e-14) for n in range(1, 10)]
)
def test_scheme(scheme, tol):
    assert scheme.points.dtype == numpy.float64, scheme.name
    assert scheme.weights.dtype == numpy.float64, scheme.name

    degree = check_degree(
        lambda poly: quadpy.e1r2.integrate(poly, scheme),
        integrate_monomial_over_enr2,
        1,
        scheme.degree + 1,
        tol=tol,
    )
    assert degree == scheme.degree, "Observed: {}   expected: {}".format(
        degree, scheme.degree
    )
    return


@pytest.mark.parametrize("scheme", [quadpy.e1r2.GaussHermite(2)])
def test_show(scheme):
    quadpy.e1r2.show(scheme)
    return


def test_hermite_mpmath():
    scheme = quadpy.e1r2.GaussHermite(4, mode="mpmath", decimal_places=51)

    tol = 1.0e-50

    x1 = mp.sqrt((3 - mp.sqrt(6)) / 2)
    x2 = mp.sqrt((3 + mp.sqrt(6)) / 2)
    assert (abs(scheme.points - [-x2, -x1, +x1, +x2]) < tol).all()

    w1 = mp.sqrt(mp.pi) / 4 / (3 - mp.sqrt(6))
    w2 = mp.sqrt(mp.pi) / 4 / (3 + mp.sqrt(6))

    assert (abs(scheme.weights - [w2, w1, w1, w2]) < tol).all()
    return


if __name__ == "__main__":
    scheme_ = quadpy.e1r2.GaussHermite(10)
    test_scheme(scheme_, 1.0e-14)
    test_show(scheme_)

# -*- coding: utf-8 -*-
#
import matplotlib.pyplot as plt
import numpy
import pytest
import quadpy

from helpers import check_degree, integrate_monomial_over_enr2


@pytest.mark.parametrize(
    "scheme",
    [
        quadpy.e3r2.stroud_e3r2_5_1(),
        quadpy.e3r2.stroud_e3r2_5_2a(),
        quadpy.e3r2.stroud_e3r2_5_2b(),
        quadpy.e3r2.stroud_e3r2_5_3(),
        quadpy.e3r2.stroud_e3r2_7_1a(),
        quadpy.e3r2.stroud_e3r2_7_1b(),
        quadpy.e3r2.stroud_e3r2_7_2a(),
        quadpy.e3r2.stroud_e3r2_7_2b(),
        quadpy.e3r2.stroud_e3r2_14_1(),
    ],
)
def test_scheme(scheme, tol=1.0e-14):
    assert scheme.points.dtype == numpy.float64, scheme.name
    assert scheme.weights.dtype == numpy.float64, scheme.name

    degree = check_degree(
        lambda poly: scheme.integrate(poly),
        integrate_monomial_over_enr2,
        3,
        scheme.degree + 1,
        tol=tol,
    )
    assert degree == scheme.degree, "Observed: {}   expected: {}".format(
        degree, scheme.degree
    )
    return


@pytest.mark.parametrize("scheme", [quadpy.e3r2.stroud_e3r2_5_1()])
def test_show(scheme, backend="mpl"):
    scheme.show(backend=backend)
    plt.close()
    return


if __name__ == "__main__":
    scheme_ = quadpy.e3r2.Stroud("7-2b")
    test_scheme(scheme_, 1.0e-14)
    test_show(scheme_, backend="vtk")

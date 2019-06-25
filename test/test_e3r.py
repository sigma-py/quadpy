# -*- coding: utf-8 -*-
#
import matplotlib.pyplot as plt
import numpy
import pytest
import accupy

import quadpy

from helpers import check_degree, integrate_monomial_over_enr


@pytest.mark.parametrize(
    "scheme",
    [
        quadpy.e3r.stroud_e3r_5_1(),
        quadpy.e3r.stroud_e3r_5_2(),
        quadpy.e3r.stroud_e3r_5_3(),
        # quadpy.e3r.stroud_e3r_7_1(),
        quadpy.e3r.stroud_e3r_7_2(),
    ],
)
def test_scheme(scheme, tol=1.0e-14):
    assert scheme.points.dtype == numpy.float64, scheme.name
    assert scheme.weights.dtype == numpy.float64, scheme.name

    degree = check_degree(
        lambda poly: scheme.integrate(poly, dot=accupy.fdot),
        integrate_monomial_over_enr,
        3,
        scheme.degree + 1,
        tol=tol,
    )
    assert degree == scheme.degree, "Observed: {}   expected: {}".format(
        degree, scheme.degree
    )
    return


@pytest.mark.parametrize("scheme", [quadpy.e3r.stroud_secrest_x()])
def test_show(scheme, backend="mpl"):
    scheme.show(backend=backend)
    plt.close()
    return


if __name__ == "__main__":
    scheme_ = quadpy.e3r.Stroud("5-1")
    test_scheme(scheme_, 1.0e-14)
    test_show(scheme_, backend="vtk")

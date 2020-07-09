import accupy
import numpy
import pytest
from helpers import check_degree
from matplotlib import pyplot as plt

import quadpy
from quadpy.enr._helpers import integrate_monomial_over_enr


@pytest.mark.parametrize(
    "scheme",
    [
        quadpy.e3r.stroud_e3r_5_1(),
        quadpy.e3r.stroud_e3r_5_2(),
        quadpy.e3r.stroud_e3r_5_3(),
        quadpy.e3r.stroud_e3r_7_1(),
        quadpy.e3r.stroud_e3r_7_2(),
    ],
)
def test_scheme(scheme):
    assert scheme.points.dtype == numpy.float64, scheme.name
    assert scheme.weights.dtype == numpy.float64, scheme.name

    print(scheme)

    degree, err = check_degree(
        lambda poly: scheme.integrate(poly, dot=accupy.fdot),
        integrate_monomial_over_enr,
        3,
        scheme.degree + 1,
        tol=scheme.test_tolerance,
    )
    assert (
        degree >= scheme.degree
    ), "{} -- observed: {}, expected: {} (max err: {:.3e})".format(
        scheme.name, degree, scheme.degree, err
    )


@pytest.mark.parametrize("scheme", [quadpy.e3r.stroud_secrest_10()])
def test_show(scheme, backend="mpl"):
    scheme.show(backend=backend)
    plt.close()


if __name__ == "__main__":
    scheme_ = quadpy.e3r.Stroud("5-1")
    test_scheme(scheme_, 1.0e-14)
    test_show(scheme_, backend="vtk")

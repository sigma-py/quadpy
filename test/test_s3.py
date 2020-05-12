import matplotlib.pyplot as plt
import numpy
import pytest

import quadpy
from helpers import check_degree
from quadpy.sn._helpers import integrate_monomial_over_unit_nball

schemes = [
    quadpy.s3.ditkin_1(),
    quadpy.s3.ditkin_2(),
    quadpy.s3.ditkin_3(),
    quadpy.s3.hammer_stroud_11_3(),
    quadpy.s3.hammer_stroud_12_3(),
    quadpy.s3.hammer_stroud_14_3(),
    quadpy.s3.hammer_stroud_15_3a(),
    quadpy.s3.hammer_stroud_15_3b(),
    quadpy.s3.mysovskih(),
    quadpy.s3.stroud_3_1(),
    quadpy.s3.stroud_5_1(),
    quadpy.s3.stroud_5_2(),
    quadpy.s3.stroud_7_1a(),
    quadpy.s3.stroud_7_1b(),
    quadpy.s3.stroud_7_2(),
    quadpy.s3.stroud_7_3(),
    quadpy.s3.stroud_7_4(),
    quadpy.s3.stroud_14_1(),
]


@pytest.mark.parametrize("scheme", schemes)
def test_scheme(scheme, tol=1.0e-14):
    assert scheme.points.dtype == numpy.float64, scheme.name
    assert scheme.weights.dtype == numpy.float64, scheme.name

    degree = check_degree(
        lambda poly: scheme.integrate(poly, [0.0, 0.0, 0.0], 1.0),
        integrate_monomial_over_unit_nball,
        3,
        scheme.degree + 1,
        tol=tol,
    )
    assert degree == scheme.degree, "Observed: {}   expected: {}".format(
        degree, scheme.degree
    )
    return


@pytest.mark.parametrize("scheme", [quadpy.s3.hammer_stroud_11_3()])
def test_show(scheme, backend="mpl"):
    scheme.show(backend=backend)
    plt.close()
    return


if __name__ == "__main__":
    # scheme_ = quadpy.s3.Stroud("S3 14-1")
    # test_scheme(scheme_, 1.0e-14)
    # test_show(scheme_, backend='vtk')
    from helpers import find_equal

    find_equal(schemes)

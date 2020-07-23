import ndim
import numpy
import pytest
from helpers import check_degree
from matplotlib import pyplot as plt

import quadpy

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
def test_scheme(scheme):
    assert scheme.points.dtype == numpy.float64, scheme.name
    assert scheme.weights.dtype == numpy.float64, scheme.name

    print(scheme)

    degree, err = check_degree(
        lambda poly: scheme.integrate(poly, [0.0, 0.0, 0.0], 1.0),
        ndim.nball.integrate_monomial,
        3,
        scheme.degree + 1,
        tol=scheme.test_tolerance,
    )
    assert (
        degree >= scheme.degree
    ), "{} -- observed: {}, expected: {} (max err: {:.3e})".format(
        scheme.name, degree, scheme.degree, err
    )


@pytest.mark.parametrize("scheme", [quadpy.s3.hammer_stroud_11_3()])
def test_show(scheme, backend="mpl"):
    scheme.show(backend=backend)
    plt.close()


if __name__ == "__main__":
    # scheme_ = quadpy.s3.Stroud("S3 14-1")
    # test_scheme(scheme_, 1.0e-14)
    # test_show(scheme_, backend='vtk')
    from helpers import find_equal

    find_equal(schemes)

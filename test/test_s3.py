import ndim
import numpy
import pytest
from helpers import check_degree
from matplotlib import pyplot as plt

import quadpy


@pytest.mark.parametrize("scheme", quadpy.s3.schemes.values())
def test_scheme(scheme):
    scheme = scheme()

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


@pytest.mark.parametrize("scheme", [quadpy.s3.schemes["hammer_stroud_11_3"]()])
def test_show(scheme, backend="mpl"):
    scheme.show(backend=backend)
    plt.close()


if __name__ == "__main__":
    # scheme_ = quadpy.s3.Stroud("S3 14-1")
    # test_scheme(scheme_, 1.0e-14)
    # test_show(scheme_, backend='vtk')
    from helpers import find_equal

    find_equal(schemes)

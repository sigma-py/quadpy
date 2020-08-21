import numpy
import orthopy
import pytest
from matplotlib import pyplot as plt

import quadpy


@pytest.mark.parametrize("scheme", quadpy.c3.schemes.values())
def test_scheme(scheme, print_degree=False):
    scheme = scheme()

    assert scheme.points.dtype in [numpy.float64, numpy.int64], scheme.name
    assert scheme.weights.dtype in [numpy.float64, numpy.int64], scheme.name

    print(scheme)

    x = [-1.0, +1.0]
    y = [-1.0, +1.0]
    z = [-1.0, +1.0]
    hexa = quadpy.c3.cube_points(x, y, z)

    evaluator = orthopy.cn.Eval(scheme.points)

    k = 0
    while True:
        approximate = scheme.integrate(lambda x: next(evaluator), hexa)
        exact = evaluator.int_p0 * 2 ** 3 if k == 0 else 0.0
        err = numpy.abs(approximate - exact)
        if numpy.any(err > scheme.test_tolerance):
            break
        k += 1

    max_err = numpy.max(err)
    assert k - 1 == scheme.degree, (
        f"{scheme.name} -- observed: {k - 1}, expected: {scheme.degree} "
        f"(max err: {max_err:.3e})"
    )


@pytest.mark.parametrize(
    "scheme", [quadpy.c3.product(quadpy.c1.newton_cotes_closed(2))]
)
def test_show(scheme):
    scheme.show(backend="mpl")
    plt.close()


if __name__ == "__main__":
    # scheme_ = Product(quadpy.c1.NewtonCotesOpen(5))
    # scheme_ = quadpy.c3.HammerStroud("6-3")
    # test_scheme(scheme_, 1.0e-14, print_degree=True)
    # test_show(scheme_)
    # scheme_.show(backend="vtk")
    from helpers import find_equal

    find_equal(schemes)

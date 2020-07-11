import numpy
import orthopy
import pytest
from matplotlib import pyplot as plt

import quadpy

schemes = (
    [quadpy.c3.product(quadpy.c1.midpoint())]
    + [quadpy.c3.product(quadpy.c1.trapezoidal())]
    + [quadpy.c3.product(quadpy.c1.gauss_legendre(k)) for k in range(1, 6)]
    + [quadpy.c3.product(quadpy.c1.newton_cotes_closed(k)) for k in range(1, 5)]
    + [quadpy.c3.product(quadpy.c1.newton_cotes_open(k)) for k in range(5)]
    + [
        quadpy.c3.hammer_stroud_1_3(),
        quadpy.c3.hammer_stroud_2_3(),
        quadpy.c3.hammer_stroud_4_3(),
        quadpy.c3.hammer_stroud_5_3a(),
        quadpy.c3.hammer_stroud_5_3b(),
        quadpy.c3.hammer_stroud_6_3(),
        quadpy.c3.hammer_wymore(),
        quadpy.c3.mustard_lyness_blatt_1(),
        quadpy.c3.mustard_lyness_blatt_2(),
        quadpy.c3.mustard_lyness_blatt_3(),
        quadpy.c3.mustard_lyness_blatt_4(),
        quadpy.c3.mustard_lyness_blatt_5(),
        quadpy.c3.mustard_lyness_blatt_6(),
        quadpy.c3.mustard_lyness_blatt_7(),
        quadpy.c3.sadowsky(),
        quadpy.c3.stroud_c3_3_1(),
        quadpy.c3.stroud_c3_3_2(),
        quadpy.c3.stroud_c3_3_3(),
        quadpy.c3.stroud_c3_3_4(),
        quadpy.c3.stroud_c3_3_5(),
        quadpy.c3.stroud_c3_3_6(),
        quadpy.c3.stroud_c3_3_7(),
        quadpy.c3.stroud_c3_5_1(),
        quadpy.c3.stroud_c3_5_2(),
        quadpy.c3.stroud_c3_5_3(),
        quadpy.c3.stroud_c3_5_4(),
        quadpy.c3.stroud_c3_5_5(),
        quadpy.c3.stroud_c3_5_6(),
        quadpy.c3.stroud_c3_5_7(),
        quadpy.c3.stroud_c3_5_8(),
        quadpy.c3.stroud_c3_7_1a(),
        quadpy.c3.stroud_c3_7_1b(),
        quadpy.c3.stroud_c3_7_2(),
        quadpy.c3.stroud_c3_7_3(),
        quadpy.c3.stroud_1967(),
        quadpy.c3.tyler_1(),
        quadpy.c3.tyler_2(),
    ]
)


@pytest.mark.parametrize("scheme", schemes)
def test_scheme(scheme, print_degree=False):
    assert scheme.points.dtype in [numpy.float64, numpy.int64], scheme.name
    assert scheme.weights.dtype in [numpy.float64, numpy.int64], scheme.name

    print(scheme)

    x = [-1.0, +1.0]
    y = [-1.0, +1.0]
    z = [-1.0, +1.0]
    hexa = quadpy.c3.cube_points(x, y, z)

    evaluator = orthopy.cn.Eval(scheme.points.T)

    k = 0
    while True:
        approximate = scheme.integrate(lambda x: next(evaluator)[0], hexa)
        exact = numpy.sqrt(2.0) ** 3 if k == 0 else 0.0
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

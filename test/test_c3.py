import numpy
import orthopy
import pytest
from matplotlib import pyplot as plt

import quadpy

from helpers import find_best_scheme


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


def test_get_good_scheme():
    degree = 0
    while True:
        best = find_best_scheme(
            quadpy.c3.schemes.values(),
            degree,
            lambda pts: numpy.all((pts >= -1) & (pts <= 1)),
            lambda keys: len(keys - set(["zero", "symm_r00", "symm_rr0", "symm_rrr"])) == 0,
        )
        if best is None:
            break

        # print(best.name)
        b = quadpy.c3.get_good_scheme(degree)
        assert best.name == b.name, f"{best.name} != {b.name}"
        degree += 1

    assert degree == 8


if __name__ == "__main__":
    test_get_good_scheme()

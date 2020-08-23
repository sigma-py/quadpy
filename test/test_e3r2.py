import numpy
import orthopy
import pytest
from matplotlib import pyplot as plt

import quadpy


@pytest.mark.parametrize("scheme", quadpy.e3r2.schemes.values())
def test_scheme(scheme, tol=1.0e-14):
    scheme = scheme()

    assert scheme.points.dtype == numpy.float64, scheme.name
    assert scheme.weights.dtype == numpy.float64, scheme.name

    print(scheme)

    evaluator = orthopy.enr2.Eval(scheme.points, "physicists")

    k = 0
    while True:
        approximate = scheme.integrate(lambda x: next(evaluator))
        exact = evaluator.int_p0 if k == 0 else 0.0
        err = numpy.abs(approximate - exact)
        if numpy.any(err > tol):
            break
        k += 1

    max_err = numpy.max(err)
    assert k - 1 == scheme.degree, (
        f"{scheme.name} -- observed: {k - 1}, expected: {scheme.degree} "
        f"(max err: {max_err:.3e})"
    )


@pytest.mark.parametrize("scheme", [quadpy.e3r2.schemes["stroud_secrest_07"]()])
def test_show(scheme, backend="mpl"):
    scheme.show(backend=backend)
    plt.close()


if __name__ == "__main__":
    scheme_ = quadpy.e3r2.Stroud("7-2b")
    test_scheme(scheme_, 1.0e-14)
    test_show(scheme_, backend="vtk")

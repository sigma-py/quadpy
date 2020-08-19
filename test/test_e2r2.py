import numpy
import orthopy
import pytest

import quadpy


@pytest.mark.parametrize("scheme", quadpy.e2r2.schemes.values())
def test_scheme(scheme):
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
        if numpy.any(err > scheme.test_tolerance * 1.1):
            break
        k += 1

    max_err = numpy.max(err)
    assert k - 1 == scheme.degree, (
        f"{scheme.name} -- observed: {k - 1}, expected: {scheme.degree} "
        f"(max err: {max_err:.3e})"
    )


@pytest.mark.parametrize("scheme", [quadpy.e2r2.schemes["rabinowitz_richter_1"]()])
def test_show(scheme):
    scheme.show()


if __name__ == "__main__":
    pass

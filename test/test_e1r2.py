import numpy
import orthopy
import pytest
from mpmath import mp

import quadpy


@pytest.mark.parametrize(
    "scheme",
    [quadpy.e1r2.gauss_hermite(n) for n in range(5, 12)]
    + [quadpy.e1r2.genz_keister(n) for n in range(8)],
)
def test_scheme(scheme):
    print(scheme.name)
    tol = 1.0e-14
    assert scheme.points.dtype == numpy.float64, scheme.name
    assert scheme.weights.dtype == numpy.float64, scheme.name

    def eval_orthopolys(x):
        return orthopy.e1r2.tree(
            x, scheme.degree + 1, standardization="normal", symbolic=False
        )

    approximate = scheme.integrate(eval_orthopolys)

    exact = numpy.zeros(approximate.shape)
    exact[0] = numpy.sqrt(numpy.sqrt(numpy.pi))

    diff = numpy.abs(approximate - exact)
    k, _ = numpy.where(diff > tol)
    assert len(k) > 0, "{} -- Degree is higher than {}.".format(
        scheme.name, scheme.degree
    )
    degree = k[0] - 1

    assert degree == scheme.degree, "{} -- Observed: {}   expected: {}".format(
        scheme.name, degree, scheme.degree
    )
    return


@pytest.mark.parametrize("scheme", [quadpy.e1r2.gauss_hermite(2)])
def test_show(scheme):
    scheme.show()
    return


def test_hermite_mpmath():
    mp.dps = 51
    scheme = quadpy.e1r2.gauss_hermite(4, mode="mpmath")

    tol = 1.0e-50

    x1 = mp.sqrt((3 - mp.sqrt(6)) / 2)
    x2 = mp.sqrt((3 + mp.sqrt(6)) / 2)
    assert (abs(scheme.points - [-x2, -x1, +x1, +x2]) < tol).all()

    w1 = mp.sqrt(mp.pi) / 4 / (3 - mp.sqrt(6))
    w2 = mp.sqrt(mp.pi) / 4 / (3 + mp.sqrt(6))

    assert (abs(scheme.weights - [w2, w1, w1, w2]) < tol).all()
    return


if __name__ == "__main__":
    scheme_ = quadpy.e1r2.gauss_hermite(10)
    test_scheme(scheme_, 1.0e-14)
    test_show(scheme_)

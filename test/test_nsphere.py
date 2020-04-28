import numpy
import pytest

import quadpy
from helpers import check_degree
from quadpy.nsphere._helpers import integrate_monomial_over_unit_nsphere


@pytest.mark.parametrize(
    "scheme",
    [quadpy.nsphere.dobrodeev_1978(n) for n in range(2, 7)]
    + [quadpy.nsphere.mysovskikh_1(n) for n in range(2, 7)]
    + [quadpy.nsphere.mysovskikh_2(n) for n in range(2, 7)]
    + [quadpy.nsphere.stroud_un_3_1(n) for n in range(2, 7)]
    + [quadpy.nsphere.stroud_un_3_2(n) for n in range(2, 7)]
    + [quadpy.nsphere.stroud_un_5_1(n) for n in range(2, 7)]
    + [quadpy.nsphere.stroud_un_5_2(n) for n in range(2, 7)]
    + [quadpy.nsphere.stroud_un_5_3(n) for n in range(2, 7)]
    + [quadpy.nsphere.stroud_un_5_4(n) for n in range(2, 7)]
    + [quadpy.nsphere.stroud_un_5_4(n) for n in range(2, 7)]
    + [quadpy.nsphere.stroud_un_7_1(n) for n in range(2, 7)]
    + [quadpy.nsphere.stroud_un_7_2(n) for n in range(2, 7)]
    # The scheme has degree 11, so don't push it too far with n. First of all,
    # the number of points increases exponentially, and so does the number of
    # polynomials of degree at most 11.
    + [quadpy.nsphere.stroud_un_11_1(n) for n in range(3, 6)]
    + [quadpy.nsphere.stroud_1967(n) for n in range(2, 7)],
)
def test_scheme(scheme, tol=1.0e-14):
    assert scheme.points.dtype in [numpy.int64, numpy.float64], scheme.name
    assert scheme.weights.dtype in [numpy.int64, numpy.float64], scheme.name

    n = scheme.dim
    center = numpy.zeros(n)
    rad = 1.0
    degree = check_degree(
        lambda poly: scheme.integrate(poly, center, rad),
        integrate_monomial_over_unit_nsphere,
        n,
        scheme.degree + 1,
        tol=tol,
    )
    assert degree >= scheme.degree, "observed: {}, expected: {}".format(
        degree, scheme.degree
    )
    return


if __name__ == "__main__":
    n_ = 5
    scheme_ = quadpy.nsphere.Stroud(n_, "Un 11-1")
    test_scheme(scheme_)

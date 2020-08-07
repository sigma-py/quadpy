import accupy
import ndim
import numpy
import pytest
from helpers import check_degree

import quadpy

schemes = [
    quadpy.e2r.cools_haegemans_9_1(),
    quadpy.e2r.cools_haegemans_9_2(),
    quadpy.e2r.cools_haegemans_13_1(),
    quadpy.e2r.haegemans_piessens_a(),
    quadpy.e2r.haegemans_piessens_b(),
    quadpy.e2r.rabinowitz_richter_1(),
    quadpy.e2r.rabinowitz_richter_2(),
    quadpy.e2r.rabinowitz_richter_3(),
    # quadpy.e2r.rabinowitz_richter_4(),
    quadpy.e2r.rabinowitz_richter_5(),
    quadpy.e2r.stroud_4_1(),
    quadpy.e2r.stroud_5_1(),
    quadpy.e2r.stroud_7_1(),
    quadpy.e2r.stroud_9_1(),
    quadpy.e2r.stroud_11_1(),
    quadpy.e2r.stroud_11_2(),
    quadpy.e2r.stroud_15_1(),
    quadpy.e2r.stroud_secrest_5(),
    quadpy.e2r.stroud_secrest_6(),
]


@pytest.mark.parametrize("scheme", schemes)
def test_scheme(scheme):
    assert scheme.points.dtype == numpy.float64, scheme.name
    assert scheme.weights.dtype == numpy.float64, scheme.name

    print(scheme)

    degree, err = check_degree(
        lambda poly: scheme.integrate(poly, dot=accupy.fdot),
        ndim.enr.integrate_monomial,
        2,
        scheme.degree + 1,
        tol=scheme.test_tolerance * 1.1,
    )
    assert degree >= scheme.degree, (
        f"{scheme.name} -- observed: {degree}, expected: {scheme.degree} "
        f"(max err: {err:.3e})"
    )


@pytest.mark.parametrize("scheme", [quadpy.e2r.rabinowitz_richter_1()])
def test_show(scheme):
    scheme.show()


if __name__ == "__main__":
    # scheme_ = quadpy.e2r.RabinowitzRichter(5)
    # test_scheme(scheme_, 1.0e-14)
    # test_show(scheme_)
    from helpers import find_equal

    find_equal(schemes)

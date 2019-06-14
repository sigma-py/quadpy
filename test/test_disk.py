# -*- coding: utf-8 -*-
#
import numpy
import pytest
import quadpy
from quadpy.nball.helpers import integrate_monomial_over_unit_nball

from helpers import check_degree


@pytest.mark.parametrize(
    "scheme,tol",
    [(scheme(), 1.0e-14) for scheme in quadpy.disk.Albrecht.values()]
    + [(quadpy.disk.CoolsHaegemans[k](), 1.0e-14) for k in range(1, 4)]
    + [(quadpy.disk.CoolsKim[k](), 1.0e-14) for k in range(1, 4)]
    + [(quadpy.disk.HaegemansPiessens(), 1.0e-14)]
    + [(scheme(), 1.0e-14) for scheme in quadpy.disk.HammerStroud.values()]
    + [(quadpy.disk.Lether(k), 1.0e-14) for k in range(1, 6)]
    + [(quadpy.disk.Peirce1957(k), 1.0e-14) for k in range(1, 6)]
    + [(quadpy.disk.PiessensHaegemans(), 1.0e-14)]
    + [(scheme(), 1.0e-14) for scheme in quadpy.disk.RabinowitzRichter.values()]
    + [(scheme(), 1.0e-14) for scheme in quadpy.disk.Stroud.values()]
    + [(quadpy.disk.WissmannBecker[k](), 1.0e-14) for k in ["6-1", "6-2", "8-1"]],
)
def test_scheme(scheme, tol):
    assert scheme.points.dtype == numpy.float64, scheme.name
    assert scheme.weights.dtype == numpy.float64, scheme.name

    degree = check_degree(
        lambda poly: quadpy.disk.integrate(poly, [0.0, 0.0], 1.0, scheme),
        integrate_monomial_over_unit_nball,
        2,
        scheme.degree + 1,
        tol=tol,
    )
    assert degree == scheme.degree, "{}  -- Observed: {}   expected: {}".format(
        scheme.name, degree, scheme.degree
    )
    return


@pytest.mark.parametrize("scheme", [quadpy.disk.Lether(3)])
def test_show(scheme):
    quadpy.disk.show(scheme)
    return


if __name__ == "__main__":
    # scheme_ = quadpy.disk.Lether(5)
    scheme_ = quadpy.disk.Albrecht[8]()
    test_scheme(scheme_, 1.0e-14)
    test_show(scheme_)

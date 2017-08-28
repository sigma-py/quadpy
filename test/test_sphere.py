# -*- coding: utf-8 -*-
#
import numpy
import pytest
import quadpy
from quadpy.nball.helpers import integrate_monomial_over_unit_nsphere
import specialpy

from helpers import check_degree


@pytest.mark.parametrize(
    'scheme',
    [quadpy.sphere.Lebedev(degree) for degree in [
        3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 35, 41, 47, 53,
        59, 65, 71, 77, 83, 89, 95, 101, 107, 113, 119, 125, 131
        ]]
    + [quadpy.sphere.Stroud(k) for k in [
        'U3 3-1',
        'U3 5-1', 'U3 5-2', 'U3 5-3', 'U3 5-4', 'U3 5-5',
        'U3 7-1', 'U3 7-2',
        'U3 8-1',
        'U3 9-1', 'U3 9-2', 'U3 9-3',
        'U3 11-1', 'U3 11-2', 'U3 11-3',
        'U3 14-1',
        ]]
    )
def test_scheme(scheme):
    # Test integration until we get to a polynomial degree `d` that can no
    # longer be integrated exactly. The scheme's degree is `d-1`.
    midpoint = numpy.array([0.0, 0.0, 0.0])
    radius = 1.0
    degree = check_degree(
            lambda poly: quadpy.sphere.integrate(
                poly, midpoint, radius, scheme, sumfun=numpy.sum
                ),
            integrate_monomial_over_unit_nsphere,
            lambda n: quadpy.helpers.partition(n, 3),
            min(30, scheme.degree + 1)
            )
    assert degree >= min(30, scheme.degree), \
        'Observed: {}, expected: {}'.format(degree, scheme.degree)
    return


@pytest.mark.parametrize(
    'scheme',
    [quadpy.sphere.Lebedev(degree) for degree in [
        3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 35, 41, 47, 53,
        # TODO speed up
        # 59, 65, 71, 77, 83, 89, 95, 101, 107, 113, 119, 125, 131
        ]]
    + [quadpy.sphere.Stroud(k) for k in [
        'U3 3-1',
        'U3 5-1', 'U3 5-2', 'U3 5-3', 'U3 5-4', 'U3 5-5',
        'U3 7-1', 'U3 7-2',
        'U3 8-1',
        'U3 9-1', 'U3 9-2', 'U3 9-3',
        'U3 11-1', 'U3 11-2', 'U3 11-3',
        'U3 14-1',
        ]]
    )
# Instead of testing exact integration against of all monomials of degree at
# most l, one can instead test exact integration of all _spherical harmonics_
# of degree at most l. While there are 2**l monomials, there are only l**2
# spherical harmonics.
def test_scheme_spherical(scheme, tol=1.0e-11):
    exact_val = numpy.zeros(scheme.degree + 1)
    exact_val[0] = numpy.sqrt(4*numpy.pi)

    vals = quadpy.sphere.integrate_spherical(
        lambda phi_theta: specialpy.sph_tree(
            scheme.degree+1, phi_theta[1], phi_theta[0]
            ),
        radius=1.0, rule=scheme, sumfun=numpy.sum
        )

    exact = numpy.zeros((scheme.degree+2, scheme.degree+2))
    exact[0, 0] = numpy.sqrt(4 * numpy.pi)

    # The exact value is sqrt(4*pi) for the Y_0^0, and 0 otherwise.
    err = vals
    err[0, 0] -= numpy.sqrt(4.0 * numpy.pi)

    # check in which level the first significant errors occur
    first_error_level = numpy.min(
        numpy.max(numpy.vstack(numpy.where(abs(err) > tol)), axis=0)
        )

    degree = first_error_level - 1

    assert degree == scheme.degree, \
        'Observed: {}, expected: {}'.format(degree, scheme.degree)
    return


@pytest.mark.parametrize(
    'scheme',
    [quadpy.sphere.Lebedev(7)]
    )
def test_show(scheme):
    quadpy.sphere.show(scheme)
    return


if __name__ == '__main__':
    # scheme_ = quadpy.sphere.Stroud('U3 11-3')
    scheme_ = quadpy.sphere.Lebedev(131)
    # test_scheme(scheme_)
    test_scheme_spherical(scheme_)
    # test_show(scheme_)

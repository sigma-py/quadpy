# -*- coding: utf-8 -*-
#
import numpy
import pytest
import quadpy
import specialpy

# Note
# ====
# Instead of testing exact integration against of all monomials of degree at
# most l, one can instead test exact integration of all _spherical harmonics_
# of degree at most l. While there are 2**l monomials, there are only l**2
# spherical harmonics.


def _cartesian_to_spherical(X):
    return numpy.stack([
        numpy.arctan2(X[:, 1], X[:, 0]),
        numpy.arccos(X[:, 2])
        ], axis=1)


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
def test_scheme_cartesian(scheme, tol=1.0e-11):
    exact_val = numpy.zeros(scheme.degree + 1)
    exact_val[0] = numpy.sqrt(4*numpy.pi)

    def sph_tree_cartesian(x):
        print(x.shape)
        phi_theta = _cartesian_to_spherical(x.T).T
        return specialpy.sph_tree(
            scheme.degree+1, phi_theta[1], phi_theta[0]
            )

    vals = quadpy.sphere.integrate(
        sph_tree_cartesian,
        center=numpy.array([0.0, 0.0, 0.0]),
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


# Test a few schemes with integrate_spherical. -- This is basically the same as
# above, no need to repeat it all in detail.
@pytest.mark.parametrize(
    'scheme',
    [quadpy.sphere.Lebedev(degree) for degree in [
        3, 5, 7, 9, 11, 13, 15, 17, 19,
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
    test_scheme_cartesian(scheme_)
    # test_show(scheme_)

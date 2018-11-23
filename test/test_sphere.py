# -*- coding: utf-8 -*-
#
import numpy
import pytest
import orthopy
import quadpy
from quadpy.sphere.helpers import cartesian_to_spherical

# Note
# ====
# Instead of testing exact integration against of all monomials of degree at
# most l, one can instead test exact integration of all _spherical harmonics_
# of degree at most l. While there are 2**l monomials, there are only l**2
# spherical harmonics.


@pytest.mark.parametrize(
    "scheme", [quadpy.sphere.Lebedev("3a"), quadpy.sphere.Stroud("U3 14-1")]
)
def test_spherical_harmonic(scheme):
    """Assert the norm of the spherical harmonic

    Y_1^1(phi, theta) = -1/2 sqrt(3/2/pi) * exp(i*phi) * sin(theta)

    is indeed 1, i.e.,

    int_0^2pi int_0^pi
        Y_1^1(phi, theta) * conj(Y_1^1(phi, theta)) * sin(theta)
        dphi dtheta = 1.
    """

    def spherical_harmonic_11(azimuthal, polar):
        # y00 = 1.0 / numpy.sqrt(4*numpy.pi)
        y11 = (
            -0.5
            * numpy.sqrt(3.0 / 2.0 / numpy.pi)
            * numpy.exp(1j * azimuthal)
            * numpy.sin(polar)
        )
        return y11 * numpy.conjugate(y11)

    val = quadpy.sphere.integrate_spherical(spherical_harmonic_11, rule=scheme)

    assert abs(val - 1.0) < 1.0e-14
    return


@pytest.mark.parametrize(
    "scheme,tol",
    [
        (quadpy.sphere.HeoXu(index), 1.0e-6)
        for index in [
            "13",
            "15",
            "17",
            "19-1",
            "19-2",
            "21-1",
            "21-2",
            "21-3",
            "21-4",
            "21-5",
            "21-6",
            "23-1",
            "23-2",
            "23-3",
            "25-1",
            "25-2",
            "27-1",
            "27-2",
            "27-3",
            "29",
            "31",
            "33",
            "35",
            "37",
            "39-1",
            "39-2",
        ]
    ]
    + [(quadpy.sphere.FliegeMaier(index), 1.0e-6) for index in ["4", "9", "16", "25"]]
    + [
        (quadpy.sphere.Lebedev(index), 1.0e-11)
        for index in [
            "3a",
            "3b",
            "3c",
            "5",
            "7",
            "9",
            "11",
            "13",
            "15",
            "17",
            "19",
            "21",
            "23",
            "25",
            "27",
            "29",
            "31",
            "35",
            "41",
            "47",
            "53",
            "59",
            "65",
            "71",
            "77",
            "83",
            "89",
            "95",
            "101",
            "107",
            "113",
            "119",
            # The highest degree formulas are too memory-intensive for circleci,
            # and the tests are oom-killed. A workaround would be to not test the
            # entire tree at once, but split it up.
            # Check <https://stackoverflow.com/q/47474140/353337>.
            # TODO reenable
            # "125", "131"
        ]
    ]
    + [
        (quadpy.sphere.Stroud(k), 1.0e-13)
        for k in [
            "U3 3-1",
            "U3 5-1",
            "U3 5-2",
            "U3 5-3",
            "U3 5-4",
            "U3 5-5",
            "U3 7-1",
            "U3 7-2",
            "U3 8-1",
            "U3 9-1",
            "U3 9-2",
            "U3 9-3",
            "U3 11-1",
            "U3 11-3",
            "U3 14-1",
        ]
    ]
    + [(quadpy.sphere.Stroud(k), 1.0e-12) for k in ["U3 11-2"]],
)
def test_scheme_cartesian(scheme, tol):
    exact_val = numpy.zeros(scheme.degree + 1)
    exact_val[0] = numpy.sqrt(4 * numpy.pi)

    def sph_tree_cartesian(x):
        flt = numpy.vectorize(float)
        azimuthal, polar = cartesian_to_spherical(flt(x).T).T
        return numpy.concatenate(
            orthopy.sphere.tree_sph(
                polar, azimuthal, scheme.degree + 1, standardization="quantum mechanic"
            )
        )

    assert scheme.points.dtype == numpy.float64, scheme.name
    assert scheme.weights.dtype == numpy.float64, scheme.name

    vals = quadpy.sphere.integrate(
        sph_tree_cartesian, center=numpy.array([0.0, 0.0, 0.0]), radius=1.0, rule=scheme
    )

    # The exact value is sqrt(4*pi) for the Y_0^0, and 0 otherwise.
    err = vals
    err[0] -= numpy.sqrt(4.0 * numpy.pi)

    # check in which level the first significant errors occur
    k = 0
    first_error_level = None
    for L in range(scheme.degree + 2):
        m = 2 * L + 1
        if numpy.any(abs(err[k : k + m]) > tol):
            first_error_level = L
            break
        k += m
    assert first_error_level is not None

    degree = first_error_level - 1

    assert degree == scheme.degree, "{}  --  Observed: {}, expected: {}".format(
        scheme.name, degree, scheme.degree
    )
    return


# Test a few schemes with integrate_spherical. -- This is basically the same as
# above, no need to repeat it all in detail.
@pytest.mark.parametrize(
    "scheme,tol",
    [
        (quadpy.sphere.HeoXu(index), 1.0e-6)
        for index in [
            "13",
            "15",
            "17",
            "19-1",
            "19-2",
            "21-1",
            "21-2",
            "21-3",
            "21-4",
            "21-5",
            "21-6",
            "23-1",
            "23-2",
            "23-3",
            "25-1",
            "25-2",
            "27-1",
            "27-2",
            "27-3",
            "29",
            "31",
            "33",
            "35",
            "37",
            "39-1",
            "39-2",
        ]
    ]
    + [
        (quadpy.sphere.Lebedev(index), 1.0e-11)
        for index in ["3a", "3b", "3c", "5", "7", "9", "11", "13", "15", "17", "19"]
    ]
    + [
        (quadpy.sphere.Stroud(k), 1.0e-13)
        for k in [
            "U3 3-1",
            "U3 5-1",
            "U3 5-2",
            "U3 5-3",
            "U3 5-4",
            "U3 5-5",
            "U3 7-1",
            "U3 7-2",
            "U3 8-1",
            "U3 9-1",
            "U3 9-2",
            "U3 9-3",
            "U3 11-1",
            "U3 11-3",
            "U3 14-1",
        ]
    ]
    + [(quadpy.sphere.Stroud(k), 1.0e-12) for k in ["U3 11-2"]],
)
def test_scheme_spherical(scheme, tol):
    exact_val = numpy.zeros(scheme.degree + 1)
    exact_val[0] = numpy.sqrt(4 * numpy.pi)

    def sph_tree(azimuthal, polar):
        return numpy.concatenate(
            orthopy.sphere.tree_sph(
                polar, azimuthal, scheme.degree + 1, standardization="quantum mechanic"
            )
        )

    vals = quadpy.sphere.integrate_spherical(sph_tree, rule=scheme)

    # The exact value is sqrt(4*pi) for the Y_0^0, and 0 otherwise.
    err = vals
    err[0] -= numpy.sqrt(4 * numpy.pi)

    # check in which level the first significant errors occur
    k = 0
    first_error_level = None
    for L in range(scheme.degree + 2):
        m = 2 * L + 1
        if numpy.any(abs(err[k : k + m]) > tol):
            first_error_level = L
            break
        k += m
    assert first_error_level is not None

    degree = first_error_level - 1

    assert degree == scheme.degree, "Observed: {}, expected: {}".format(
        degree, scheme.degree
    )
    return


@pytest.mark.parametrize("scheme", [quadpy.sphere.Lebedev("7")])
def test_show(scheme):
    quadpy.sphere.show(scheme)
    return


if __name__ == "__main__":
    scheme_ = quadpy.sphere.Stroud("U3 5-2")
    # test_scheme(scheme_)
    test_scheme_spherical(scheme_, tol=1.0e-7)
    # test_show(scheme_)

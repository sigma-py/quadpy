import numpy
import pytest
import sympy
from helpers import check_degree

import quadpy
from quadpy.tn._helpers import integrate_monomial_over_unit_simplex


def _integrate_exact(f, t3):
    #
    # Note that
    #
    #     \int_T f(x) dx = \int_T0 |J(xi)| f(P(xi)) dxi
    #
    # with
    #
    #     P(xi) = x0 * (1-xi[0]-xi[1]) + x1 * xi[0] + x2 * xi[1].
    #
    # and T0 being the reference t3 [(0.0, 0.0), (1.0, 0.0), (0.0,
    # 1.0)].
    # The determinant of the transformation matrix J equals twice the volume of
    # the t3. (See, e.g.,
    # <http://math2.uncc.edu/~shaodeng/TEACHING/math5172/Lectures/Lect_15.PDF>).
    #
    xi = sympy.DeferredVector("xi")
    x_xi = (
        +t3[0] * (1 - xi[0] - xi[1] - xi[2])
        + t3[1] * xi[0]
        + t3[2] * xi[1]
        + t3[3] * xi[2]
    )
    abs_det_J = 6 * quadpy.t3.volume(t3)
    exact = sympy.integrate(
        sympy.integrate(
            sympy.integrate(abs_det_J * f(x_xi), (xi[2], 0, 1 - xi[0] - xi[1])),
            (xi[1], 0, 1 - xi[0]),
        ),
        (xi[0], 0, 1),
    )
    return float(exact)


@pytest.mark.parametrize("scheme", quadpy.t3.schemes.values())
def test_scheme(scheme):
    print(scheme)
    scheme = scheme()

    assert scheme.points.dtype in [numpy.float64, numpy.int64], scheme.name
    assert scheme.weights.dtype in [numpy.float64, numpy.int64], scheme.name

    print(scheme)

    # Test integration until we get to a polynomial degree `d` that can no
    # longer be integrated exactly. The scheme's degree is `d-1`.
    t3 = numpy.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    )

    degree, err = check_degree(
        lambda poly: scheme.integrate(poly, t3),
        integrate_monomial_over_unit_simplex,
        3,
        scheme.degree + 1,
        scheme.test_tolerance * 1.1,
    )

    assert (
        degree >= scheme.degree
    ), "{} -- observed: {}, expected: {} (max err: {:.3e})".format(
        scheme.name, degree, scheme.degree, err
    )


@pytest.mark.skip(reason="gh-actions's python cannot use system vtk")
@pytest.mark.parametrize("scheme", [quadpy.t3.schemes["hammer_marlowe_stroud_3"]()])
def test_show(scheme):
    tet = numpy.array(
        [
            [numpy.cos(0.5 * numpy.pi), numpy.sin(0.5 * numpy.pi), -0.5],
            [numpy.cos(7.0 / 6.0 * numpy.pi), numpy.sin(7.0 / 6.0 * numpy.pi), -0.5],
            [numpy.cos(11.0 / 6.0 * numpy.pi), numpy.sin(11.0 / 6.0 * numpy.pi), -0.5],
            [0.0, 0.0, 1.0],
        ]
    )
    scheme.show(tet, render=False)


if __name__ == "__main__":
    # scheme_ = quadpy.t3.Stroud("T3 7-1")
    # test_scheme(scheme_)
    # # test_show(scheme_)
    # quadpy.t3.show(scheme_, backend="vtk")
    from helpers import find_equal

    find_equal(quadpy.t3.schemes)

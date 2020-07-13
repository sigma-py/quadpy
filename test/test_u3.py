import numpy
import orthopy
import pytest
from matplotlib import pyplot as plt

import quadpy

# Note
# ====
# Instead of testing exact integration against of all monomials of degree at most l, one
# can instead test exact integration of all _spherical harmonics_ of degree at most l.
# While there are 2 ** l monomials, there are only l ** 2 spherical harmonics.


@pytest.mark.parametrize(
    "scheme", [quadpy.u3.lebedev_003a(), quadpy.u3.stroud_u3_14_1()]
)
def test_spherical_harmonic(scheme):
    """Assert the norm of the spherical harmonic

    Y_1^1(phi, theta) = -1/2 sqrt(3/2/pi) * exp(i*phi) * sin(theta)

    is indeed 1, i.e.,

    int_0^2pi int_0^pi
        Y_1^1(phi, theta) * conj(Y_1^1(phi, theta)) * sin(theta)
        dphi dtheta = 1.
    """

    def spherical_harmonic_11(theta_phi):
        # y00 = 1.0 / numpy.sqrt(4*numpy.pi)
        theta, phi = theta_phi
        y11 = (
            -0.5
            * numpy.sqrt(3.0 / 2.0 / numpy.pi)
            * numpy.exp(1j * phi)
            * numpy.sin(theta)
        )
        return y11 * numpy.conjugate(y11)

    val = scheme.integrate_spherical(spherical_harmonic_11)

    assert abs(val - 1.0) < 1.0e-14


@pytest.mark.parametrize(
    "scheme",
    [quadpy.u3.bazant_oh_09(), quadpy.u3.bazant_oh_11(), quadpy.u3.bazant_oh_13()]
    + [
        quadpy.u3.heo_xu_13(),
        quadpy.u3.heo_xu_15(),
        quadpy.u3.heo_xu_17(),
        quadpy.u3.heo_xu_19a(),
        quadpy.u3.heo_xu_19_2(),
        quadpy.u3.heo_xu_21_1(),
        quadpy.u3.heo_xu_21_2(),
        quadpy.u3.heo_xu_21_3(),
        quadpy.u3.heo_xu_21_4(),
        quadpy.u3.heo_xu_21_5(),
        quadpy.u3.heo_xu_21_6(),
        quadpy.u3.heo_xu_23_1(),
        quadpy.u3.heo_xu_23_2(),
        quadpy.u3.heo_xu_23_3(),
        quadpy.u3.heo_xu_25_1(),
        quadpy.u3.heo_xu_25_2(),
        quadpy.u3.heo_xu_27_1(),
        quadpy.u3.heo_xu_27_2(),
        quadpy.u3.heo_xu_27_3(),
        quadpy.u3.heo_xu_29(),
        quadpy.u3.heo_xu_31(),
        quadpy.u3.heo_xu_33(),
        quadpy.u3.heo_xu_35(),
        quadpy.u3.heo_xu_37(),
        quadpy.u3.heo_xu_39_1(),
        quadpy.u3.heo_xu_39_2(),
    ]
    # + [
    #     quadpy.u3.fliege_maier_04(),
    #     quadpy.u3.fliege_maier_09(),
    #     quadpy.u3.fliege_maier_16(),
    #     quadpy.u3.fliege_maier_25(),
    # ]
    # + [
    #     quadpy.u3.lebedev_003a(),
    #     quadpy.u3.lebedev_003b(),
    #     quadpy.u3.lebedev_003c(),
    #     quadpy.u3.lebedev_005(),
    #     quadpy.u3.lebedev_007(),
    #     quadpy.u3.lebedev_009(),
    #     quadpy.u3.lebedev_011(),
    #     quadpy.u3.lebedev_013(),
    #     quadpy.u3.lebedev_015(),
    #     quadpy.u3.lebedev_017(),
    #     quadpy.u3.lebedev_019(),
    #     quadpy.u3.lebedev_021(),
    #     quadpy.u3.lebedev_023(),
    #     quadpy.u3.lebedev_025(),
    #     quadpy.u3.lebedev_027(),
    #     quadpy.u3.lebedev_029(),
    #     quadpy.u3.lebedev_031(),
    #     quadpy.u3.lebedev_035(),
    #     quadpy.u3.lebedev_041(),
    #     quadpy.u3.lebedev_047(),
    #     quadpy.u3.lebedev_053(),
    #     quadpy.u3.lebedev_059(),
    #     quadpy.u3.lebedev_065(),
    #     quadpy.u3.lebedev_071(),
    #     quadpy.u3.lebedev_077(),
    #     quadpy.u3.lebedev_083(),
    #     quadpy.u3.lebedev_089(),
    #     quadpy.u3.lebedev_095(),
    #     quadpy.u3.lebedev_101(),
    #     quadpy.u3.lebedev_107(),
    #     quadpy.u3.lebedev_113(),
    #     quadpy.u3.lebedev_119(),
    #     quadpy.u3.lebedev_125(),
    #     quadpy.u3.lebedev_131(),
    # ]
    # + [
    #     quadpy.u3.stroud_u3_3_1(),
    #     quadpy.u3.stroud_u3_5_1(),
    #     quadpy.u3.stroud_u3_5_2(),
    #     quadpy.u3.stroud_u3_5_3(),
    #     quadpy.u3.stroud_u3_5_4(),
    #     quadpy.u3.stroud_u3_5_5(),
    #     quadpy.u3.stroud_u3_7_1(),
    #     quadpy.u3.stroud_u3_7_2(),
    #     quadpy.u3.stroud_u3_8_1(),
    #     quadpy.u3.stroud_u3_9_1(),
    #     quadpy.u3.stroud_u3_9_2(),
    #     quadpy.u3.stroud_u3_9_3(),
    #     quadpy.u3.stroud_u3_11_1(),
    #     quadpy.u3.stroud_u3_11_2(),
    #     quadpy.u3.stroud_u3_11_3(),
    #     quadpy.u3.stroud_u3_14_1(),
    # ],
)
def test_scheme_cartesian(scheme):
    assert scheme.points.dtype == numpy.float64, scheme.name
    assert scheme.weights.dtype == numpy.float64, scheme.name

    print(scheme)

    # We're using the iterator here; it's much less memory-intensive than computing the
    # full tree at once.
    evaluator = orthopy.u3.EvalCartesian(scheme.points.T, "quantum mechanic")

    k = 0
    while True:
        approximate = scheme.integrate(lambda x: next(evaluator), [0.0, 0.0, 0.0], 1.0)
        exact = numpy.sqrt(4 * numpy.pi) if k == 0 else 0.0
        err = numpy.abs(approximate - exact)
        if numpy.any(err > scheme.test_tolerance):
            break
        k += 1
    max_err = numpy.max(err)

    # find the max error across all polynomials
    max_err = 0.0
    evaluator = orthopy.u3.EvalCartesian(scheme.points.T, "quantum mechanic")
    for i in range(scheme.degree + 1):
        approximate = scheme.integrate(lambda x: next(evaluator), [0.0, 0.0, 0.0], 1.0)
        exact = numpy.sqrt(4 * numpy.pi) if i == 0 else 0.0
        err = numpy.abs(approximate - exact)
        max_err = max(max_err, numpy.max(err))

    assert k - 1 == scheme.degree, (
        f"{scheme.name} -- observed: {k - 1}, expected: {scheme.degree} "
        f"(max err: {max_err:.3e})"
    )


# Test a few schemes with integrate_spherical. -- This is basically the same as above,
# no need to repeat it all in detail.
@pytest.mark.parametrize(
    "scheme",
    [
        quadpy.u3.heo_xu_13(),
        quadpy.u3.heo_xu_15(),
        quadpy.u3.heo_xu_17(),
        quadpy.u3.heo_xu_19a(),
        quadpy.u3.heo_xu_19_2(),
        quadpy.u3.heo_xu_21_1(),
        quadpy.u3.heo_xu_21_2(),
        quadpy.u3.heo_xu_21_3(),
        quadpy.u3.heo_xu_21_4(),
        quadpy.u3.heo_xu_21_5(),
        quadpy.u3.heo_xu_21_6(),
        quadpy.u3.heo_xu_23_1(),
        quadpy.u3.heo_xu_23_2(),
        quadpy.u3.heo_xu_23_3(),
        quadpy.u3.heo_xu_25_1(),
        quadpy.u3.heo_xu_25_2(),
        quadpy.u3.heo_xu_27_1(),
        quadpy.u3.heo_xu_27_2(),
        quadpy.u3.heo_xu_27_3(),
        quadpy.u3.heo_xu_29(),
        quadpy.u3.heo_xu_31(),
        quadpy.u3.heo_xu_33(),
        quadpy.u3.heo_xu_35(),
        quadpy.u3.heo_xu_37(),
        quadpy.u3.heo_xu_39_1(),
        quadpy.u3.heo_xu_39_2(),
    ]
    + [
        quadpy.u3.lebedev_003a(),
        quadpy.u3.lebedev_003b(),
        quadpy.u3.lebedev_003c(),
        quadpy.u3.lebedev_005(),
        quadpy.u3.lebedev_007(),
        quadpy.u3.lebedev_009(),
        quadpy.u3.lebedev_011(),
        quadpy.u3.lebedev_013(),
        quadpy.u3.lebedev_015(),
        quadpy.u3.lebedev_017(),
        quadpy.u3.lebedev_019(),
    ]
    + [
        quadpy.u3.stroud_u3_3_1(),
        quadpy.u3.stroud_u3_5_1(),
        quadpy.u3.stroud_u3_5_2(),
        quadpy.u3.stroud_u3_5_3(),
        quadpy.u3.stroud_u3_5_4(),
        quadpy.u3.stroud_u3_5_5(),
        quadpy.u3.stroud_u3_7_1(),
        quadpy.u3.stroud_u3_7_2(),
        quadpy.u3.stroud_u3_8_1(),
        quadpy.u3.stroud_u3_9_1(),
        quadpy.u3.stroud_u3_9_2(),
        quadpy.u3.stroud_u3_9_3(),
        quadpy.u3.stroud_u3_11_1(),
        quadpy.u3.stroud_u3_11_2(),
        quadpy.u3.stroud_u3_11_3(),
        quadpy.u3.stroud_u3_14_1(),
    ],
)
def test_scheme_spherical(scheme):
    print(scheme)

    evaluator = orthopy.u3.EvalSpherical(scheme.theta_phi.T, "quantum mechanic")

    k = 0
    while True:
        approximate = scheme.integrate_spherical(lambda theta_phi: next(evaluator))
        exact = numpy.sqrt(4 * numpy.pi) if k == 0 else 0.0
        err = numpy.abs(approximate - exact)
        if numpy.any(err > scheme.test_tolerance):
            break
        k += 1
    max_err = numpy.max(err)

    assert k - 1 == scheme.degree, (
        f"{scheme.name} -- observed: {k - 1}, expected: {scheme.degree} "
        f"(max err: {max_err:.3e})"
    )


@pytest.mark.parametrize("scheme", [quadpy.u3.lebedev_007()])
def test_show(scheme):
    scheme.show()
    plt.close()


if __name__ == "__main__":
    scheme_ = quadpy.u3.Stroud("U3 5-2")
    # test_scheme(scheme_)
    test_scheme_spherical(scheme_, tol=1.0e-7)
    # test_show(scheme_)

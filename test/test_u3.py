import numpy
import orthopy
import pytest
from helpers import find_best_scheme
from matplotlib import pyplot as plt

import quadpy

# Note
# ====
# Instead of testing exact integration against of all monomials of degree at most l, one
# can instead test exact integration of all _spherical harmonics_ of degree at most l.
# While there are 2 ** l monomials, there are only l ** 2 spherical harmonics.


@pytest.mark.parametrize(
    "scheme", [quadpy.u3.schemes["lebedev_003a"](), quadpy.u3.schemes["mclaren_10"]()],
)
def test_spherical_harmonic(scheme):
    """Assert the norm of the spherical harmonic

    Y_1^1(phi, theta) = -1/2 sqrt(3/2/pi) * exp(i*phi) * sin(theta)

    is indeed

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


@pytest.mark.parametrize("scheme", quadpy.u3.schemes.values())
def test_scheme_cartesian(scheme):
    # initialize
    scheme = scheme()

    assert scheme.points.dtype == numpy.float64, scheme.name
    assert scheme.weights.dtype == numpy.float64, scheme.name

    print(scheme)

    # assert contiguous x
    def f(x):
        assert x.flags["C_CONTIGUOUS"]
        assert x.shape[0] == 3
        return numpy.ones(x.shape[1:])

    scheme.integrate(f, [0, 0, 0], 1)

    res = scheme.compute_residuals(scheme.degree + 1)
    deg = numpy.where(res > scheme.test_tolerance * 1.1)[0][0] - 1
    if deg != scheme.degree:
        max_res = max(res[:-1])
        raise AssertionError(
            f"{scheme.name} -- observed: {deg}, expected: {scheme.degree} "
            f"(max residual: {max_res})"
        )


@pytest.mark.parametrize("scheme", quadpy.u3.schemes.values())
def test_scheme_spherical(scheme):
    # initialize
    scheme = scheme()

    print(scheme)

    evaluator = orthopy.u3.EvalSpherical(scheme.theta_phi, "quantum mechanic")

    k = 0
    while True:
        approximate = scheme.integrate_spherical(lambda theta_phi: next(evaluator))
        exact = numpy.sqrt(4 * numpy.pi) if k == 0 else 0.0
        err = numpy.abs(approximate - exact)
        if numpy.any(err > scheme.test_tolerance * 1.1):
            break
        k += 1
    max_err = numpy.max(err)

    assert k - 1 == scheme.degree, (
        f"{scheme.name} -- observed: {k - 1}, expected: {scheme.degree} "
        f"(max err: {max_err:.3e})"
    )


@pytest.mark.parametrize("scheme", [quadpy.u3.schemes["lebedev_007"]()])
def test_show(scheme):
    scheme.show()
    plt.close()


def test_get_good_scheme():
    degree = 0
    while True:
        best = find_best_scheme(
            quadpy.u3.schemes.values(),
            degree,
            # lambda pts: True,
            lambda pts: numpy.all(
                numpy.abs(pts[0] ** 2 + pts[1] ** 2 + pts[2] ** 2 - 1.0) < 1.0e-13
            ),
            lambda keys: len(
                keys
                - set(["a1", "a2", "a3", "rs0", "llm2", "rsw2", "pq0", "llm", "rsw"])
            )
            == 0,
        )
        if best is None:
            break

        b = quadpy.u3.get_good_scheme(degree)

        assert best.name == b.name, f"{best.name} != {b.name}"
        degree += 1

    assert degree == 48


if __name__ == "__main__":
    test_get_good_scheme()

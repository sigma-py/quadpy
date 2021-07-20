import ndim
import numpy as np
import pytest
import sympy
from helpers import check_degree

import quadpy


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
    # the t3.
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

    assert scheme.points.dtype in [np.float64, np.int64], scheme.name
    assert scheme.weights.dtype in [np.float64, np.int64], scheme.name

    print(scheme)

    # Test integration until we get to a polynomial degree `d` that can no
    # longer be integrated exactly. The scheme's degree is `d-1`.
    t3 = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])

    degree, err = check_degree(
        lambda poly: scheme.integrate(poly, t3),
        ndim.nsimplex.integrate_monomial,
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
    tet = np.array(
        [
            [np.cos(0.5 * np.pi), np.sin(0.5 * np.pi), -0.5],
            [np.cos(7.0 / 6.0 * np.pi), np.sin(7.0 / 6.0 * np.pi), -0.5],
            [np.cos(11.0 / 6.0 * np.pi), np.sin(11.0 / 6.0 * np.pi), -0.5],
            [0.0, 0.0, 1.0],
        ]
    )
    scheme.show(tet, backend="vtk", render=False)


@pytest.mark.skip()
def test_get_good_scheme():
    for degree in range(51):
        best = None
        for scheme in quadpy.t3.schemes.values():
            scheme = scheme()  # initialize

            # filter schemes for eligibility
            if scheme.degree < degree:
                continue

            # allow only positive weights
            if any(scheme.weights < 0):
                continue

            # disallow points outside of the domain
            if np.any(scheme.points < 0):
                continue

            if scheme.test_tolerance > 1.0e-13:
                continue

            try:
                keys = set(scheme.symmetry_data.keys())
            except AttributeError:
                continue

            if len(keys - {"s4", "s31", "s22", "s211", "s1111"}) > 0:
                continue

            # okay, now compare the scheme with `best`
            if best is None:
                best = scheme
                continue

            if len(scheme.weights) > len(best.weights):
                continue
            elif len(scheme.weights) < len(best.weights):
                best = scheme
                continue
            else:  # len(scheme.weights) == len(best.weights):
                abs_weights = np.abs(scheme.weights)
                ratio = max(abs_weights) / min(abs_weights)
                bratio = max(np.abs(best.weights)) / min(np.abs(best.weights))
                if ratio < bratio:
                    best = scheme
                    continue
                elif ratio > bratio:
                    continue
                else:  # ratio == bratio
                    # # check if it's actually the same scheme
                    # if np.all(np.abs(scheme.points - best.points) < 1.0e-12):
                    #     print("DUP", best.name, scheme.name)
                    #     # pick the older one

                    # for all intents and purposes, the schemes are equal; take the
                    # older one
                    scheme_year = "0" if scheme.source is None else scheme.source.year
                    best_year = "0" if best.source is None else best.source.year
                    if scheme_year < best_year:
                        best = scheme
                        continue
                    elif scheme_year > best_year:
                        continue
                    else:  # years are equal
                        pass

        print(degree, best.name)

        # print(best)
    return


if __name__ == "__main__":
    test_get_good_scheme()

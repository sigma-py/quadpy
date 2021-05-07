import numpy as np
import orthopy
import pytest
from helpers import find_best_scheme

import quadpy


@pytest.mark.parametrize("scheme", quadpy.t2.schemes.values())
def test_scheme(scheme):
    try:
        scheme = scheme()  # initialize
    except TypeError:
        scheme = scheme(5)

    assert scheme.points.dtype in [np.float64, np.int64], scheme.name
    assert scheme.weights.dtype in [np.float64, np.int64], scheme.name

    print(scheme)

    triangle = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])

    evaluator = orthopy.t2.Eval(scheme.points, "normal")

    # assert contiguous x
    def f(x):
        assert x.flags["C_CONTIGUOUS"]
        assert x.shape[0] == 2
        return np.ones(x.shape[1:])

    approximate = scheme.integrate(f, triangle)

    k = 0
    max_err = 0.0
    while True:
        approximate = scheme.integrate(lambda x: next(evaluator), triangle)
        exact = evaluator.int_p0 * 0.5 if k == 0 else 0.0
        err = np.abs(approximate - exact)
        max_err = max(max_err, np.max(err))
        if np.any(err > scheme.test_tolerance * 1.1):
            break
        k += 1

    if k - 1 != scheme.degree:
        # find the max error across all polynomials
        for i in range(k + 1, scheme.degree + 1):
            approximate = scheme.integrate(lambda x: next(evaluator), triangle)
            exact = evaluator.int_p0 * 0.5 if i == 0 else 0.0
            err = np.abs(approximate - exact)
            max_err = max(max_err, np.max(err))

        raise AssertionError(
            f"{scheme.name} -- observed: {k - 1}, expected: {scheme.degree} "
            f"(max err: {max_err:.3e})"
        )


@pytest.mark.parametrize("scheme", [quadpy.t2.get_good_scheme(10)])
def test_show(scheme):
    triangle = np.array(
        [
            [np.cos(0.5 * np.pi), np.sin(0.5 * np.pi)],
            [np.cos(7.0 / 6.0 * np.pi), np.sin(7.0 / 6.0 * np.pi)],
            [np.cos(11.0 / 6.0 * np.pi), np.sin(11.0 / 6.0 * np.pi)],
        ]
    )
    scheme.show(triangle)


def test_volume():
    # Assert computation of triangle volume in 3D is correct
    triangle = np.array([[0.0, 0.0, 0.0], [1.0, 2.0, 3.0], [0.7, 0.4, 1.1]])
    ref = np.sqrt(3.0) / 2.0
    assert abs(quadpy.t2.get_vol(triangle) - ref) < 1.0e-14 * ref

    triangle = np.array([[0.0, 0.0, 0.0], [0.3, 0.4, 0.5], [0.7, 0.4, 1.1]])
    ref = np.sqrt(0.0209)
    assert abs(quadpy.t2.get_vol(triangle) - ref) < 1.0e-14 * ref


def test_multidim():
    scheme = quadpy.t2.schemes["dunavant_05"]()

    np.random.seed(0)
    # simple scalar integration
    tri = np.random.rand(3, 2)
    val = scheme.integrate(lambda x: np.sin(x[0]), tri)
    assert val.shape == ()

    # scalar integration on 4 subdomains
    tri = np.random.rand(3, 4, 2)
    val = scheme.integrate(lambda x: np.sin(x[0]), tri)
    assert val.shape == (4,)

    # scalar integration in 4D
    tri = np.random.rand(3, 4)
    val = scheme.integrate(lambda x: np.sin(x[0]), tri)
    assert val.shape == ()

    # vector-valued integration on 4 subdomains
    tri = np.random.rand(3, 4, 2)
    val = scheme.integrate(lambda x: [np.sin(x[0]), np.cos(x[1])], tri)
    assert val.shape == (2, 4)

    # vector-valued integration in 4D
    tri = np.random.rand(3, 4)
    val = scheme.integrate(lambda x: [np.sin(x[0]), np.cos(x[1])], tri)
    assert val.shape == (2,)

    # # another vector-valued integration in 3D
    # # This is one case where the integration routine may not properly recognize the
    # # dimensionality of the domain. Use the `dim` parameter.
    # val = scheme.integrate(
    #     lambda x: [
    #         x[0] + np.sin(x[1]),
    #         np.cos(x[0]) * x[2],
    #         np.sin(x[0]) + x[1] + x[2],
    #     ],
    #     [[0.0, 1.0, 2.0], [1.0, 2.0, 3.0]],
    #     dim=1,
    # )
    # assert val.shape == (3,)


def test_get_good_scheme():
    degree = 0
    while True:
        best = find_best_scheme(
            quadpy.t2.schemes.values(),
            degree,
            lambda pts: np.all(pts >= 0),
            lambda keys: len(keys - {"d3_aa", "d3_ab", "centroid", "vertex"}) == 0,
        )
        if best is None:
            break

        b = quadpy.t2.get_good_scheme(degree)

        assert best.name == b.name, f"{best.name} != {b.name}"
        degree += 1

    assert degree == 51


if __name__ == "__main__":
    test_get_good_scheme()
    # test_multidim()
    # scheme_ = quadpy.t2.WandzuraXiao(3)
    # test_scheme(scheme_, 1.0e-14)
    # test_show(scheme_)
    # from helpers import find_equal
    # schemes_ = [scheme[0] for scheme in schemes_tol]
    # find_equal(schemes_)

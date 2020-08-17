import numpy
import orthopy
import pytest

import quadpy


@pytest.mark.parametrize("scheme", quadpy.t2.schemes.values())
def test_scheme(scheme):
    try:
        scheme = scheme()  # initialize
    except TypeError:
        scheme = scheme(5)

    assert scheme.points.dtype in [numpy.float64, numpy.int64], scheme.name
    assert scheme.weights.dtype in [numpy.float64, numpy.int64], scheme.name

    print(scheme)

    triangle = numpy.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])

    evaluator = orthopy.t2.Eval(scheme.points, "normal")

    # assert contiguous x
    def f(x):
        assert x.flags["C_CONTIGUOUS"]
        assert x.shape[0] == 2
        return numpy.ones(x.shape[1:])

    approximate = scheme.integrate(f, triangle)

    k = 0
    max_err = 0.0
    while True:
        approximate = scheme.integrate(lambda x: next(evaluator), triangle)
        exact = evaluator.int_p0 * 0.5 if k == 0 else 0.0
        err = numpy.abs(approximate - exact)
        max_err = max(max_err, numpy.max(err))
        if numpy.any(err > scheme.test_tolerance * 1.1):
            break
        k += 1

    if k - 1 != scheme.degree:
        # find the max error across all polynomials
        for i in range(k + 1, scheme.degree + 1):
            approximate = scheme.integrate(lambda x: next(evaluator), triangle)
            exact = evaluator.int_p0 * 0.5 if i == 0 else 0.0
            err = numpy.abs(approximate - exact)
            max_err = max(max_err, numpy.max(err))

        raise AssertionError(
            f"{scheme.name} -- observed: {k - 1}, expected: {scheme.degree} "
            f"(max err: {max_err:.3e})"
        )


@pytest.mark.parametrize("scheme", [quadpy.t2.get_good_scheme(10)])
def test_show(scheme):
    triangle = numpy.array(
        [
            [numpy.cos(0.5 * numpy.pi), numpy.sin(0.5 * numpy.pi)],
            [numpy.cos(7.0 / 6.0 * numpy.pi), numpy.sin(7.0 / 6.0 * numpy.pi)],
            [numpy.cos(11.0 / 6.0 * numpy.pi), numpy.sin(11.0 / 6.0 * numpy.pi)],
        ]
    )
    scheme.show(triangle)


def test_volume():
    # Assert computation of triangle volume in 3D is correct
    triangle = numpy.array([[0.0, 0.0, 0.0], [1.0, 2.0, 3.0], [0.7, 0.4, 1.1]])
    ref = numpy.sqrt(3.0) / 2.0
    assert abs(quadpy.t2.get_vol(triangle) - ref) < 1.0e-14 * ref

    triangle = numpy.array([[0.0, 0.0, 0.0], [0.3, 0.4, 0.5], [0.7, 0.4, 1.1]])
    ref = numpy.sqrt(0.0209)
    assert abs(quadpy.t2.get_vol(triangle) - ref) < 1.0e-14 * ref


def test_multidim():
    scheme = quadpy.t2.schemes["dunavant_05"]()

    numpy.random.seed(0)
    # simple scalar integration
    tri = numpy.random.rand(3, 2)
    val = scheme.integrate(lambda x: numpy.sin(x[0]), tri)
    assert val.shape == ()

    # scalar integration on 4 subdomains
    tri = numpy.random.rand(3, 4, 2)
    val = scheme.integrate(lambda x: numpy.sin(x[0]), tri)
    assert val.shape == (4,)

    # scalar integration in 4D
    tri = numpy.random.rand(3, 4)
    val = scheme.integrate(lambda x: numpy.sin(x[0]), tri)
    assert val.shape == ()

    # vector-valued integration on 4 subdomains
    tri = numpy.random.rand(3, 4, 2)
    val = scheme.integrate(lambda x: [numpy.sin(x[0]), numpy.cos(x[1])], tri)
    assert val.shape == (2, 4)

    # vector-valued integration in 4D
    tri = numpy.random.rand(3, 4)
    val = scheme.integrate(lambda x: [numpy.sin(x[0]), numpy.cos(x[1])], tri)
    assert val.shape == (2,)

    # # another vector-valued integration in 3D
    # # This is one case where the integration routine may not properly recognize the
    # # dimensionality of the domain. Use the `dim` parameter.
    # val = scheme.integrate(
    #     lambda x: [
    #         x[0] + numpy.sin(x[1]),
    #         numpy.cos(x[0]) * x[2],
    #         numpy.sin(x[0]) + x[1] + x[2],
    #     ],
    #     [[0.0, 1.0, 2.0], [1.0, 2.0, 3.0]],
    #     dim=1,
    # )
    # assert val.shape == (3,)


@pytest.mark.skip()
def test_get_good_scheme():
    for degree in range(51):
        best = None
        for scheme in quadpy.t2.schemes.values():
            try:
                scheme = scheme()  # initialize
            except TypeError:
                scheme = scheme(5)

            # filter schemes for eligibility
            if scheme.degree < degree:
                continue

            # allow only positive weights
            if any(scheme.weights < 0):
                continue

            # disallow points outside of the domain
            if numpy.any(scheme.points < 0):
                continue

            keys = set(scheme.symmetry_data.keys())
            if len(keys - set(["d3_aa", "d3_ab", "centroid", "vertex"])) > 0:
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
                abs_weights = numpy.abs(scheme.weights)
                ratio = max(abs_weights) / min(abs_weights)
                bratio = max(numpy.abs(best.weights)) / min(numpy.abs(best.weights))
                if ratio < bratio:
                    best = scheme
                    continue
                elif ratio > bratio:
                    continue
                else:  # ratio == bratio
                    # # check if it's actually the same scheme
                    # if numpy.all(numpy.abs(scheme.points - best.points) < 1.0e-12):
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
                    else:   # years are equal
                        pass

        print(degree, best.name)

        # print(best)
    return


if __name__ == "__main__":
    test_get_good_scheme()
    # test_multidim()
    # scheme_ = quadpy.t2.WandzuraXiao(3)
    # test_scheme(scheme_, 1.0e-14)
    # test_show(scheme_)
    # from helpers import find_equal
    # schemes_ = [scheme[0] for scheme in schemes_tol]
    # find_equal(schemes_)

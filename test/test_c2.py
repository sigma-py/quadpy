import numpy
import orthopy
import pytest

import quadpy

schemes = (
    list(quadpy.c2.schemes.values())
    + [quadpy.c2.product(quadpy.c1.midpoint())]
    + [quadpy.c2.product(quadpy.c1.trapezoidal())]
    + [quadpy.c2.product(quadpy.c1.gauss_legendre(k)) for k in range(1, 5)]
    + [quadpy.c2.product(quadpy.c1.newton_cotes_closed(k)) for k in range(1, 5)]
    + [quadpy.c2.product(quadpy.c1.newton_cotes_open(k)) for k in range(1, 6)]
)


@pytest.mark.parametrize("scheme", schemes)
def test_scheme(scheme):
    # instantiate
    try:
        scheme = scheme()
    except TypeError:
        pass

    assert scheme.points.dtype in [numpy.float64, numpy.int64], scheme.name
    assert scheme.weights.dtype in [numpy.float64, numpy.int64], scheme.name

    print(scheme)

    quad = quadpy.c2.rectangle_points([-1.0, +1.0], [-1.0, +1.0])

    evaluator = orthopy.cn.Eval(scheme.points)

    k = 0
    max_err = 0.0
    while True:
        approximate = scheme.integrate(lambda x: next(evaluator), quad)
        exact = evaluator.int_p0 * 4 if k == 0 else 0.0
        err = numpy.abs(approximate - exact)
        max_err = max(max_err, numpy.max(err))
        if numpy.any(err > scheme.test_tolerance * 1.1):
            break
        k += 1

    if k - 1 != scheme.degree:
        # find the max error across all polynomials
        for i in range(k + 1, scheme.degree + 1):
            approximate = scheme.integrate(lambda x: next(evaluator), quad)
            exact = 2.0 if i == 0 else 0.0
            err = numpy.abs(approximate - exact)
            max_err = max(max_err, numpy.max(err))

        raise AssertionError(
            f"{scheme.name} -- observed: {k - 1}, expected: {scheme.degree} "
            f"(max err: {max_err:.3e})"
        )


@pytest.mark.parametrize("scheme", [quadpy.c2.product(quadpy.c1.gauss_legendre(5))])
def test_show(scheme):
    scheme.show()


# def test_get_good_scheme():
#     for degree in range(22):
#         best = None
#         for scheme in quadpy.c2.schemes.values():
#             try:
#                 scheme = scheme()  # initialize
#             except TypeError:
#                 scheme = scheme(5)
#
#             if scheme.degree < degree:
#                 continue
#
#             # allow only positive weights
#             if any(scheme.weights < 0):
#                 continue
#
#             # disallow points outside of the domain
#             if numpy.any((scheme.points < -1) | (scheme.points > 1)):
#                 continue
#
#             if scheme.test_tolerance > 1.0e-13:
#                 continue
#
#             # TODO force symmetry data for all schemes
#             try:
#                 keys = set(scheme.symmetry_data.keys())
#             except AttributeError:
#                 continue
#
#             # filter out disallowed (unsymmetrical) keys
#             if len(keys - set(["c4_a0", "c4_aa", "d4", "zero"])) > 0:
#                 continue
#
#             if best is not None:
#                 if len(scheme.weights) > len(best.weights):
#                     continue
#                 if len(scheme.weights) == len(best.weights):
#                     ratio = max(numpy.abs(scheme.weights)) / min(
#                         numpy.abs(scheme.weights)
#                     )
#                     bratio = max(numpy.abs(best.weights)) / min(numpy.abs(best.weights))
#                     if ratio > bratio:
#                         continue
#                     # check if it's actually the same scheme
#                     if numpy.all(numpy.abs(scheme.points - best.points) < 1.0e-12):
#                         continue
#
#             # okay, looks like we found a better one!
#             best = scheme
#
#         print(degree, best.name)
#         # print(best)
#     return


if __name__ == "__main__":
    test_get_good_scheme()
    # scheme_ = Product(quadpy.c1.gauss_legendre(6))
    # scheme_ = quadpy.c2.HammerStroud("3-2")
    # scheme_ = quadpy.c2.Stroud["C2 3-2"]()
    # test_show(scheme_)
    # test_scheme(scheme_, 1.0e-14)
    # from helpers import find_equal
    # schemes_ = [scheme[0] for scheme in schemes]
    # find_equal(schemes_)

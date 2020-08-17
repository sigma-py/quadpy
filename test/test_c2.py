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


def test_get_good_scheme():
    for degree in range(22):
        best = None
        for scheme in quadpy.c2.schemes.values():
            # filter schemes for eligibility
            try:
                scheme = scheme()  # initialize
            except TypeError:
                scheme = scheme(5)

            if scheme.degree < degree:
                continue

            # allow only positive weights
            if any(scheme.weights < 0):
                continue

            # disallow points outside of the domain
            if numpy.any((scheme.points < -1) | (scheme.points > 1)):
                continue

            if scheme.test_tolerance > 1.0e-13:
                continue

            # TODO force symmetry data for all schemes
            try:
                keys = set(scheme.symmetry_data.keys())
            except AttributeError:
                continue

            # filter out disallowed (unsymmetrical) keys
            if len(keys - set(["c4_a0", "c4_aa", "d4", "zero"])) > 0:
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

            # okay, looks like we found a better one!
            best = scheme

        print(degree, best.name)
        # print(best)
    return


if __name__ == "__main__":
    test_get_good_scheme()

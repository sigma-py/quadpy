import numpy
import orthopy
import pytest

import quadpy


@pytest.mark.parametrize("scheme", quadpy.s2.schemes.values())
def test_scheme(scheme):
    try:
        scheme = scheme()
    except TypeError:
        scheme = scheme(1)

    assert scheme.points.dtype == numpy.float64, scheme.name
    assert scheme.weights.dtype == numpy.float64, scheme.name

    print(scheme)

    evaluator = orthopy.s2.xu.Eval(scheme.points, "normal")

    k = 0
    max_err = 0.0
    while True:
        approximate = scheme.integrate(lambda x: next(evaluator), [0.0, 0.0], 1.0)
        exact = evaluator.int_p0 if k == 0 else 0.0
        err = numpy.abs(approximate - exact)
        max_err = max(max_err, numpy.max(err))
        if numpy.any(err > scheme.test_tolerance * 1.1):
            break
        k += 1

    if k - 1 != scheme.degree:
        # find the max error across all polynomials
        for i in range(k + 1, scheme.degree + 1):
            approximate = scheme.integrate(lambda x: next(evaluator), [0.0, 0.0], 1.0)
            exact = evaluator.int_p0 if i == 0 else 0.0
            err = numpy.abs(approximate - exact)
            max_err = max(max_err, numpy.max(err))

        raise AssertionError(
            f"{scheme.name} -- observed: {k - 1}, expected: {scheme.degree} "
            f"(max err: {max_err:.3e})"
        )


@pytest.mark.parametrize("scheme", [quadpy.s2.schemes["lether"](3)])
def test_show(scheme):
    scheme.show()


@pytest.mark.skip()
def test_get_good_scheme():
    for degree in range(51):
        best = None
        for scheme in quadpy.s2.schemes.values():
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
            if numpy.any(scheme.points[0] ** 2 + scheme.points[1] ** 2 > 1):
                continue

            if scheme.test_tolerance > 1.0e-13:
                continue

            keys = set(scheme.symmetry_data.keys())
            if "plain" in keys:
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
                    else:  # years are equal
                        pass

        print(degree, best.name if best is not None else None)

        # print(best)
    return


if __name__ == "__main__":
    test_get_good_scheme()

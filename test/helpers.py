import numpy

from quadpy.helpers import get_all_exponents


def check_degree(quadrature, exact, dim, max_degree, tol):
    exponents = get_all_exponents(dim, max_degree)
    # flatten list
    exponents = numpy.array([item for sublist in exponents for item in sublist])

    flt = numpy.vectorize(float)
    exact_vals = flt([exact(k) for k in exponents])

    def evaluate_all_monomials(x):
        # Evaluate monomials.
        # There's a more complex, faster implementation using matmul, exp, log.
        # However, this only works for strictly positive `x`, and requires some
        # tinkering. See below and
        # <https://stackoverflow.com/a/45421128/353337>.
        return numpy.prod(x[..., None] ** exponents.T[:, None], axis=0).T

    vals = quadrature(evaluate_all_monomials)

    # check relative error
    err = abs(exact_vals - vals)
    is_smaller = err < (1 + abs(exact_vals)) * tol

    if numpy.all(is_smaller):
        return max_degree, numpy.max(err / (1 + abs(exact_vals)))

    k = numpy.where(~is_smaller)[0]
    # Return the max error for all exponents that are one smaller than the max_degree.
    # This is because this functions is usually called with target_degree + 1.
    idx = numpy.sum(exponents, axis=1) < max_degree
    return (
        numpy.sum(exponents[k[0]]) - 1,
        numpy.max(err[idx] / (1 + abs(exact_vals[idx]))),
    )


def find_equal(schemes):
    tol = 1.0e-13
    n = len(schemes)
    for i in range(n):
        found_equal = False
        for j in range(n):
            if schemes[i].name == schemes[j].name:
                continue
            if len(schemes[i].weights) != len(schemes[j].weights):
                continue
            # Check if the point sets are equal
            x = numpy.vstack([schemes[i].weights, schemes[i].points])
            y = numpy.vstack([schemes[j].weights, schemes[j].points])
            is_equal = True
            for x_i in x:
                diff = y - x_i
                diff = numpy.min(numpy.sum(diff ** 2, axis=-1))
                if diff > tol:
                    is_equal = False
                    break
            if is_equal:
                found_equal = True
                a = f"'{schemes[i].name}'"
                try:
                    a += f" ({schemes[i].citation.year})"
                except AttributeError:
                    pass
                b = f"'{schemes[j].name}'"
                try:
                    b += f" ({schemes[j].citation.year})"
                except AttributeError:
                    pass
                print(f"Schemes {a} and {b} are equal.")
        if found_equal:
            print()


def find_best_scheme(schemes, degree, is_points_okay, is_symmetries_okay):
    best = None
    for scheme in schemes:
        try:
            scheme = scheme()  # initialize
        except TypeError:
            continue

        # filter schemes for eligibility
        if scheme.degree < degree:
            continue

        # allow only positive weights
        if any(scheme.weights < 0):
            continue

        # disallow points outside of the domain
        if not is_points_okay(scheme.points):
            continue

        if scheme.test_tolerance > 1.0e-13:
            continue

        # TODO force symmetry data for all schemes
        try:
            keys = set(scheme.symmetry_data.keys())
        except AttributeError:
            continue

        # filter out disallowed (unsymmetrical) keys
        if not is_symmetries_okay(keys):
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

        # okay, looks like we found a better one!
        best = scheme

    return best

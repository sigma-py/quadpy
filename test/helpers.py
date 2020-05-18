import numpy

from quadpy.helpers import get_all_exponents


def check_degree_1d(quadrature, exact, max_degree, tol=1.0e-14):
    val = quadrature(
        lambda x: [x ** degree for degree in range(max_degree + 1)]
    ).flatten()
    exact_val = numpy.array([exact(degree) for degree in range(max_degree + 1)])
    eps = numpy.finfo(float).eps
    # check relative error
    # Allow 1e1 over machine precision.
    alpha = abs(exact_val) * tol + (1e1 + tol + exact_val) * eps
    # check where the error is larger than alpha
    is_larger = (exact_val - val) > alpha
    return numpy.where(is_larger)[0] - 1 if any(is_larger) else max_degree


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
    mytol = (1 + abs(exact_vals)) * tol
    err = abs(exact_vals - vals)
    is_smaller = abs(exact_vals - vals) < mytol

    if numpy.all(is_smaller):
        return max_degree, numpy.max(err)

    k = numpy.where(numpy.logical_not(is_smaller))[0]
    return numpy.sum(exponents[k[0]]) - 1, numpy.max(err)


def check_degree_ortho(approximate, exact, abs_tol=1.0e-14):
    # get the maximum error in all levels execpt the last
    max_err = 0.0
    for k in range(len(approximate) - 1):
        ex = exact[k]
        approx = approximate[k]
        err = numpy.max(numpy.abs(ex - approx))
        if err > max_err:
            max_err = err

    # check absolute error
    for degree, (approx, ex) in enumerate(zip(approximate, exact)):
        err = abs(ex - approx)
        if not numpy.all(err < abs_tol):
            return degree - 1, max_err

    # All values are equal; the degree is at least this.
    return len(approximate), max_err


def find_equal(schemes):
    tol = 1.0e-13
    n = len(schemes)
    for i in range(n):
        found_equal = False
        for j in range(n):
            if schemes[i].name == schemes[j].name:
                continue
            if len(schemes[i].points) != len(schemes[j].points):
                continue
            # Check if the point sets are equal
            x = numpy.column_stack([schemes[i].weights, schemes[i].points])
            y = numpy.column_stack([schemes[j].weights, schemes[j].points])
            is_equal = True
            for x_i in x:
                diff = y - x_i
                diff = numpy.min(numpy.sum(diff ** 2, axis=-1))
                if diff > tol:
                    is_equal = False
                    break
            if is_equal:
                found_equal = True
                a = "'{}'".format(schemes[i].name)
                try:
                    a += " ({})".format(schemes[i].citation.year)
                except AttributeError:
                    pass
                b = "'{}'".format(schemes[j].name)
                try:
                    b += " ({})".format(schemes[j].citation.year)
                except AttributeError:
                    pass
                print(f"Schemes {a} and {b} are equal.")
        if found_equal:
            print()
    return

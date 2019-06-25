# -*- coding: utf-8 -*-
#
import numpy

import orthopy
import quadpy
from quadpy.nsimplex.helpers import integrate_monomial_over_unit_simplex
from quadpy.triangle.helpers import _rot_ab, _s3, _s21, _s111ab


def with_monomials(degree):
    exponents = numpy.concatenate(
        [quadpy.helpers.partition(d, 2) for d in range(degree + 1)]
    )

    exact_vals = numpy.array(
        [integrate_monomial_over_unit_simplex(k) for k in exponents]
    )

    def fun(x):
        k = exponents.T
        # <https://stackoverflow.com/a/46689653/353337>
        s = x.shape[1:] + k.shape[1:]
        return (
            (x.reshape(x.shape[0], -1, 1) ** k.reshape(k.shape[0], 1, -1))
            .prod(0)
            .reshape(s)
        )

    return fun, exact_vals


def triangle_compute_weights(bary, degree):
    """Computes weights from points for a given scheme. Useful for cross-checking scheme
    integrity.
    """
    out = orthopy.triangle.tree(bary, degree, "normal")
    A = numpy.vstack(out)
    exact_vals = numpy.zeros(len(A))
    exact_vals[0] = numpy.sqrt(2) / 2
    return numpy.linalg.lstsq(A, exact_vals, rcond=None)


def quad_compute_weights(points, degree):
    out = orthopy.quadrilateral.tree(points, degree)
    A = numpy.vstack(out)
    exact_vals = numpy.zeros(len(A))
    exact_vals[0] = 4.0
    return numpy.linalg.lstsq(A, exact_vals, rcond=None)


if __name__ == "__main__":
    # scheme = quadpy.triangle.SevenPoint()
    # # scheme = quadpy.triangle.Papanicolopulos("rot", 8)
    # factor = 2  # depends on the convention of the scheme
    # x, res, rank, sv = triangle_compute_weights(scheme.bary.T, scheme.degree)

    scheme = quadpy.quadrilateral.Schmid(2)
    # scheme = quadpy.triangle.Papanicolopulos("rot", 8)
    factor = 2  # depends on the convention of the scheme
    x, res, rank, sv = quad_compute_weights(scheme.points.T, scheme.degree)

    print("num unknowns: {}".format(len(x)))
    print("rank A: {}".format(rank))
    print("residual: {}".format(res))
    assert abs(res[0]) < 1.0e-14
    print("singular values:")
    for val in sv:
        print("  {:.15e}".format(val))

    print()
    print("solution:")
    for val in x:
        print("  {:.15e}".format(factor * val))

    print()
    print("hardcoded weights:")
    for val in scheme.weights:
        print("  {:.15e}".format(val))

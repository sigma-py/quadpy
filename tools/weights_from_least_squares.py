import numpy

import orthopy
import quadpy
from quadpy.enr._helpers import integrate_monomial_over_enr

# from quadpy.nsimplex.helpers import integrate_monomial_over_unit_simplex
# from quadpy.triangle.helpers import _rot_ab, _s3, _s21, _s111ab


def partition(boxes, balls):
    """Create all nonnegative tuples of length d which sum up to n.
    """
    # <https://stackoverflow.com/a/36748940/353337>
    # See <https://stackoverflow.com/a/45348441/353337> for an alterative
    # solution.
    def rec(boxes, balls, parent=tuple()):
        if boxes > 1:
            for i in range(balls + 1):
                yield from rec(boxes - 1, i, parent + (balls - i,))
        else:
            yield parent + (balls,)

    return list(rec(boxes, balls))


# def simplex_monomials(degree):
#     exponents = numpy.concatenate(
#         [quadpy.helpers.partition(d, 2) for d in range(degree + 1)]
#     )
#
#     exact_vals = numpy.array(
#         [integrate_monomial_over_unit_simplex(k) for k in exponents]
#     )
#
#     def fun(x):
#         k = exponents.T
#         # <https://stackoverflow.com/a/46689653/353337>
#         s = x.shape[1:] + k.shape[1:]
#         return (
#             (x.reshape(x.shape[0], -1, 1) ** k.reshape(k.shape[0], 1, -1))
#             .prod(0)
#             .reshape(s)
#         )
#
#     return fun, exact_vals


def e2r_compute_weights(points, degree):
    """Using monomials.
    """
    exponents = numpy.concatenate([partition(2, d) for d in range(degree + 1)])

    exact_vals = numpy.array([integrate_monomial_over_enr(k) for k in exponents])

    def fun(x):
        k = exponents.T
        # <https://stackoverflow.com/a/46689653/353337>
        s = x.shape[1:] + k.shape[1:]
        return (
            (x.reshape(x.shape[0], -1, 1) ** k.reshape(k.shape[0], 1, -1))
            .prod(0)
            .reshape(s)
        )

    A = fun(points).T
    out = numpy.linalg.lstsq(A, exact_vals, rcond=None)
    res = numpy.dot(A, out[0]) - exact_vals
    return out, res


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

    # scheme = quadpy.quadrilateral.Schmid(2)
    # scheme = quadpy.triangle.Papanicolopulos("rot", 8)
    # factor = 2  # depends on the convention of the scheme
    # x, res, rank, sv = quad_compute_weights(scheme.points.T, scheme.degree)

    scheme = quadpy.e2r.rabinowitz_richter_4()
    factor = 1
    (x, _, rank, sv), res = e2r_compute_weights(scheme.points.T, 13)

    print("num unknowns: {}".format(len(x)))
    print(f"rank A: {rank}")
    res_norm = numpy.sqrt(numpy.dot(res, res))
    print(f"res norm: {res_norm}")
    assert res_norm < 1.0e-14
    print("singular values:")
    for val in sv:
        print(f"  {val:.15e}")

    print()
    print("solution:")
    for val in x:
        print("  {:.15e}".format(factor * val))

    print()
    print("hardcoded weights:")
    for val in scheme.weights:
        print(f"  {val:.15e}")

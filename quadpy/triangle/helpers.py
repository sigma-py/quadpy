# -*- coding: utf-8 -*-
#
import numpy
import orthopy
from sympy import Rational

from .. import helpers


def _s3():
    return numpy.full((1, 3), Rational(1, 3))


def _s21(a):
    a = numpy.array(a)
    b = 1 - 2*a
    return numpy.array([
        [a, a, b],
        [a, b, a],
        [b, a, a],
        ])


def _s111(a, b):
    c = 1 - a - b
    return numpy.array([
        [a, b, c],
        [c, a, b],
        [b, c, a],
        [b, a, c],
        [c, b, a],
        [a, c, b],
        ])


def _s111ab(a, b):
    c = 1 - a - b
    out = numpy.array([
        [a, b, c],
        [c, a, b],
        [b, c, a],
        [b, a, c],
        [c, b, a],
        [a, c, b],
        ])
    out = numpy.swapaxes(out, 0, 1)
    return out


def _rot(a, b):
    c = 1 - a - b
    return numpy.array([
        [a, b, c],
        [c, a, b],
        [b, c, a],
        ])


def _rot_ab(a, b):
    c = 1 - a - b
    out = numpy.array([
        [a, b, c],
        [c, a, b],
        [b, c, a],
        ])
    out = numpy.swapaxes(out, 0, 1)
    return out


def _collapse0(a):
    '''Collapse all dimensions of `a` except the first.
    '''
    return numpy.reshape(a, (a.shape[0], numpy.prod(a.shape[1:])))


def untangle2(data):
    bary = []
    weights = []

    if 's3' in data:
        data['s3'] = numpy.array(data['s3']).T
        bary.append(_s3().T)
        weights.append(numpy.tile(data['s3'][0], 1))

    if 's2' in data:
        data['s2'] = numpy.array(data['s2']).T
        s2_data = _s21(data['s2'][1])
        bary.append(_collapse0(s2_data))
        weights.append(numpy.tile(data['s2'][0], 3))

    if 's1' in data:
        data['s1'] = numpy.array(data['s1']).T
        s1_data = _s111ab(*data['s1'][1:])
        bary.append(_collapse0(s1_data))
        weights.append(numpy.tile(data['s1'][0], 6))

    if 'rot' in data:
        data['rot'] = numpy.array(data['rot']).T
        rot_data = _rot_ab(*data['rot'][1:])
        bary.append(_collapse0(rot_data))
        weights.append(numpy.tile(data['rot'][0], 3))

    bary = numpy.column_stack(bary).T
    weights = numpy.concatenate(weights)
    return bary, weights


def untangle3(point_data, weight_data):
    bary = []
    weights = []

    k = 0

    if 's3' in point_data:
        point_data['s3'] = numpy.array(point_data['s3']).T
        bary.append(_s3().T)
        n = point_data['s3'].shape[1]
        weights.append(numpy.tile(weight_data[k:k+n], 1))
        k += n

    if 's2' in point_data:
        point_data['s2'] = numpy.array(point_data['s2']).T
        s2_data = _s21(point_data['s2'][0])
        bary.append(_collapse0(s2_data))
        n = point_data['s2'].shape[1]
        weights.append(numpy.tile(weight_data[k:k+n], 3))
        k += n

    if 's1' in point_data:
        point_data['s1'] = numpy.array(point_data['s1']).T
        s1_data = _s111ab(*point_data['s1'])
        bary.append(_collapse0(s1_data))
        n = point_data['s1'].shape[1]
        weights.append(numpy.tile(weight_data[k:k+n], 6))
        k += n

    if 'rot' in point_data:
        point_data['rot'] = numpy.array(point_data['rot']).T
        rot_data = _rot_ab(*point_data['rot'])
        bary.append(_collapse0(rot_data))
        n = point_data['rot'].shape[1]
        weights.append(numpy.tile(weight_data[k:k+n], 3))
        k += n

    bary = numpy.column_stack(bary).T
    weights = numpy.concatenate(weights)

    return bary, weights


def weights_from_points(point_data, degree):
    '''In a quadrature scheme, the weights really only depend on the points and
    the information up to which degree the scheme is supposed to be exact. This
    method solves a least-squares problem with polynomials orthogonal on the
    triangle, leading to a well-conditioned system.
    '''
    # from quadpy.nsimplex.helpers import integrate_monomial_over_unit_simplex
    # import scipy.special

    # s = quadpy.triangle.Dunavant(20)
    # s = quadpy.triangle.Cubtri()

    # WITH MONOMIALS
    # ==============
    # exponents = numpy.concatenate([
    #     quadpy.helpers.partition(d, 2)
    #     for d in range(s.degree+1)
    #     ])
    #
    # exact_vals = numpy.array([
    #     integrate_monomial_over_unit_simplex(k) for k in exponents
    #     ])
    #
    #
    # def fun(x):
    #     k = exponents.T
    #     # <https://stackoverflow.com/a/46689653/353337>
    #     s = x.shape[1:] + k.shape[1:]
    #     return (
    #         x.reshape(x.shape[0], -1, 1)**k.reshape(k.shape[0], 1, -1)
    #         ).prod(0).reshape(s)

    # WITH ORTHOGONAL POLYNOMIALS
    # ===========================
    exponents = numpy.concatenate([
        helpers.partition(d, 2)
        for d in range(degree+1)
        ])
    # The exact integrals of the orthogonal polynomials over the triangle are
    # 0, except for the one with degree 0 for which we have sqrt(2)/2.
    exact_vals = numpy.zeros(len(exponents))
    exact_vals[0] = numpy.sqrt(2) / 2

    def eval_orthpolys(x):
        '''Evaluate all orthogonal polynomials at x.
        See, e.g.,

        S.-A. Papanicolopulos,
        New fully symmetric and rotationally symmetric cubature rules on the
        triangle using minimal orthonormal bases,
        <https://arxiv.org/pdf/1411.5631.pdf>.
        '''
        def f(i, j, x, y):
            p0, a1, b1, c1 = \
                orthopy.recurrence_coefficients.jacobi(i, 0, 0, 'monic')
            val1 = orthopy.tools.evaluate_orthogonal_polynomial(
                    (x-y)/(x+y), p0, a1, b1, c1
                    )
            # val1 = numpy.polyval(scipy.special.jacobi(i, 0, 0), (x-y)/(x+y))

            p0, a2, b2, c2 = orthopy.recurrence_coefficients.jacobi(
                    j, 2*i+1, 0, 'monic'
                    )
            val2 = orthopy.tools.evaluate_orthogonal_polynomial(
                    1-2*(x+y), p0, a1, b2, c2
                    )
            # val2 = \
            #     numpy.polyval(scipy.special.jacobi(j, 2*i+1, 0), 1-2*(x+y))

            return (
                numpy.sqrt(2*i + 1) * val1 * (x+y)**i
                * numpy.sqrt(2*j + 2*i + 2) * val2
                )

        return numpy.array([
            f(ij[0], ij[1], x[0], x[1])
            for ij in exponents
            ]).T

    a_data = []
    if 's3' in point_data:
        a_data.append(eval_orthpolys(_s3().T[1:]))

    if 's2' in point_data:
        point_data['s2'] = numpy.array(point_data['s2'])
        s2_data = _s21(*point_data['s2'].T)
        a_data.append(numpy.sum(eval_orthpolys(s2_data[1:]), axis=1))

    if 's1' in point_data:
        point_data['s1'] = numpy.array(point_data['s1'])
        s1_data = _s111ab(*point_data['s1'].T)
        a_data.append(numpy.sum(eval_orthpolys(s1_data[1:]), axis=1))

    A = numpy.concatenate(a_data).T

    x, res, rank, sv = numpy.linalg.lstsq(A, exact_vals)
    assert numpy.all(abs(res) < 1.0e-15)
    return 2*x

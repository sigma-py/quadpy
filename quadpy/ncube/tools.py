# -*- coding: utf-8 -*-
#
import numpy

from .. import helpers


def transform(xi, cube):
    '''Transform the points `xi` from the reference cube to `cube`.
    '''
    # For d==2, the result used to be computed with
    #
    # out = (
    #     + outer(0.25*(1.0-xi[0])*(1.0-xi[1]), cube[0, 0])
    #     + outer(0.25*(1.0+xi[0])*(1.0-xi[1]), cube[1, 0])
    #     + outer(0.25*(1.0-xi[0])*(1.0+xi[1]), cube[0, 1])
    #     + outer(0.25*(1.0+xi[0])*(1.0+xi[1]), cube[1, 1])
    #     )
    #
    # This array of multiplications and additions is reminiscent of dot(), and
    # indeed tensordot() can handle the situation. We just need to compute the
    # `1+-xi` products and align them with `cube`.
    one_mp_xi = numpy.stack([
        0.5 * (1.0 - xi),
        0.5 * (1.0 + xi),
        ], axis=1)
    a = helpers.n_outer(one_mp_xi)

    # TODO kahan tensordot
    # <https://stackoverflow.com/q/45372098/353337>
    d = xi.shape[0]
    return numpy.tensordot(a, cube, axes=(range(d), range(d)))


def get_detJ(xi, cube):
    '''Get the determinant of the transformation matrix.
    '''
    # For d==2, the result can be computed with
    # ```
    # J0 = (
    #     - numpy.multiply.outer(0.25*(1-xi[1]), quad[0, 0])
    #     + numpy.multiply.outer(0.25*(1-xi[1]), quad[1, 0])
    #     - numpy.multiply.outer(0.25*(1+xi[1]), quad[0, 1])
    #     + numpy.multiply.outer(0.25*(1+xi[1]), quad[1, 1])
    #     ).T
    # J1 = (
    #     - numpy.multiply.outer(0.25*(1-xi[0]), quad[0, 0])
    #     - numpy.multiply.outer(0.25*(1+xi[0]), quad[1, 0])
    #     + numpy.multiply.outer(0.25*(1-xi[0]), quad[0, 1])
    #     + numpy.multiply.outer(0.25*(1+xi[0]), quad[1, 1])
    #     ).T
    # out = J0[0]*J1[1] - J1[0]*J0[1]
    # ```
    # Like transform(), simplify here and form the determinant explicitly.
    d = xi.shape[0]

    one_mp_xi = numpy.stack([
        0.5 * (1.0 - xi),
        0.5 * (1.0 + xi),
        ], axis=1)

    # Build the Jacobi matrix row by row.
    J = []
    for k in range(d):
        a = one_mp_xi.copy()
        a[k, 0, :] = -0.5
        a[k, 1, :] = +0.5
        a0 = helpers.n_outer(a)
        J.append(numpy.tensordot(a0, cube, axes=(range(d), range(d))).T)

    # `det` needs the square at the end. Fortran...
    # For d==2 or d==3, we could avoid this copy and compute the determinant
    # with their elementary formulas, i.e.,
    #
    #     + J[0][0]*J[1][1] - J[1][0]*J[0][1];
    #
    #     + J0[0]*J1[1]*J2[2] + J1[0]*J2[1]*J0[2] + J2[0]*J0[1]*J1[2]
    #     - J0[2]*J1[1]*J2[0] - J1[2]*J2[1]*J0[0] - J2[2]*J0[1]*J1[0].
    #
    J = numpy.array(J)
    J = numpy.moveaxis(J, (0, 1), (-2, -1))
    out = numpy.linalg.det(J)
    return out


def integrate(f, ncube, scheme, dot=numpy.dot):
    flt = numpy.vectorize(float, otypes=[float])
    x = transform(flt(scheme.points).T, ncube).T
    detJ = get_detJ(flt(scheme.points).T, ncube)
    return dot(f(x)*abs(detJ), flt(scheme.weights))


def ncube_points(*xyz):
    '''Given the end points of an n-cube aligned with the coordinate axes, this
    returns the corner points of the cube in the correct data structure.
    '''
    return numpy.moveaxis(
            numpy.array(numpy.meshgrid(*xyz, indexing='ij')),
            0, -1
            )

# -*- coding: utf-8 -*-
#
import numpy

from .. import helpers


def transform(xi, cube):
    '''Transform the points xi from the reference cube to the cube `cube`.
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

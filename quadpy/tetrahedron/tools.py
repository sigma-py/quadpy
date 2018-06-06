# -*- coding: utf-8 -*-
#
import numpy

from .. import helpers
from ..nsimplex import transform, get_vol


def show(
    scheme,
    tet=numpy.array(
        [
            [+1, 0, -1.0 / numpy.sqrt(2.0)],
            [-1, 0, -1.0 / numpy.sqrt(2.0)],
            [0, +1, +1.0 / numpy.sqrt(2.0)],
            [0, -1, +1.0 / numpy.sqrt(2.0)],
        ]
    ),
    backend="mpl",
):
    edges = numpy.array([[tet[i], tet[j]] for i in range(4) for j in range(i)])
    edges = numpy.moveaxis(edges, 1, 2)
    helpers.backend_to_function[backend](
        transform(scheme.points.T, tet.T).T, scheme.weights, get_vol(tet), edges
    )
    return

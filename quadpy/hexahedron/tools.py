# -*- coding: utf-8 -*-
#
import numpy

from .. import helpers
from ..ncube import transform, integrate
from ..ncube import ncube_points as cube_points


def show(
        scheme,
        hexa=cube_points([0.0, 1.0], [0.0, 1.0], [0.0, 1.0]),
        backend='mpl'
        ):
    '''Shows the quadrature points on a given hexahedron. The size of the
    balls around the points coincides with their weights.
    '''
    backend_to_function = {
        'mayavi': helpers.show_mayavi,
        'mpl': helpers.show_mpl,
        'vtk': helpers.show_vtk,
        }
    backend_to_function[backend](
            transform(scheme.points.T, hexa),
            scheme.weights,
            integrate(lambda x: 1.0, hexa, scheme),
            _get_edges(hexa)
            )
    return


def _get_edges(hexa):
    edges = numpy.array([
        [hexa[0, 0, 0], hexa[1, 0, 0]],
        [hexa[1, 0, 0], hexa[1, 1, 0]],
        [hexa[1, 1, 0], hexa[0, 1, 0]],
        [hexa[0, 1, 0], hexa[0, 0, 0]],
        #
        [hexa[0, 0, 1], hexa[1, 0, 1]],
        [hexa[1, 0, 1], hexa[1, 1, 1]],
        [hexa[1, 1, 1], hexa[0, 1, 1]],
        [hexa[0, 1, 1], hexa[0, 0, 1]],
        #
        [hexa[0, 0, 0], hexa[0, 0, 1]],
        [hexa[1, 0, 0], hexa[1, 0, 1]],
        [hexa[1, 1, 0], hexa[1, 1, 1]],
        [hexa[0, 1, 0], hexa[0, 1, 1]],
        ])
    edges = numpy.moveaxis(edges, 1, 2)
    return edges

# -*- coding: utf-8 -*-
#
import numpy

from .. import helpers
from ..ncube import NCubeScheme
from ..ncube import ncube_points as rectangle_points
from ..ncube import transform
from ..nsimplex import get_vol


class QuadrilateralScheme(NCubeScheme):
    def __init__(self, name, weights, points, degree, citation=None):
        self.name = name
        self.weights = weights
        self.points = points
        self.degree = degree
        self.citation = citation
        return

    def show(self, *args, **kwargs):
        import matplotlib.pyplot as plt

        self.plot(*args, **kwargs)
        plt.show()
        return

    def plot(self, quad=rectangle_points([0.0, 1.0], [0.0, 1.0]), show_axes=False):
        """Shows the quadrature points on a given quad. The area of the disks
        around the points coincides with their weights.
        """
        import matplotlib.pyplot as plt

        plt.plot(quad[0][0], quad[1][0], "-k")
        plt.plot(quad[1][0], quad[1][1], "-k")
        plt.plot(quad[1][1], quad[0][1], "-k")
        plt.plot(quad[0][1], quad[0][0], "-k")

        plt.axis("equal")

        if not show_axes:
            plt.gca().set_axis_off()

        transformed_pts = transform(self.points.T, quad)

        # compute volume by splitting it in two triangles
        vol = get_vol(numpy.array([quad[0][0], quad[1][0], quad[0][1]])) + get_vol(
            numpy.array([quad[0][0], quad[0][1], quad[1][1]])
        )
        helpers.plot_disks(plt, transformed_pts, self.weights, vol)
        return


def zero(weight):
    return numpy.array([weight]), numpy.array([[0, 0]])


def pmx(*data):
    w, x, y = numpy.array(data).T
    points = _stack_first_last([[+x, y], [-x, y]])
    weights = numpy.tile(w, 2)
    return weights, points


def pmy(*data):
    w, x, y = numpy.array(data).T
    points = _stack_first_last([[x, +y], [x, -y]])
    weights = numpy.tile(w, 2)
    return weights, points


def pm(*data):
    w, s, t = numpy.array(data).T
    points = _stack_first_last([[+s, +t], [-s, -t]])
    weights = numpy.tile(w, 2)
    return weights, points


def pm2(*data):
    w, x, y = numpy.array(data).T
    points = _stack_first_last([[+x, +y], [+x, -y], [-x, +y], [-x, -y]])
    weights = numpy.tile(w, 4)
    return weights, points


def symm_r0(*data):
    data = numpy.array(data)
    w, r = data.T
    zero = numpy.zeros(w.shape[0], dtype=r.dtype)
    points = _stack_first_last([[+r, zero], [-r, zero], [zero, +r], [zero, -r]])
    weights = numpy.tile(w, 4)
    return weights, points


def symm_s(*data):
    w, s = numpy.array(data).T
    points = _stack_first_last([[+s, +s], [+s, -s], [-s, +s], [-s, -s]])
    weights = numpy.tile(w, 4)
    return weights, points


def symm_s_t(*data):
    w, s, t = numpy.array(data).T
    points = _stack_first_last(
        [[+s, +t], [-s, +t], [+s, -t], [-s, -t], [+t, +s], [-t, +s], [+t, -s], [-t, -s]]
    )
    weights = numpy.tile(w, 8)
    return weights, points


def s4(*data):
    w, a, b = numpy.array(data).T
    points = _stack_first_last([[+a, +b], [-a, -b], [-b, +a], [+b, -a]])
    weights = numpy.tile(w, 4)
    return weights, points


def s4a(*data):
    w, a = numpy.array(data).T
    points = _stack_first_last([[+a, +a], [+a, -a], [-a, +a], [-a, -a]])
    weights = numpy.tile(w, 4)
    return weights, points


def _stack_first_last(arr):
    """Stacks an input array of shape (i, j, k) such that the output array is of shape
    (i*k, j).
    """
    arr = numpy.swapaxes(arr, 0, 1)
    return arr.reshape(arr.shape[0], -1).T


def concat(*data):
    weights = numpy.concatenate([t[0] for t in data])
    points = numpy.vstack([t[1] for t in data])
    return weights, points

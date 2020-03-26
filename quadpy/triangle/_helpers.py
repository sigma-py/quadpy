import numpy
import sympy

from ..helpers import plot_disks
from ..nsimplex import NSimplexScheme, get_vol, transform


class TriangleScheme(NSimplexScheme):
    def __init__(self, name, weights, points, degree, citation=None):
        self.name = name
        self.degree = degree
        self.citation = citation

        if weights.dtype == numpy.float64:
            self.weights = weights
        else:
            assert weights.dtype in [numpy.dtype("O"), numpy.int_]
            self.weights = weights.astype(numpy.float64)
            self.weights_symbolic = weights

        if points.dtype == numpy.float64:
            self.points = points
        else:
            assert points.dtype in [numpy.dtype("O"), numpy.int_]
            self.points = points.astype(numpy.float64)
            self.points_symbolic = points
        return

    def show(self, *args, **kwargs):
        import matplotlib.pyplot as plt

        self.plot(*args, **kwargs)
        plt.show()
        return

    def plot(
        self,
        triangle=numpy.array([[-0.5, 0.0], [+0.5, 0.0], [0, 0.5 * (numpy.sqrt(3))]]),
        show_axes=False,
    ):
        """Shows the quadrature points on a given triangle. The size of the circles
        around the points coincides with their weights.
        """
        import matplotlib.pyplot as plt

        plt.plot(triangle[:, 0], triangle[:, 1], "-k")
        plt.plot(
            [triangle[-1, 0], triangle[0, 0]], [triangle[-1, 1], triangle[0, 1]], "-k"
        )

        if not show_axes:
            plt.gca().set_axis_off()

        transformed_pts = transform(self.points.T, triangle.T).T

        vol = get_vol(triangle)
        plot_disks(plt, transformed_pts, self.weights, vol)

        plt.axis("equal")
        return


def _s3(symbolic):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    return numpy.full((1, 3), frac(1, 3))


def _s21(a):
    a = numpy.array(a)
    b = 1 - 2 * a
    return numpy.array([[a, a, b], [a, b, a], [b, a, a]])


def _s111ab(a, b):
    c = 1 - a - b
    out = numpy.array(
        [[a, b, c], [c, a, b], [b, c, a], [b, a, c], [c, b, a], [a, c, b]]
    )
    out = numpy.swapaxes(out, 0, 1)
    return out


def _rot_ab(a, b):
    c = 1 - a - b
    out = numpy.array([[a, b, c], [c, a, b], [b, c, a]])
    out = numpy.swapaxes(out, 0, 1)
    return out


def _collapse0(a):
    """Collapse all dimensions of `a` except the first.
    """
    return a.reshape(a.shape[0], -1)


def untangle2(data, symbolic=False):
    points = []
    weights = []

    if "s3" in data:
        d = numpy.array(data["s3"]).T
        points.append(_s3(symbolic).T)
        weights.append(numpy.tile(d[0], 1))

    if "s2" in data:
        d = numpy.array(data["s2"]).T
        s2_data = _s21(d[1])
        points.append(_collapse0(s2_data))
        weights.append(numpy.tile(d[0], 3))

    if "s1" in data:
        d = numpy.array(data["s1"]).T
        s1_data = _s111ab(*d[1:])
        points.append(_collapse0(s1_data))
        weights.append(numpy.tile(d[0], 6))

    if "rot" in data:
        d = numpy.array(data["rot"]).T
        rot_data = _rot_ab(*d[1:])
        points.append(_collapse0(rot_data))
        weights.append(numpy.tile(d[0], 3))

    points = numpy.column_stack(points).T
    weights = numpy.concatenate(weights)
    return points, weights


def s3(weight):
    symbolic = isinstance(weight, float)
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    return numpy.array([weight]), numpy.full((1, 3), frac(1, 3))


def s2(*data):
    w, a = numpy.array(data).T
    b = 1 - 2 * a
    points = _stack_first_last([[a, a, b], [a, b, a], [b, a, a]])
    weights = numpy.tile(w, 3)
    return weights, points


def s1(*data):
    w, a, b = numpy.array(data).T
    c = 1 - a - b
    points = _stack_first_last(
        [[a, b, c], [c, a, b], [b, c, a], [b, a, c], [c, b, a], [a, c, b]]
    )
    weights = numpy.tile(w, 6)
    return weights, points


def r(*data):
    w, r = numpy.array(data).T
    a = r + (1 - r) / 3
    b = (1 - a) / 2
    points = _stack_first_last([[a, b, b], [b, a, b], [b, b, a]])
    weights = numpy.tile(w, 3)
    return weights, points


def rot_ab(*data):
    w, a, b = numpy.array(data).T
    c = 1 - a - b
    points = _stack_first_last([[a, b, c], [c, a, b], [b, c, a]])
    weights = numpy.tile(w, 3)
    return weights, points


def mirror(*data):
    w, a, b = numpy.array(data).T
    c = 1 - a - b
    points = _stack_first_last([[a, b, c], [b, a, c]])
    weights = numpy.tile(w, 2)
    return weights, points


def alpha(*data):
    """From the article Liu-Vinokur:

    mu_i = (1 + (n-1) alpha) / n,
    mu_j = (1 - alpha) / n    for j!=i,

    where n is the number of vertices
    """
    w, alpha = numpy.array(data).T
    a = (1 + 2 * alpha) / 3
    b = (1 - alpha) / 3
    points = _stack_first_last([[a, b, b], [b, a, b], [b, b, a]])
    weights = numpy.tile(w, 3)
    return weights, points


def gamma_delta(*data):
    """From the article Liu-Vinokur:

    mu_i = (1 + (n-1) gamma - delta) / n,
    mu_j = (1 + (n-1) delta - gamma) / n,
    mu_k = (1 - gamma - delta) / n    for k!=i, k!=j,

    where n is the number of vertices
    """
    w, gamma, delta = numpy.array(data).T
    a = (1 + 2 * gamma - delta) / 3
    b = (1 + 2 * delta - gamma) / 3
    c = (1 - gamma - delta) / 3
    points = _stack_first_last(
        [[a, b, c], [c, a, b], [b, c, a], [a, c, b], [b, a, c], [c, b, a]]
    )
    weights = numpy.tile(w, 6)
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

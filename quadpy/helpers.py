# -*- coding: utf-8 -*-
#
import math
import numpy


def n_outer(a):
    '''Given a list (tuple, array) of arrays, this method computes their outer
    product. If the dimension of the input arrays is larger than one, the
    product is formed across the first dimension; all other dimensions must
    coincide in size.

    Examples:
    n_outer([np.ones(4), np.ones(5)]).shape == (4, 5)
    n_outer([np.ones(4), np.ones(5), np.ones(6)]).shape == (4, 5, 6)
    n_outer([np.ones(4, 3, 7), np.ones(5, 3, 7)]).shape == (4, 5, 3, 7)
    '''
    # <https://stackoverflow.com/a/45376730/353337>
    d = len(a)

    # If the elements are more than one-dimensional, assert that the extra
    # dimensions are all equal.
    s0 = a[0].shape
    for arr in a:
        assert s0[1:] == arr.shape[1:]

    out = a[0]
    for k in range(1, d):
        # Basically outer products. Checkout `numpy.outer`'s implementation for
        # comparison.
        out = numpy.multiply(
                # Insert a newaxis after k `:`
                out[(slice(None),) * k + (numpy.newaxis,)],
                # Insert a newaxis at the beginning
                a[k][numpy.newaxis],
                )
    return out


def kahan_sum(a, axis=0):
    '''Kahan summation of the numpy array `a` along axis `axis`.
    '''
    # See <https://en.wikipedia.org/wiki/Kahan_summation_algorithm> for
    # details.
    k = axis % len(a.shape)
    s = numpy.zeros(a.shape[:axis] + a.shape[k+1:])
    c = numpy.zeros(s.shape)
    for i in range(a.shape[axis]):
        # http://stackoverflow.com/a/42817610/353337
        y = a[(slice(None),) * k + (i,)] - c
        t = s + y
        c = (t - s) - y
        s = t.copy()
    return s


def plot_disks(plt, pts, weights, total_area):
    '''Plot a circles at quadrature points according to weights.
    '''
    sum_weights = math.fsum(weights)
    for tp, weight in zip(pts, weights):
        # use matplotlib 2.0's color scheme
        color = '#1f77b4' if weight >= 0 else '#d62728'
        # highlight circle center
        plt.plot(
            [tp[0]], [tp[1]],
            linestyle='None', marker='.', color=color
            )
        # Choose radius such that the sum of areas of the circles equals
        # total_area.
        radius = math.sqrt(abs(weight)/sum_weights * total_area/math.pi)
        circ = plt.Circle((tp[0], tp[1]), radius, color=color, alpha=0.5)
        plt.gca().add_artist(circ)

    return


# pylint: disable=too-many-locals
def plot_spheres(
        plt, ax, pts, weights, total_volume
        ):
    h = 1.0e-2

    sum_weights = math.fsum(weights)
    for tp, weight in zip(pts, weights):
        # Choose radius such that the sum of volumes of the balls equals
        # total_volume.
        r = (
            abs(weight)/sum_weights * total_volume/(4.0/3.0 * numpy.pi)
            )**(1.0/3.0)

        # http://matplotlib.org/examples/mplot3d/surface3d_demo2.html
        # Compute sphere for every point anew. This is more costly on the
        # numerical side, but gives the flexibility of drawing sphere of
        # different size with different number of points. Another options would
        # be to precomoute x, y, z before the loop, but this can be heavy on
        # the graphics output. See
        # <https://stackoverflow.com/q/45324258/353337>.
        u = numpy.linspace(0, 2 * numpy.pi, int(2*numpy.pi/h*r) + 1)
        v = numpy.linspace(0, numpy.pi, int(numpy.pi/h*r) + 1)
        x = numpy.outer(numpy.cos(u), numpy.sin(v))
        y = numpy.outer(numpy.sin(u), numpy.sin(v))
        z = numpy.outer(numpy.ones(numpy.size(u)), numpy.cos(v))

        color = '#1f77b4' if weight >= 0 else '#d62728'
        # highlight ball center
        plt.plot(
            [tp[0]], [tp[1]], [tp[2]],
            linestyle='None', marker='.', color=color
            )

        ax.plot_surface(
            r*x + tp[0], r*y + tp[1], r*z + tp[2],
            color=color,
            alpha=0.3,
            linewidth=1
            )

    ax.set_axis_off()
    return


def partition(balls, boxes):
    '''Create all nonnegative tuples of length d which sum up to n.
    '''
    # <https://stackoverflow.com/a/36748940/353337>
    # See <https://stackoverflow.com/a/45348441/353337> for an alterantive
    # solution.
    def rec(boxes, balls, parent=tuple()):
        if boxes > 1:
            for i in range(balls + 1):
                for x in rec(boxes - 1, i, parent + (balls - i,)):
                    yield x
        else:
            yield parent + (balls,)

    return list(rec(boxes, balls))

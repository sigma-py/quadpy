# -*- coding: utf-8 -*-
#
import matplotlib.pyplot as plt
import numpy

from .. import helpers
from .gauss_kronrod import _gauss_kronrod_integrate


def integrate(f, interval, scheme, sumfun=helpers.kahan_sum):
    xi = scheme.points
    x = \
        + numpy.multiply.outer(0.5 * (1.0 - xi), interval[0]) \
        + numpy.multiply.outer(0.5 * (1.0 + xi), interval[1])
    x = x.T
    diff = interval[1] - interval[0]
    len_intervals = numpy.sqrt(numpy.einsum('...j,...j->...', diff, diff))
    # The factor 0.5 is from the length of the reference line [-1, 1].
    return sumfun(
        numpy.moveaxis(scheme.weights * f(x), -1, 0) * 0.5 * len_intervals
        )


def _numpy_all_except(a, axis=-1):
    axes = numpy.arange(a.ndim)
    axes = numpy.delete(axes, axis)
    return numpy.all(a, axis=tuple(axes))


# pylint: disable=too-many-arguments
# pylint: disable=too-many-locals
def adaptive_integrate(
        f, intervals, eps,
        kronrod_degree=7,
        minimum_interval_length=None,
        sumfun=helpers.kahan_sum
        ):
    intervals = numpy.array(intervals)
    if len(intervals.shape) == 1:
        intervals = intervals[..., None]

    lengths = abs(intervals[1] - intervals[0])
    total_length = sumfun(lengths)

    if minimum_interval_length is None:
        minimum_interval_length = total_length / 2**10

    # Use Gauss-Kronrod 7/15 scheme for error estimation and adaptivity.
    _, val_g, error_estimate = \
        _gauss_kronrod_integrate(kronrod_degree, f, intervals, sumfun=sumfun)

    # Mark intervals with acceptable approximations. For this, take all()
    # across every dimension except the last one, which is the interval index.
    is_good = _numpy_all_except(
            error_estimate < eps * lengths / total_length,
            axis=-1
            )
    # add values from good intervals to sum
    quad_sum = sumfun(val_g[..., is_good], axis=-1)
    global_error_estimate = sumfun(error_estimate[..., is_good], axis=-1)

    is_bad = numpy.logical_not(is_good)
    while any(is_bad):
        # split the bad intervals in half
        intervals = intervals[..., is_bad]
        midpoints = 0.5 * (intervals[0] + intervals[1])
        intervals = numpy.array([
            numpy.concatenate([intervals[0], midpoints]),
            numpy.concatenate([midpoints, intervals[1]]),
            ])
        # compute values and error estimates for the new intervals
        _, val_g, error_estimate = _gauss_kronrod_integrate(
                kronrod_degree, f, intervals, sumfun=sumfun
                )
        # mark good intervals, gather values and error estimates
        lengths = abs(intervals[1] - intervals[0])
        assert all(lengths > minimum_interval_length)
        is_good = _numpy_all_except(
                error_estimate < eps * lengths / total_length,
                axis=-1
                )
        # add values from good intervals to sum
        quad_sum += sumfun(val_g[..., is_good], axis=-1)
        global_error_estimate += sumfun(error_estimate[..., is_good], axis=-1)
        is_bad = numpy.logical_not(is_good)

    return quad_sum, global_error_estimate


def show(*args, **kwargs):
    plot(*args, **kwargs)
    plt.show()
    return


def plot(scheme, interval=numpy.array([[-1.0], [1.0]]), show_axes=False):
    # change default range so that new disks will work
    plt.axis('equal')
    # ax.set_xlim((-1.5, 1.5))
    # ax.set_ylim((-1.5, 1.5))

    if not show_axes:
        plt.gca().set_axis_off()

    plt.plot(interval, [0, 0], color='k')

    pts = numpy.column_stack([scheme.points, numpy.zeros(len(scheme.points))])

    # The total area is used to gauge the disk radii. This is only meaningful
    # for 2D manifolds, not for the circle. What we do instead is choose the
    # total_area such that the sum of the disk radii equals b-a.
    length = interval[1] - interval[0]
    total_area = 0.25 * length**2 * numpy.pi * sum(scheme.weights) \
        / sum(numpy.sqrt(abs(scheme.weights)))**2

    helpers.plot_disks(
        plt, pts, scheme.weights, total_area
        )
    return

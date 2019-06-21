# -*- coding: utf-8 -*-
#
import numpy

from ._gauss_kronrod import _gauss_kronrod_integrate


def _numpy_all_except(a, axis=-1):
    axes = numpy.arange(a.ndim)
    axes = numpy.delete(axes, axis)
    return numpy.all(a, axis=tuple(axes))


class IntegrationError(Exception):
    pass


def integrate_adaptive(
    f, intervals, eps, kronrod_degree=7, minimum_interval_length=0.0, dot=numpy.dot
):
    sumfun = numpy.sum

    intervals = numpy.array(intervals)
    if len(intervals.shape) == 1:
        intervals = intervals[..., None]

    lengths = abs(intervals[1] - intervals[0])
    total_length = sumfun(lengths)

    # Use Gauss-Kronrod 7/15 scheme for error estimation and adaptivity.
    _, val_g, error_estimate = _gauss_kronrod_integrate(
        kronrod_degree, f, intervals, dot=dot
    )

    # Mark intervals with acceptable approximations. For this, take all() across every
    # dimension except the last one, which is the interval index.
    is_good = _numpy_all_except(error_estimate < eps * lengths / total_length, axis=-1)
    # add values from good intervals to sum
    quad_sum = sumfun(val_g[..., is_good], axis=-1)
    global_error_estimate = sumfun(error_estimate[..., is_good], axis=-1)

    is_bad = numpy.logical_not(is_good)
    while any(is_bad):
        # split the bad intervals in half
        intervals = intervals[..., is_bad]
        midpoints = 0.5 * (intervals[0] + intervals[1])
        intervals = numpy.array(
            [
                numpy.concatenate([intervals[0], midpoints]),
                numpy.concatenate([midpoints, intervals[1]]),
            ]
        )
        # compute values and error estimates for the new intervals
        _, val_g, error_estimate = _gauss_kronrod_integrate(
            kronrod_degree, f, intervals, dot=dot
        )
        # mark good intervals, gather values and error estimates
        lengths = abs(intervals[1] - intervals[0])
        if any(lengths < minimum_interval_length):
            raise IntegrationError(
                "Tolerance ({}) could not be reached with the minimum_interval_length (= {}).".format(
                    eps, minimum_interval_length
                )
            )
        is_good = _numpy_all_except(
            error_estimate < eps * lengths / total_length, axis=-1
        )
        # add values from good intervals to sum
        quad_sum += sumfun(val_g[..., is_good], axis=-1)
        global_error_estimate += sumfun(error_estimate[..., is_good], axis=-1)
        is_bad = numpy.logical_not(is_good)

    return quad_sum, global_error_estimate

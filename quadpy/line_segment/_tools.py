import numpy

from ._gauss_kronrod import _gauss_kronrod_integrate


def _numpy_all_except_last(a):
    return numpy.all(a, axis=tuple(range(len(a.shape) - 1)))


class IntegrationError(Exception):
    pass


def integrate_adaptive(
    f,
    intervals,
    eps_abs=1.0e-10,
    eps_rel=1.0e-10,
    kronrod_degree=7,
    minimum_interval_length=0.0,
    max_num_subintervals=numpy.inf,
    dot=numpy.dot,
    domain_shape=None,
    range_shape=None,
):
    intervals = numpy.asarray(intervals)
    assert intervals.shape[0] == 2

    assert (
        eps_abs is not None or eps_rel is not None
    ), "One of eps_abs, eps_rel must be specified."

    # Use Gauss-Kronrod scheme for error estimation and adaptivity.
    # The method also returns guesses for the domain_shape and range_shape (if any of
    # them is None).
    _, val, a, error_estimate, domain_shape, range_shape = _gauss_kronrod_integrate(
        kronrod_degree,
        f,
        intervals,
        dot=dot,
        domain_shape=domain_shape,
        range_shape=range_shape,
    )

    # Flatten the list of intervals so we can do good-bad bookkeeping via a list.
    intervals = intervals.reshape((2,) + domain_shape + (-1,))
    num_subintervals = 1
    a_orig = a.reshape(-1)
    val_shape = val.shape
    val = val.reshape(range_shape + (-1,))
    error_estimate = error_estimate.reshape(range_shape + (-1,))

    total_val = numpy.zeros(val.shape, dtype=val.dtype)
    total_error_estimate = numpy.zeros(error_estimate.shape)

    # Mark intervals with acceptable approximations. For this, take all() across every
    # dimension except the last one (which is the interval index).
    is_good = numpy.ones(error_estimate.shape[-1], dtype=bool)
    if eps_abs is not None:
        is_good = numpy.logical_and(
            is_good, _numpy_all_except_last(error_estimate < eps_abs)
        )
    if eps_rel is not None:
        is_good = numpy.logical_and(
            is_good, _numpy_all_except_last(error_estimate < eps_rel * numpy.abs(val)),
        )
    idx = numpy.arange(intervals.shape[-1])

    total_val[..., is_good] += val[..., is_good]
    total_error_estimate[..., is_good] += error_estimate[..., is_good]

    k = 0
    while numpy.any(~is_good):
        # split the bad intervals in half
        intervals = intervals[..., ~is_good]
        midpoints = 0.5 * (intervals[0] + intervals[1])
        intervals = numpy.array(
            [
                numpy.concatenate([intervals[0], midpoints], axis=-1),
                numpy.concatenate([midpoints, intervals[1]], axis=-1),
            ]
        )
        idx = numpy.concatenate([idx[~is_good], idx[~is_good]])
        num_subintervals += numpy.sum(~is_good)

        if num_subintervals > max_num_subintervals:
            raise IntegrationError(
                f"Tolerances (abs: {eps_abs}, rel: {eps_rel}) could not be reached "
                f"with the given max_num_subintervals (= {max_num_subintervals})."
            )

        # compute values and error estimates for the new intervals
        _, val, a, error_estimate, _, _ = _gauss_kronrod_integrate(
            kronrod_degree,
            f,
            intervals,
            dot=dot,
            domain_shape=domain_shape,
            range_shape=range_shape,
        )

        # relative interval lenghts
        b = a / a_orig

        # mark good intervals, gather values and error estimates
        if numpy.any(a < minimum_interval_length):
            raise IntegrationError(
                f"Tolerances (abs: {eps_abs}, rel: {eps_rel}) could not be reached "
                f"with the given minimum_interval_length (= {minimum_interval_length})."
            )

        # TODO speed up
        # tentative total value (as if all intervals were good)
        ttv = total_val.copy()
        for j, i in enumerate(idx):
            ttv[..., i] += val[..., j]

        is_good = numpy.ones(error_estimate.shape[-1], dtype=bool)
        if eps_abs is not None:
            is_good = numpy.logical_and(
                is_good, _numpy_all_except_last(error_estimate < eps_abs * b),
            )
        if eps_rel is not None:
            is_good = numpy.logical_and(
                is_good,
                _numpy_all_except_last(error_estimate < eps_rel * b * numpy.abs(ttv)),
            )

        # TODO speed up
        for j, i in enumerate(idx[is_good]):
            total_val[..., i] += val[..., is_good][..., j]
            total_error_estimate[..., i] += error_estimate[..., is_good][..., j]
        k += 1

    total_val = total_val.reshape(val_shape)
    total_error_estimate = total_error_estimate.reshape(val_shape)
    return total_val, total_error_estimate

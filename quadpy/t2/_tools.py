import numpy as np

from ..tn import get_vol
from ._dunavant import dunavant_05, dunavant_10


def _numpy_all_except(a, axis=-1):
    axes = np.arange(a.ndim)
    axes = np.delete(axes, axis)
    return np.all(a, axis=tuple(axes))


def integrate_adaptive(
    f,
    triangles,
    eps,
    minimum_triangle_area=None,
    scheme1=dunavant_05(),
    scheme2=dunavant_10(),
    dot=np.dot,
):
    sumfun = np.sum

    triangles = np.array(triangles)
    if len(triangles.shape) == 2:
        # add dimension in the second-to-last place
        triangles = np.expand_dims(triangles, -2)

    areas = get_vol(triangles)
    total_area = sumfun(areas)

    if minimum_triangle_area is None:
        minimum_triangle_area = total_area * 0.25 ** 10

    val1 = scheme1.integrate(f, triangles, dot=dot)
    val2 = scheme2.integrate(f, triangles, dot=dot)
    error_estimate = abs(val1 - val2)

    # Mark intervals with acceptable approximations. For this, take all()
    # across every dimension except the last one, which is the interval index.
    is_good = _numpy_all_except(error_estimate < eps * areas / total_area, axis=-1)

    # add values from good intervals to sum
    quad_sum = sumfun(val1[..., is_good], axis=-1)
    global_error_estimate = sumfun(error_estimate[..., is_good], axis=-1)

    is_bad = np.logical_not(is_good)
    while any(is_bad):
        # split the bad triangles into four #triforce
        #
        #         /\
        #        /__\
        #       /\  /\
        #      /__\/__\
        #
        triangles = triangles[..., is_bad, :]
        midpoints = [
            0.5 * (triangles[1] + triangles[2]),
            0.5 * (triangles[2] + triangles[0]),
            0.5 * (triangles[0] + triangles[1]),
        ]
        triangles = np.array(
            [
                np.concatenate(
                    [triangles[0], triangles[1], triangles[2], midpoints[0]]
                ),
                np.concatenate(
                    [midpoints[1], midpoints[2], midpoints[0], midpoints[1]]
                ),
                np.concatenate(
                    [midpoints[2], midpoints[0], midpoints[1], midpoints[2]]
                ),
            ]
        )
        areas = get_vol(triangles)
        assert all(areas > minimum_triangle_area)

        # compute values and error estimates for the new intervals
        val1 = scheme1.integrate(f, triangles, dot=dot)
        val2 = scheme2.integrate(f, triangles, dot=dot)
        error_estimate = abs(val1 - val2)

        # mark good intervals, gather values and error estimates
        is_good = _numpy_all_except(error_estimate < eps * areas / total_area, axis=-1)
        # add values from good intervals to sum
        quad_sum += sumfun(val1[..., is_good], axis=-1)
        global_error_estimate += sumfun(error_estimate[..., is_good], axis=-1)
        is_bad = np.logical_not(is_good)

    return quad_sum, global_error_estimate

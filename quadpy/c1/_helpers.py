import numpy as np

from ..helpers import QuadratureScheme, plot_disks_1d


def _find_shapes(fx, intervals, x, domain_shape=None, range_shape=None):
    # This following logic is based on the below requirements
    #
    # assert intervals.shape == (2,) + domain_shape + interval_set_shape
    # assert fx.shape == range_shape + interval_set_shape + x.shape
    #
    assert len(x.shape) == 1

    if intervals.shape[0] != 2:
        raise ValueError(
            f"Expected intervals to bbe of shape (2, ...), "
            f"but got shape {intervals.shape} instead."
        )

    if fx.shape[-1] != x.shape[0]:
        raise ValueError(
            f"Expected the function return value to be of shape (..., {x.shape[0]}), "
            f"but got shape {fx.shape} instead."
        )

    if domain_shape is not None and range_shape is not None:
        # Only needed for some assertions
        interval_set_shape = intervals.shape[1 + len(domain_shape) :]
    elif domain_shape is not None:
        interval_set_shape = intervals.shape[1 + len(domain_shape) :]
        range_shape = fx.shape[: -len(interval_set_shape) - 1]
    elif range_shape is not None:
        interval_set_shape = fx.shape[len(range_shape) : -1]
        if len(interval_set_shape) == 0:
            domain_shape = intervals.shape[1:]
        else:
            domain_shape = intervals.shape[1 : -len(interval_set_shape)]
    else:
        # Find the common tail of fx.shape[:-1] and intervals.shape. That is the
        # interval_set_shape unless the tails of domain_shape and range_shape coincide
        # to some degree. (For this case, one can specify domain_shape, range_shape
        # explicitly.)
        interval_set_shape = []
        for k in range(min(len(intervals.shape) - 1, len(fx.shape) - 1)):
            d = intervals.shape[-k - 1]
            if fx.shape[-k - 2] != d:
                break
            interval_set_shape.append(d)
        interval_set_shape = tuple(reversed(interval_set_shape))

        if interval_set_shape == ():
            domain_shape = intervals.shape[1:]
        else:
            domain_shape = intervals.shape[1 : -len(interval_set_shape)]
        range_shape = fx.shape[: -len(interval_set_shape) - 1]

    expected_shape_ivals = (2,) + domain_shape + interval_set_shape
    if intervals.shape != expected_shape_ivals:
        raise ValueError(
            f"Expected intervals to be of shape {expected_shape_ivals}, "
            f"but got shape {intervals.shape} instead."
        )
    expected_shape_fx = range_shape + interval_set_shape + x.shape
    if fx.shape != expected_shape_fx:
        raise ValueError(
            f"Expected the function return value to be of shape {expected_shape_fx}, "
            f"but got shape {fx.shape} instead."
        )
    return domain_shape, range_shape, interval_set_shape


class C1Scheme(QuadratureScheme):
    def __init__(self, name, degree, weights, points, source=None):
        self.domain = "C1"
        super().__init__(name, weights, points, degree, source)

    def integrate(self, f, intervals, domain_shape=None, range_shape=None, dot=np.dot):
        iv = np.asarray(intervals)
        x0 = 0.5 * (1.0 - self.points)
        x1 = 0.5 * (1.0 + self.points)
        x = np.multiply.outer(iv[0], x0) + np.multiply.outer(iv[1], x1)
        fx = np.asarray(f(x))

        domain_shape, range_shape, _ = _find_shapes(
            fx, iv, self.points, domain_shape=domain_shape, range_shape=range_shape
        )

        # np.sum is slower than dot() and friends, but allows for scalar input.
        diff = iv[1] - iv[0]
        len_intervals = np.sqrt(np.sum(diff ** 2, axis=tuple(range(len(domain_shape)))))
        # The factor 0.5 is from the length of the reference line [-1, 1].
        return 0.5 * len_intervals * dot(fx, self.weights)

    def plot(self, interval=np.array([[-1.0], [1.0]]), show_axes=False):
        from matplotlib import pyplot as plt

        # change default range so that new disks will work
        plt.axis("equal")
        # ax.set_xlim((-1.5, 1.5))
        # ax.set_ylim((-1.5, 1.5))

        if not show_axes:
            plt.gca().set_axis_off()

        plt.plot(interval, [0, 0], color="k")

        pts = np.column_stack([self.points, np.zeros(len(self.points))])
        total_area = interval[1] - interval[0]
        plot_disks_1d(plt, pts, self.weights, total_area)

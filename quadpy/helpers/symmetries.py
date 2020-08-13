import numpy


def _zero(data):
    return numpy.array([[0.0], [0.0]])


def _d4(data):
    """dihedral symmetry d4.
    """
    s, t = numpy.array(data)
    points = numpy.array(
        [[+s, +t], [-s, +t], [+s, -t], [-s, -t], [+t, +s], [-t, +s], [+t, -s], [-t, -s]]
    )
    points = numpy.moveaxis(points, 0, 1)
    return points


def _symm_s(a):
    points = numpy.array([[+a, +a], [+a, -a], [-a, +a], [-a, -a]])
    points = numpy.moveaxis(points, 0, 1)
    return points


def _pma(data):
    a = numpy.asarray(data)
    points = numpy.array([[+a, +a], [-a, +a], [+a, -a], [-a, -a]])
    points = numpy.moveaxis(points, 0, 1)
    return points


def _symm_r0(r):
    zero = numpy.zeros_like(r)
    points = numpy.array([[+r, zero], [-r, zero], [zero, +r], [zero, -r]])
    points = numpy.moveaxis(points, 0, 1)
    return points


def _s40(data):
    a = numpy.asarray(data)
    zero = numpy.zeros_like(a)
    points = numpy.array([[+a, zero], [-a, zero], [zero, +a], [zero, -a]])
    points = numpy.moveaxis(points, 0, 1)
    return points


def _s4(data):
    a, b = data
    points = numpy.array([[+a, +b], [-a, -b], [-b, +a], [+b, -a]])
    points = numpy.moveaxis(points, 0, 1)
    return points


def _pm2(data):
    x, y = data
    points = numpy.array([[+x, +y], [+x, -y], [-x, +y], [-x, -y]])
    points = numpy.moveaxis(points, 0, 1)
    return points


def _pm(data):
    x, y = data
    points = numpy.array([[+x, +y], [-x, -y]])
    points = numpy.moveaxis(points, 0, 1)
    return points


def _pmx2(data):
    x, y = data
    points = numpy.array([[+x, y], [-x, y]])
    points = numpy.moveaxis(points, 0, 1)
    return points


def _pmx(r):
    zero = numpy.zeros_like(r)
    points = numpy.array([[+r, zero], [-r, zero]])
    points = numpy.moveaxis(points, 0, 1)
    return points


def _pmy(data):
    a = numpy.asarray(data)
    zero = numpy.zeros_like(a)
    points = numpy.array([[zero, +a], [zero, -a]])
    points = numpy.moveaxis(points, 0, 1)
    return points


def expand_symmetries_points_only(data):
    points = []
    counts = []

    for key, points_raw in data.items():
        fun = {
            "d4": _d4,
            "symm_s": _symm_s,
            "symm_r0": _symm_r0,
            "s4": _s4,
            "zero": _zero,
            "pm2": _pm2,
            "pm": _pm,
            "pmx": _pmx,
            "pmx2": _pmx2,
            "pmy": _pmy,
            "pma": _pma,
            "s40": _s40,
            "plain": lambda vals: vals.reshape(2, 1, -1),
        }[key]
        pts = fun(numpy.asarray(points_raw))

        counts.append(pts.shape[1])
        pts = pts.reshape(pts.shape[0], -1)
        points.append(pts)

    points = numpy.ascontiguousarray(numpy.concatenate(points, axis=1))
    return points, counts


def expand_symmetries(data):
    # separate points and weights
    points_raw = {}
    weights_raw = []
    for key, values in data.items():
        weights_raw.append(values[0])
        points_raw[key] = values[1:]

    points, counts = expand_symmetries_points_only(points_raw)
    weights = numpy.concatenate(
        [numpy.tile(values, count) for count, values in zip(counts, weights_raw)]
    )
    return points, weights

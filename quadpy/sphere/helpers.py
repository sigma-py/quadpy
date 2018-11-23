# -*- coding: utf-8 -*-
#
import numpy
import sympy


def cartesian_to_spherical(X):
    return numpy.stack([numpy.arctan2(X[:, 1], X[:, 0]), numpy.arccos(X[:, 2])], axis=1)


def _atan2_0(X):
    """Like sympy.atan2, but return 0 for x=y=0. Mathematically, the value is
    undefined, so sympy returns NaN, but for the sake of the coordinate
    conversion, its value doesn't matter. NaNs, however, produce NaNs down the
    line.
    """
    out = numpy.array([sympy.atan2(X[k, 1], X[k, 0]) for k in range(len(X))])
    out[out == sympy.nan] = 0
    return out


def cartesian_to_spherical_sympy(X):
    vacos = numpy.vectorize(sympy.acos)
    return numpy.stack([_atan2_0(X), vacos(X[:, 2])], axis=1)


def _a1():
    return numpy.array(
        [
            [+1.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, +1.0, 0.0],
            [0.0, -1.0, 0.0],
            [0.0, 0.0, +1.0],
            [0.0, 0.0, -1.0],
        ]
    )


def _a2():
    return numpy.array(
        [
            [+1.0, +1.0, 0.0],
            [+1.0, -1.0, 0.0],
            [-1.0, +1.0, 0.0],
            [-1.0, -1.0, 0.0],
            #
            [+1.0, 0.0, +1.0],
            [+1.0, 0.0, -1.0],
            [-1.0, 0.0, +1.0],
            [-1.0, 0.0, -1.0],
            #
            [0.0, +1.0, +1.0],
            [0.0, +1.0, -1.0],
            [0.0, -1.0, +1.0],
            [0.0, -1.0, -1.0],
        ]
    ) / numpy.sqrt(2.0)


def _a3():
    return numpy.array(
        [
            [+1.0, +1.0, +1.0],
            [+1.0, +1.0, -1.0],
            [+1.0, -1.0, +1.0],
            [+1.0, -1.0, -1.0],
            [-1.0, +1.0, +1.0],
            [-1.0, +1.0, -1.0],
            [-1.0, -1.0, +1.0],
            [-1.0, -1.0, -1.0],
        ]
    ) / numpy.sqrt(3.0)


def _pq0(alpha):
    a = numpy.sin(alpha * numpy.pi)
    b = numpy.cos(alpha * numpy.pi)
    zero = numpy.zeros_like(alpha)
    return numpy.array(
        [
            [+a, +b, zero],
            [-a, +b, zero],
            [-a, -b, zero],
            [+a, -b, zero],
            #
            [+b, +a, zero],
            [-b, +a, zero],
            [-b, -a, zero],
            [+b, -a, zero],
            #
            [+a, zero, +b],
            [-a, zero, +b],
            [-a, zero, -b],
            [+a, zero, -b],
            #
            [+b, zero, +a],
            [-b, zero, +a],
            [-b, zero, -a],
            [+b, zero, -a],
            #
            [zero, +a, +b],
            [zero, -a, +b],
            [zero, -a, -b],
            [zero, +a, -b],
            #
            [zero, +b, +a],
            [zero, -b, +a],
            [zero, -b, -a],
            [zero, +b, -a],
        ]
    )


def _pq02(a):
    b = numpy.sqrt(1 - a ** 2)
    zero = numpy.zeros_like(a)
    return numpy.array(
        [
            [+a, +b, zero],
            [-a, +b, zero],
            [-a, -b, zero],
            [+a, -b, zero],
            #
            [+b, +a, zero],
            [-b, +a, zero],
            [-b, -a, zero],
            [+b, -a, zero],
            #
            [+a, zero, +b],
            [-a, zero, +b],
            [-a, zero, -b],
            [+a, zero, -b],
            #
            [+b, zero, +a],
            [-b, zero, +a],
            [-b, zero, -a],
            [+b, zero, -a],
            #
            [zero, +a, +b],
            [zero, -a, +b],
            [zero, -a, -b],
            [zero, +a, -b],
            #
            [zero, +b, +a],
            [zero, -b, +a],
            [zero, -b, -a],
            [zero, +b, -a],
        ]
    )


def _llm(beta):
    # translate the point into cartesian coords; note that phi=pi/4.
    beta *= numpy.pi
    L = numpy.sin(beta) / numpy.sqrt(2)
    m = numpy.cos(beta)
    return numpy.array(
        [
            [+L, +L, +m],
            [-L, +L, +m],
            [+L, -L, +m],
            [-L, -L, +m],
            [+L, +L, -m],
            [-L, +L, -m],
            [+L, -L, -m],
            [-L, -L, -m],
            #
            [+L, +m, +L],
            [-L, +m, +L],
            [+L, +m, -L],
            [-L, +m, -L],
            [+L, -m, +L],
            [-L, -m, +L],
            [+L, -m, -L],
            [-L, -m, -L],
            #
            [+m, +L, +L],
            [+m, -L, +L],
            [+m, +L, -L],
            [+m, -L, -L],
            [-m, +L, +L],
            [-m, -L, +L],
            [-m, +L, -L],
            [-m, -L, -L],
        ]
    )


def _llm2(L):
    m = numpy.sqrt(1 - 2 * L ** 2)
    return numpy.array(
        [
            [+L, +L, +m],
            [-L, +L, +m],
            [+L, -L, +m],
            [-L, -L, +m],
            [+L, +L, -m],
            [-L, +L, -m],
            [+L, -L, -m],
            [-L, -L, -m],
            #
            [+L, +m, +L],
            [-L, +m, +L],
            [+L, +m, -L],
            [-L, +m, -L],
            [+L, -m, +L],
            [-L, -m, +L],
            [+L, -m, -L],
            [-L, -m, -L],
            #
            [+m, +L, +L],
            [+m, -L, +L],
            [+m, +L, -L],
            [+m, -L, -L],
            [-m, +L, +L],
            [-m, -L, +L],
            [-m, +L, -L],
            [-m, -L, -L],
        ]
    )


def _rsw(azimuthal, polar):
    # translate the point into cartesian coords; note that phi=pi/4.
    azimuthal *= numpy.pi
    polar *= numpy.pi

    sin_polar = numpy.sin(polar)
    cos_polar = numpy.cos(polar)
    sin_azimuthal = numpy.sin(azimuthal)
    cos_azimuthal = numpy.cos(azimuthal)

    r = sin_polar * cos_azimuthal
    s = sin_polar * sin_azimuthal
    w = cos_polar

    return numpy.array(
        [
            [+r, +s, +w],
            [+w, +r, +s],
            [+s, +w, +r],
            [+s, +r, +w],
            [+w, +s, +r],
            [+r, +w, +s],
            #
            [-r, +s, +w],
            [+w, -r, +s],
            [+s, +w, -r],
            [+s, -r, +w],
            [+w, +s, -r],
            [-r, +w, +s],
            #
            [+r, -s, +w],
            [+w, +r, -s],
            [-s, +w, +r],
            [-s, +r, +w],
            [+w, -s, +r],
            [+r, +w, -s],
            #
            [+r, +s, -w],
            [-w, +r, +s],
            [+s, -w, +r],
            [+s, +r, -w],
            [-w, +s, +r],
            [+r, -w, +s],
            #
            [-r, -s, +w],
            [+w, -r, -s],
            [-s, +w, -r],
            [-s, -r, +w],
            [+w, -s, -r],
            [-r, +w, -s],
            #
            [-r, +s, -w],
            [-w, -r, +s],
            [+s, -w, -r],
            [+s, -r, -w],
            [-w, +s, -r],
            [-r, -w, +s],
            #
            [+r, -s, -w],
            [-w, +r, -s],
            [-s, -w, +r],
            [-s, +r, -w],
            [-w, -s, +r],
            [+r, -w, -s],
            #
            [-r, -s, -w],
            [-w, -r, -s],
            [-s, -w, -r],
            [-s, -r, -w],
            [-w, -s, -r],
            [-r, -w, -s],
        ]
    )


def untangle2(data):
    points = []
    weights = []
    if "a1" in data:
        assert len(data["a1"]) == 1
        points.append(_a1())
        w = data["a1"][0]
        weights.append(numpy.full(6, w))

    if "a2" in data:
        assert len(data["a2"]) == 1
        points.append(_a2())
        w = data["a2"][0]
        weights.append(numpy.full(12, w))

    if "a3" in data:
        assert len(data["a3"]) == 1
        points.append(_a3())
        w = data["a3"][0]
        weights.append(numpy.full(8, w))

    if "llm" in data:
        beta = numpy.array(data["llm"])[:, 1]
        out = _collapse0(numpy.moveaxis(_llm(beta), 0, 1)).T
        points.append(out)
        w = numpy.array(data["llm"])[:, 0]
        weights.append(numpy.tile(w, 24))

    if "llm2" in data:
        beta = numpy.array(data["llm2"])[:, 1]
        out = _collapse0(numpy.moveaxis(_llm2(beta), 0, 1)).T
        points.append(out)
        w = numpy.array(data["llm2"])[:, 0]
        weights.append(numpy.tile(w, 24))

    if "pq0" in data:
        beta = numpy.array(data["pq0"])[:, 1]
        out = _collapse0(numpy.moveaxis(_pq0(beta), 0, 1)).T
        points.append(out)
        w = numpy.array(data["pq0"])[:, 0]
        weights.append(numpy.tile(w, 24))

    if "pq02" in data:
        beta = numpy.array(data["pq02"])[:, 1]
        out = _collapse0(numpy.moveaxis(_pq02(beta), 0, 1)).T
        points.append(out)
        w = numpy.array(data["pq02"])[:, 0]
        weights.append(numpy.tile(w, 24))

    if "rsw" in data:
        beta = numpy.array(data["rsw"])[:, 1:].T
        out = _collapse0(numpy.moveaxis(_rsw(*beta), 0, 1)).T
        points.append(out)
        w = numpy.array(data["rsw"])[:, 0]
        weights.append(numpy.tile(w, 48))

    points = numpy.concatenate(points)
    weights = numpy.concatenate(weights)
    return points, weights


def _collapse0(a):
    """Collapse all dimensions of `a` except the first.
    """
    return a.reshape(a.shape[0], -1)

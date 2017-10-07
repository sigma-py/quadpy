# -*- coding: utf-8 -*-
#
import matplotlib.pyplot as plt
import numpy

from .. import helpers


def area(radius):
    return 4*numpy.pi * numpy.array(radius)**2


def show(*args, **kwargs):
    plot(*args, **kwargs)
    plt.show()
    return


def plot(scheme, backend='mpl'):
    backend_to_function = {
        'mpl': _plot_mpl
        }
    backend_to_function[backend](scheme)
    return


def _plot_mpl(scheme):
    # pylint: disable=relative-import, unused-variable
    from mpl_toolkits.mplot3d import Axes3D

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.set_aspect('equal')

    flt = numpy.vectorize(float)
    pts = flt(scheme.points)
    wgs = flt(scheme.weights)

    for p, w in zip(pts, wgs):
        # <https://en.wikipedia.org/wiki/Spherical_cap>
        w *= 4 * numpy.pi
        theta = numpy.arccos(1.0 - abs(w) / (2*numpy.pi))
        color = '#1f77b4' if w >= 0 else '#d62728'
        _plot_spherical_cap_mpl(ax, p, theta, color)

    ax.set_axis_off()
    return


def _plot_spherical_cap_mpl(ax, b, opening_angle, color, elevation=1.01):
    r = elevation
    phi = numpy.linspace(0, 2 * numpy.pi, 30)
    theta = numpy.linspace(0, opening_angle, 20)
    X = r * numpy.stack([
        numpy.outer(numpy.cos(phi), numpy.sin(theta)),
        numpy.outer(numpy.sin(phi), numpy.sin(theta)),
        numpy.outer(numpy.ones(numpy.size(phi)), numpy.cos(theta)),
        ], axis=-1)

    # rotate X such that [0, 0, 1] gets rotated to `c`;
    # <https://math.stackexchange.com/a/476311/36678>.
    a = numpy.array([0.0, 0.0, 1.0])
    a_x_b = numpy.cross(a, b)
    a_dot_b = numpy.dot(a, b)
    if a_dot_b == -1.0:
        X_rot = -X
    else:
        X_rot = (
            X +
            numpy.cross(a_x_b, X) +
            numpy.cross(a_x_b, numpy.cross(a_x_b, X)) / (1.0 + a_dot_b)
            )

    ax.plot_surface(
            X_rot[..., 0], X_rot[..., 1], X_rot[..., 2],
            rstride=3, cstride=3,
            color=color,
            alpha=0.5,
            linewidth=0
            )
    return


def integrate(f, center, radius, rule, sumfun=helpers.kahan_sum):
    '''Quadrature where `f` is defined in Cartesian coordinates.
    '''
    center = numpy.array(center)
    rr = numpy.multiply.outer(radius, rule.points)
    rr = numpy.swapaxes(rr, 0, -2)
    ff = numpy.array(f((rr + center).T))
    return area(radius) * sumfun(rule.weights * ff, axis=-1)


def integrate_spherical(f, radius, rule, sumfun=helpers.kahan_sum):
    '''Quadrature where `f` is defined in spherical coordinates `phi`, `theta`.

    This is using the ISO 31-11 convection of `theta` being the polar, `phi`
    the azimuthal angle. Note that this meaning is often reversed in physics,
    cf. <http://mathworld.wolfram.com/SphericalCoordinates.html>.
    '''
    rr = numpy.swapaxes(rule.phi_theta, 0, -2)
    ff = numpy.array(f(rr.T))
    return area(radius) * sumfun(rule.weights * ff, axis=-1)

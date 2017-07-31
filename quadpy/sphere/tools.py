# -*- coding: utf-8 -*-
#
import numpy

from .. import helpers


def area(radius):
    return 4*numpy.pi * numpy.array(radius)**2


def show(scheme):
    from matplotlib import pyplot as plt
    # pylint: disable=relative-import, unused-variable
    from mpl_toolkits.mplot3d import Axes3D

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.set_aspect('equal')

    # http://matplotlib.org/examples/mplot3d/surface3d_demo2.html
    u = numpy.linspace(0, 2 * numpy.pi, 100)
    v = numpy.linspace(0, numpy.pi, 100)
    x = numpy.outer(numpy.cos(u), numpy.sin(v))
    y = numpy.outer(numpy.sin(u), numpy.sin(v))
    z = numpy.outer(numpy.ones(numpy.size(u)), numpy.cos(v))

    ax.plot_surface(
            x, y, z,
            rstride=3, cstride=3,
            color='0.9',
            alpha=1.0,
            linewidth=0
            )

    ax.scatter(
        1.05 * scheme.points[:, 0],
        1.05 * scheme.points[:, 1],
        1.05 * scheme.points[:, 2],
        color='#1f77b4',
        s=60
        )

    ax.set_axis_off()
    plt.show()
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

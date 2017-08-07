# -*- coding: utf-8 -*-
#
from __future__ import division

import itertools
import math
from math import sqrt, factorial as fact
import numpy


def untangle(data):
    weights, points = zip(*data)
    return (
        numpy.concatenate(points),
        numpy.repeat(weights, [len(grp) for grp in points])
        )


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


def show_mpl(points, weights, volume, edges):
    import matplotlib.pyplot as plt
    # pylint: disable=relative-import, unused-variable
    from mpl_toolkits.mplot3d import Axes3D

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
            # different size with different number of points. Another options
            # would be to precompute x, y, z before the loop, but this can be
            # heavy on the graphics output. See
            # <https://stackoverflow.com/q/45324258/353337>.
            u = numpy.linspace(0, 2 * numpy.pi, int(2*numpy.pi/h*r) + 1)
            v = numpy.linspace(0, numpy.pi, int(numpy.pi/h*r) + 1)
            _x = numpy.outer(numpy.cos(u), numpy.sin(v))
            _y = numpy.outer(numpy.sin(u), numpy.sin(v))
            _z = numpy.outer(numpy.ones(numpy.size(u)), numpy.cos(v))

            color = '#1f77b4' if weight >= 0 else '#d62728'
            # highlight ball center
            plt.plot(
                [tp[0]], [tp[1]], [tp[2]],
                linestyle='None', marker='.', color=color
                )

            ax.plot_surface(
                r*_x + tp[0], r*_y + tp[1], r*_z + tp[2],
                color=color,
                alpha=0.3,
                linewidth=1
                )

        ax.set_axis_off()
        return

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.set_aspect('equal')
    ax.set_axis_off()

    for edge in edges:
        plt.plot(*edge, color='k', linestyle='-')

    plot_spheres(plt, ax, points, weights, volume)
    plt.show()
    return


# pylint: disable=too-many-locals
def show_mayavi(points, weights, volume, edges):
    # pylint: disable=import-error
    import mayavi.mlab as mlab

    mlab.figure(bgcolor=(1.0, 1.0, 1.0))

    for edge in edges:
        mlab.plot3d(*edge, tube_radius=0.5e-2, color=(0.0, 0.0, 0.0))

    blue = (31./255., 119.0/255., 180./255.)
    red = (84./255., 15.0/255., 16./255.)

    h = 1.0e-2
    sum_weights = math.fsum(weights)
    for tp, weight in zip(points, weights):
        # Choose radius such that the sum of volumes of the balls equals
        # total_volume.
        r = (
            abs(weight)/sum_weights * volume/(4.0/3.0 * numpy.pi)
            )**(1.0/3.0)

        # Create a sphere
        u = numpy.linspace(0, 2 * numpy.pi, int(2*numpy.pi/h*r) + 1)
        v = numpy.linspace(0, numpy.pi, int(numpy.pi/h*r) + 1)
        sin_u, cos_u = numpy.sin(u), numpy.cos(u)
        sin_v, cos_v = numpy.sin(v), numpy.cos(v)
        _x = numpy.outer(cos_u, sin_v)
        _y = numpy.outer(sin_u, sin_v)
        _z = numpy.outer(numpy.ones(numpy.size(u)), cos_v)

        mlab.mesh(
            r*_x + tp[0], r*_y + tp[1], r*_z + tp[2],
            color=blue if weight >= 0 else red,
            opacity=1.0
            )
    mlab.show()
    return


# pylint: disable=too-many-locals
def show_vtk(points, weights, volume, edges):
    # pylint: disable=import-error
    import vtk

    def get_line_actor(x0, x1):
        source = vtk.vtkLineSource()
        source.SetPoint1(x0)
        source.SetPoint2(x1)
        # mapper
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(source.GetOutputPort())
        # actor
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        # color actor
        actor.GetProperty().SetColor(0, 0, 0)
        return actor

    def get_sphere_actor(x0, r, color):
        # Generate polygon data for a sphere
        sphere = vtk.vtkSphereSource()

        sphere.SetCenter(x0)
        sphere.SetRadius(r)

        sphere.SetPhiResolution(100)
        sphere.SetThetaResolution(100)

        # Create a mapper for the sphere data
        sphere_mapper = vtk.vtkPolyDataMapper()
        # sphere_mapper.SetInput(sphere.GetOutput())
        sphere_mapper.SetInputConnection(sphere.GetOutputPort())

        # Connect the mapper to an actor
        sphere_actor = vtk.vtkActor()
        sphere_actor.SetMapper(sphere_mapper)
        sphere_actor.GetProperty().SetColor(color)
        # sphere_actor.GetProperty().SetOpacity(0.9)
        return sphere_actor

    line_actors = [get_line_actor(edge[:, 0], edge[:, 1]) for edge in edges]

    blue = numpy.array([31.0, 119.0, 180.0]) / 255.0
    red = numpy.array([84.0, 15.0, 16.0]) / 255.0

    sum_weights = math.fsum(weights)
    sphere_actors = [
        get_sphere_actor(
            pt,
            numpy.cbrt(abs(weight)/sum_weights * volume/(4.0/3.0 * numpy.pi)),
            color=blue if weight > 0.0 else red
            )
        for pt, weight in zip(points, weights)
        ]

    # Create a renderer and add the sphere actor to it
    renderer = vtk.vtkRenderer()
    renderer.SetBackground(1.0, 1.0, 1.0)
    for sphere_actor in sphere_actors:
        renderer.AddActor(sphere_actor)
    for line_actor in line_actors:
        renderer.AddActor(line_actor)

    # Create a render window
    render_window = vtk.vtkRenderWindow()
    render_window.AddRenderer(renderer)

    # Create an interactor
    interactor = vtk.vtkRenderWindowInteractor()
    interactor.SetRenderWindow(render_window)
    # Initialize the interactor and start the rendering loop
    interactor.Initialize()
    render_window.Render()
    interactor.Start()
    return


backend_to_function = {
    'mayavi': show_mayavi,
    'mpl': show_mpl,
    'vtk': show_vtk,
    }


def z(n):
    return numpy.zeros((1, n))


def fsd(n, r, d):
    '''Get all permutations of [+-r, +-r, 0, ..., 0] of length n, where +-r
    occurs d times.
    len(out) == 2**d * (n over d).
    n==1:  2*n
    n==2:  2*n*(n-1)
    n==3:  4*n*(n-1)*(n-2) / 3
    '''
    assert 0 <= d <= n
    return combine([[+r, -r]] * d + [[0.0]] * (n-d))


def fsd2(n, r, s, i, j):
    '''Get all permutations of [+-r, +-r, +-s, +-s, 0, ..., 0] of length n,
    with i times the number r and and j times the number s.
    '''
    assert 2 <= i+j <= n
    return combine([[+r, -r]] * i + [[+s, -s]] * j + [[0.0]] * (n-i-j))


def combine(pools):
    '''Given an input array with lists of options, e.g.,

    [[a, b], [c], [d]],

    this methods returns all combinations with one element from each
    subset, e.g.,

    [a, c, d], [a, d, c], [c, d, a], ...
    [b, c, d], [b, d, c], [c, d, b], ...
    '''
    # https://stackoverflow.com/a/45322199/353337
    return numpy.array(list(set(itertools.chain.from_iterable([
        itertools.permutations(x) for x in itertools.product(*pools)
        ]))))


def pm(n, a):
    '''Return all combinations of [+a, -a] with length n (with repetition).
    len(out) == 2**n.
    '''
    return numpy.array(list(itertools.product([+a, -a], repeat=n)))


# pylint: disable=too-many-arguments
def compute_dobrodeev(n, I0, I2, I22, I4, pm_type, i, j, k):
    '''Compute some helper quantities used in

    L.N. Dobrodeev,
    Cubature rules with equal coefficients for integrating functions with
    respect to symmetric domains,
    USSR Computational Mathematics and Mathematical Physics,
    Volume 18, Issue 4, 1978, Pages 27-34,
    <https://doi.org/10.1016/0041-5553(78)90064-2>.
    '''
    t = 1 if pm_type == 'I' else -1

    L = fact(n) // (fact(i) * fact(n-i)) * 2**i
    M = fact(n) // (fact(j) * fact(k) * fact(n-j-k)) * 2**(j+k)
    N = L + M
    F = I22/I0 - I2**2/I0**2 + (I4/I0 - I22/I0) / n
    R = -(j+k-i) / i * I2**2/I0**2 + (j+k-1)/n * I4/I0 - (n-1)/n * I22/I0
    H = 1/i * (
        (j+k-i) * I2**2/I0**2 + (j+k)/n * ((i-1) * I4/I0 - (n-1)*I22/I0)
        )
    Q = L/M*R + H - t * 2*I2/I0 * (j+k-i)/i * sqrt(L/M*F)

    G = 1/N
    a = sqrt(n/i * (I2/I0 + t * sqrt(M/L*F)))
    b = sqrt(n/(j+k) * (I2/I0 - t * sqrt(L/M*F) + t * sqrt(k/j*Q)))
    c = sqrt(n/(j+k) * (I2/I0 - t * sqrt(L/M*F) - t * sqrt(j/k*Q)))
    return G, a, b, c

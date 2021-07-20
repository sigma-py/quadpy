import math

import numpy as np

__all__ = ["plot_disks_1d", "plot_disks"]


def plot_disks_1d(plt, pts, weights, total_area):
    """Plot a circles at quadrature points according to weights. The diameters
    sum up to the total area.
    """
    radii = 0.5 * abs(weights) / math.fsum(weights) * total_area
    colors = ["tab:blue" if weight >= 0 else "tab:red" for weight in weights]
    _plot_disks_helpers(plt, pts, radii, colors)


def plot_disks(plt, pts, weights, total_area):
    """Plot a circles at quadrature points according to weights."""
    flt = np.vectorize(float)
    pts = flt(pts)
    weights = flt(weights)
    radii = np.sqrt(abs(weights) / math.fsum(weights) * total_area / math.pi)
    colors = [
        # use matplotlib 2.0's color scheme
        "tab:blue" if weight >= 0 else "tab:red"
        for weight in weights
    ]
    _plot_disks_helpers(plt, pts, radii, colors)


def _plot_disks_helpers(plt, pts, radii, colors):
    for pt, radius, color in zip(pts, radii, colors):
        # highlight circle center
        plt.plot([pt[0]], [pt[1]], linestyle="None", marker=".", color=color)
        # Choose radius such that the sum of areas of the circles equals total_area.
        # Make sure to set the line width to 0,
        # <https://github.com/matplotlib/matplotlib/issues/17421>.
        circ = plt.Circle((pt[0], pt[1]), radius, color=color, alpha=0.5, linewidth=0)
        plt.gca().add_patch(circ)


def show_mpl(points, weights, volume, edges, balls=None):
    from matplotlib import pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    flt = np.vectorize(float)
    points = flt(points)
    weights = flt(weights)

    def plot_spheres(plt, ax, pts, radii, colors):
        h = 1.0e-2

        for tp, r, color in zip(pts, radii, colors):
            # https://matplotlib.org/examples/mplot3d/surface3d_demo2.html
            # Compute sphere for every point anew. This is more costly on the
            # numerical side, but gives the flexibility of drawing sphere of
            # different size with different number of points. Another options
            # would be to precompute x, y, z before the loop, but this can be
            # heavy on the graphics output. See
            # <https://stackoverflow.com/q/45324258/353337>.
            u = np.linspace(0, 2 * np.pi, int(2 * np.pi / h * r) + 1)
            v = np.linspace(0, np.pi, int(np.pi / h * r) + 1)
            _x = np.outer(np.cos(u), np.sin(v))
            _y = np.outer(np.sin(u), np.sin(v))
            _z = np.outer(np.ones(np.size(u)), np.cos(v))

            # highlight ball center
            plt.plot(
                [tp[0]], [tp[1]], [tp[2]], linestyle="None", marker=".", color=color
            )

            ax.plot_surface(
                r * _x + tp[0],
                r * _y + tp[1],
                r * _z + tp[2],
                color=color,
                alpha=0.3,
                linewidth=1,
            )

        ax.set_axis_off()
        return

    balls = [] if balls is None else balls

    ax = plt.axes(projection=Axes3D.name)
    # ax.set_aspect("equal")
    ax.set_axis_off()

    for edge in edges:
        plt.plot(*edge, color="k", linestyle="-")

    plot_spheres(
        plt,
        ax,
        points,
        # Choose radius such that the sum of volumes of the balls equals
        # total_volume.
        radii=np.cbrt(abs(weights) / math.fsum(weights) * volume / (4.0 / 3.0 * np.pi)),
        colors=["tab:blue" if weight >= 0 else "tab:red" for weight in weights],
    )

    for ball in balls:
        plot_spheres(plt, ax, [ball[0]], [ball[1]], ["#dddddd"])

    plt.show()
    return plt


# def show_mayavi(points, weights, volume, edges, balls=None):
#     import mayavi.mlab as mlab
#
#     mlab.figure(bgcolor=(1.0, 1.0, 1.0))
#
#     for edge in edges:
#         mlab.plot3d(*edge, tube_radius=0.5e-2, color=(0.0, 0.0, 0.0))
#
#     blue = (31.0 / 255.0, 119.0 / 255.0, 180.0 / 255.0)
#     red = (84.0 / 255.0, 15.0 / 255.0, 16.0 / 255.0)
#
#     h = 1.0e-2
#     sum_weights = math.fsum(weights)
#     for tp, weight in zip(points, weights):
#         # Choose radius such that the sum of volumes of the balls equals
#         # total_volume.
#         r = (abs(weight) / sum_weights * volume / (4.0 / 3.0 * np.pi)) ** (1.0 / 3.0)
#
#         # Create a sphere
#         u = np.linspace(0, 2 * np.pi, int(2 * np.pi / h * r) + 1)
#         v = np.linspace(0, np.pi, int(np.pi / h * r) + 1)
#         sin_u, cos_u = np.sin(u), np.cos(u)
#         sin_v, cos_v = np.sin(v), np.cos(v)
#         _x = np.outer(cos_u, sin_v)
#         _y = np.outer(sin_u, sin_v)
#         _z = np.outer(np.ones(np.size(u)), cos_v)
#
#         mlab.mesh(
#             r * _x + tp[0],
#             r * _y + tp[1],
#             r * _z + tp[2],
#             color=blue if weight >= 0 else red,
#             opacity=1.0,
#         )
#
#     balls = [] if balls is None else balls
#     for ball in balls:
#         tp = ball[0]
#         r = ball[1]
#
#         # Create a sphere
#         u = np.linspace(0, 2 * np.pi, int(2 * np.pi / h * r) + 1)
#         v = np.linspace(0, np.pi, int(np.pi / h * r) + 1)
#         sin_u, cos_u = np.sin(u), np.cos(u)
#         sin_v, cos_v = np.sin(v), np.cos(v)
#         _x = np.outer(cos_u, sin_v)
#         _y = np.outer(sin_u, sin_v)
#         _z = np.outer(np.ones(np.size(u)), cos_v)
#
#         mlab.mesh(
#             r * _x + tp[0], r * _y + tp[1], r * _z + tp[2], color=[0, 0, 0], opacity=1.0
#         )
#
#     mlab.show()
#     return


def show_vtk(points, weights, volume, edges, balls=None, render=True):
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

    def get_sphere_actor(x0, r, color, opacity=1.0):
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
        sphere_actor.GetProperty().SetOpacity(opacity)
        return sphere_actor

    flt = np.vectorize(float)
    points = flt(points)
    weights = flt(weights)

    balls = [] if balls is None else balls

    line_actors = [get_line_actor(edge[:, 0], edge[:, 1]) for edge in edges]

    blue = np.array([31.0, 119.0, 180.0]) / 255.0
    red = np.array([84.0, 15.0, 16.0]) / 255.0

    radii = np.cbrt(abs(weights) / math.fsum(weights) * volume / (4.0 / 3.0 * np.pi))
    sphere_actors = [
        get_sphere_actor(pt, radius, color=blue if weight > 0.0 else red)
        for pt, weight, radius in zip(points.T, weights, radii)
    ]

    sphere_actors.extend(
        [
            get_sphere_actor(
                np.array(ball[0]),
                ball[1],
                color=np.array([0.0, 0.0, 0.0]) / 255.0,
                opacity=0.5,
            )
            for ball in balls
        ]
    )

    # Create a renderer and add the sphere actor to it
    renderer = vtk.vtkRenderer()
    renderer.SetBackground(1.0, 1.0, 1.0)
    # Available in more recent versions of VTK
    # <https://vtk.org/doc/nightly/html/classvtkViewport.html#aed4374e05dbbea1692f7c9c865407664>
    # renderer.SetBackgroundAlpha(1.0)
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

    if render:
        # Initialize the interactor and start the rendering loop
        interactor.Initialize()
        render_window.Render()
        interactor.Start()

    # # Screenshot
    # TODO transparent background
    # w2if = vtk.vtkWindowToImageFilter()
    # w2if.SetInput(render_window)
    # w2if.Update()
    # writer = vtk.vtkPNGWriter()
    # writer.SetFileName('screenshot.png')
    # writer.SetInputConnection(w2if.GetOutputPort())
    # writer.Write()
    return


backend_to_function = {
    # "mayavi": show_mayavi,
    "mpl": show_mpl,
    "vtk": show_vtk,
}

import numpy
import quadpy


def test_circle():
    quadpy.circle.integrate(
        lambda x: numpy.exp(x[0]),
        numpy.array([0.0, 0.3]), 0.7,
        quadpy.circle.Equidistant(3)
        )
    quadpy.circle.integrate(
        lambda x: [numpy.exp(x[0]), numpy.exp(x[0])],
        numpy.array([[1.0, 1.0], [0.0, 0.3], [2.0, 2.0]]),
        [1.0, 0.7, 0.333],
        quadpy.circle.Equidistant(5)
        )
    return


def test_disk():
    quadpy.disk.integrate(
        lambda x: numpy.exp(x[0]),
        numpy.array([0.0, 0.3]), 0.7,
        quadpy.disk.Peirce(5)
        )
    quadpy.disk.integrate(
        lambda x: [numpy.exp(x[0]), numpy.exp(x[1])],
        numpy.array([[1.0, 1.0], [0.0, 0.3], [2.0, 2.0]]),
        [1.0, 0.7, 0.333],
        quadpy.disk.Peirce(5)
        )
    return


def test_hexahedron():
    quadpy.hexahedron.integrate(
            lambda x: numpy.exp(x[0]),
            numpy.array([
                [0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
                [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1],
                ]),
            quadpy.hexahedron.Product(quadpy.line_segment.NewtonCotesClosed(3))
            )
    quadpy.hexahedron.integrate(
            lambda x: [numpy.exp(x[0]), numpy.exp(x[1])],
            numpy.array([
                [[0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0]],
                [[1, 0, 0], [1, 0, 0], [1, 0, 0], [1, 0, 0], [1, 0, 0]],
                [[1, 1, 0], [1, 1, 0], [1, 1, 0], [1, 1, 0], [1, 1, 0]],
                [[0, 1, 0], [0, 1, 0], [0, 1, 0], [0, 1, 0], [0, 1, 0]],
                #
                [[0, 0, 1], [0, 0, 1], [0, 0, 1], [0, 0, 1], [0, 0, 1]],
                [[1, 0, 1], [1, 0, 1], [1, 0, 1], [1, 0, 1], [1, 0, 1]],
                [[1, 1, 1], [1, 1, 1], [1, 1, 1], [1, 1, 1], [1, 1, 1]],
                [[0, 1, 1], [0, 1, 1], [0, 1, 1], [0, 1, 1], [0, 1, 1]],
                ]),
            quadpy.hexahedron.Product(quadpy.line_segment.NewtonCotesClosed(3))
            )
    return


def test_line_segment():
    quadpy.line_segment.integrate(
            lambda x: numpy.exp(x),
            numpy.array([0.0, 1.0]),
            quadpy.line_segment.NewtonCotesClosed(4)
            )
    quadpy.line_segment.integrate(
            lambda x: [numpy.exp(x), numpy.sin(x), numpy.cos(x)],
            numpy.array([
                [0.0, 1.0],
                [1.0, 2.0],
                ]),
            quadpy.line_segment.NewtonCotesClosed(4)
            )
    return


def test_pyramid():
    quadpy.pyramid.integrate(
            lambda x: numpy.exp(x[0]),
            numpy.array([
                [0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
                [0, 0, 1],
                ]),
            quadpy.pyramid.Felippa(3)
            )

    quadpy.pyramid.integrate(
            lambda x: [numpy.exp(x[0]), numpy.exp(x[1])],
            numpy.array([
                [[0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0]],
                [[1, 0, 0], [1, 0, 0], [1, 0, 0], [1, 0, 0], [1, 0, 0]],
                [[1, 1, 0], [1, 1, 0], [1, 1, 0], [1, 1, 0], [1, 1, 0]],
                [[0, 1, 0], [0, 1, 0], [0, 1, 0], [0, 1, 0], [0, 1, 0]],
                [[0, 0, 1], [0, 0, 1], [0, 0, 1], [0, 0, 1], [0, 0, 1]],
                ]),
            quadpy.pyramid.Felippa(5)
            )
    return


def test_quadrilateral():
    quadpy.quadrilateral.integrate(
            lambda x: numpy.exp(x[0]),
            numpy.array([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]),
            quadpy.quadrilateral.Stroud(6)
            )

    quadpy.quadrilateral.integrate(
            lambda x: [numpy.exp(x[0]), numpy.exp(x[1])],
            numpy.array([
                [[0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0]],
                [[1.0, 0.0], [1.0, 0.0], [1.0, 0.0], [1.0, 0.0], [1.0, 0.0]],
                [[1.0, 1.0], [1.0, 1.0], [1.0, 1.0], [1.0, 1.0], [1.0, 1.0]],
                [[0.0, 1.0], [0.0, 1.0], [0.0, 1.0], [0.0, 1.0], [0.0, 1.0]],
                ]),
            quadpy.quadrilateral.Stroud(6)
            )
    return


def test_sphere():
    quadpy.sphere.integrate(
            lambda x: numpy.exp(x[0]),
            numpy.array([0.0, 0.3, 0.0]), 0.7,
            quadpy.sphere.Lebedev(3)
            )

    quadpy.sphere.integrate(
            lambda x: [numpy.exp(x[0]), numpy.exp(x[1])],
            numpy.array([[1.0, 1.0, 0.0], [0.0, 0.3, 0.0], [2.0, 2.0, 0.0]]),
            [1.0, 0.7, 0.333],
            quadpy.sphere.Lebedev(3)
            )
    return


def test_tetrahedron():
    quadpy.tetrahedron.integrate(
            lambda x: numpy.exp(x[0]),
            numpy.array([
                [0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]
                ], dtype=float),
            quadpy.tetrahedron.ShunnHam(3)
            )

    quadpy.tetrahedron.integrate(
            lambda x: [numpy.exp(x[0]), numpy.exp(x[1])],
            numpy.array([
                [[0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0]],
                [[1, 0, 0], [1, 0, 0], [1, 0, 0], [1, 0, 0], [1, 0, 0]],
                [[0, 1, 0], [0, 1, 0], [0, 1, 0], [0, 1, 0], [0, 1, 0]],
                [[0, 0, 1], [0, 0, 1], [0, 0, 1], [0, 0, 1], [0, 0, 1]],
                ], dtype=float),
            quadpy.tetrahedron.ShunnHam(3)
            )
    return


def test_triangle():
    quadpy.triangle.integrate(
            lambda x: numpy.exp(x[0]),
            numpy.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]),
            quadpy.triangle.Cubtri()
            )

    quadpy.triangle.integrate(
            lambda x: [numpy.exp(x[0]), numpy.exp(x[1])],
            numpy.array([
                [[0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0]],
                [[1.0, 0.0], [1.0, 0.0], [1.0, 0.0], [1.0, 0.0], [1.0, 0.0]],
                [[0.0, 1.0], [0.0, 1.0], [0.0, 1.0], [0.0, 1.0], [0.0, 1.0]],
                ]),
            quadpy.triangle.Cubtri()
            )
    return


def test_wedge():
    quadpy.wedge.integrate(
            lambda x: numpy.exp(x[0]),
            numpy.array([
                [0, 0, 0], [1, 0, 0], [0, 1, 0],
                [0, 0, 1], [1, 0, 1], [0, 1, 1],
                ]),
            quadpy.wedge.Felippa(4)
            )

    quadpy.wedge.integrate(
            lambda x: [numpy.exp(x[0]), numpy.exp(x[1])],
            numpy.array([
                [[0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0], [0, 0, 0]],
                [[1, 0, 0], [1, 0, 0], [1, 0, 0], [1, 0, 0], [1, 0, 0]],
                [[0, 1, 0], [0, 1, 0], [0, 1, 0], [0, 1, 0], [0, 1, 0]],
                [[0, 0, 1], [0, 0, 1], [0, 0, 1], [0, 0, 1], [0, 0, 1]],
                [[1, 0, 1], [1, 0, 1], [1, 0, 1], [1, 0, 1], [1, 0, 1]],
                [[0, 1, 1], [0, 1, 1], [0, 1, 1], [0, 1, 1], [0, 1, 1]],
                ]),
            quadpy.wedge.Felippa(4)
            )
    return

# -*- coding: utf-8 -*-
#
import numpy


def volume(tet):
    '''Computes the center of the circumsphere of
    '''
    a = tet[1, :] - tet[0, :]
    b = tet[2, :] - tet[0, :]
    c = tet[3, :] - tet[0, :]

    omega = numpy.dot(a, numpy.cross(b, c))

    # https://en.wikipedia.org/wiki/Tetrahedron#Volume
    cell_volumes = abs(omega) / 6.0
    return cell_volumes


def _transform_to_unit_tetrahedron(f, tetrahedron):
    '''Transformation

      x = x0 * N0(xi, eta, zeta) \
        + x1 * N1(xi, eta, zeta) \
        + x2 * N2(xi, eta, zeta) \
        + x3 * N2(xi, eta, zeta)

    with

      N0(xi, eta) = 1 - xi - eta - zeta,
      N1(xi, eta) = xi,
      N2(xi, eta) = eta.
      N3(xi, eta) = zeta.
    '''
    return lambda xi: f(
        + tetrahedron[0] * (1.0 - xi[0] - xi[1] - xi[2])
        + tetrahedron[1] * xi[0]
        + tetrahedron[2] * xi[1]
        + tetrahedron[3] * xi[2]
        )


def integrate(f, tetrahedron, rule):
    g = _transform_to_unit_tetrahedron(f, tetrahedron)
    out = 0.0
    for point, weight in zip(rule.points, rule.weights):
        out += weight * g(point)
    return volume(tetrahedron) * out


class Keast(object):
    '''
    https://people.sc.fsu.edu/~jburkardt/datasets/quadrature_rules_tet/quadrature_rules_tet.html
    '''
    def __init__(self, index):
        if index == 0:
            self.weights = [1.0]
            self.points = numpy.array([
                [0.25, 0.25, 0.25]
                ])
            self.degree = 0
        else:
            raise ValueError('Illegal Keast index')
        return

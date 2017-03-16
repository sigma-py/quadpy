# -*- coding: utf-8 -*-
#
import math
import numpy

from . import helpers


def show(
        scheme,
        wedge=numpy.array([
            [0, 0, 0], [1, 0, 0], [0, 1, 0],
            [0, 0, 1], [1, 0, 1], [0, 1, 1],
            ]),
        show_axes=False
        ):
    '''Shows the quadrature points on a given wedge. The size of the
    balls around the points coincides with their weights.
    '''
    from matplotlib import pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.set_aspect('equal')

    if not show_axes:
        ax.set_axis_off()

    edges = numpy.array([
        [wedge[0], wedge[1]],
        [wedge[1], wedge[2]],
        [wedge[0], wedge[2]],
        #
        [wedge[3], wedge[4]],
        [wedge[4], wedge[5]],
        [wedge[5], wedge[3]],
        #
        [wedge[0], wedge[3]],
        [wedge[1], wedge[4]],
        [wedge[2], wedge[5]],
        ])
    for edge in edges:
        plt.plot(edge[:, 0], edge[:, 1], edge[:, 2], '-k')

    xi = scheme.points[:, 0]
    eta = scheme.points[:, 1]
    zeta = scheme.points[:, 2]
    transformed_pts = \
        + numpy.outer(0.5 * (1.0 - xi - eta)*(1.0 - zeta), wedge[0]) \
        + numpy.outer(0.5 * xi * (1.0 - zeta), wedge[1]) \
        + numpy.outer(0.5 * eta * (1.0 - zeta), wedge[2]) \
        + numpy.outer(0.5 * (1.0 - xi - eta)*(1.0 - zeta), wedge[3]) \
        + numpy.outer(0.5 * xi * (1.0 + zeta), wedge[4]) \
        + numpy.outer(0.5 * eta * (1.0 + zeta), wedge[5])

    vol = integrate(lambda x: numpy.ones(1), wedge, Felippa(1))
    helpers.plot_spheres(
        plt, ax, transformed_pts, scheme.weights, vol,
        wedge[:, 0].min(), wedge[:, 0].max(),
        wedge[:, 1].min(), wedge[:, 1].max(),
        wedge[:, 2].min(), wedge[:, 2].max(),
        )
    plt.show()
    return


def integrate(f, wedge, scheme, sum=helpers.kahan_sum):
    xi = scheme.points.T
    mo = numpy.multiply.outer
    x = \
        + mo(0.5 * (1.0-xi[0]-xi[1]) * (1.0-xi[2]), wedge[0]) \
        + mo(0.5 * xi[0] * (1.0-xi[2]), wedge[1]) \
        + mo(0.5 * xi[1] * (1.0-xi[2]), wedge[2]) \
        + mo(0.5 * (1.0-xi[0]-xi[1]) * (1.0+xi[2]), wedge[3]) \
        + mo(0.5 * xi[0] * (1.0+xi[2]), wedge[4]) \
        + mo(0.5 * xi[1] * (1.0+xi[2]), wedge[5])
    x = x.T
    J0 = \
        - mo(0.5*(1.0 - xi[2]), wedge[0]) \
        + mo(0.5*(1.0 - xi[2]), wedge[2]) \
        - mo(0.5*(1.0 + xi[2]), wedge[3]) \
        + mo(0.5*(1.0 + xi[2]), wedge[5])
    J0 = J0.T
    J1 = \
        - mo(0.5*(1.0 - xi[2]), wedge[0]) \
        + mo(0.5*(1.0 - xi[2]), wedge[1]) \
        - mo(0.5*(1.0 + xi[2]), wedge[3]) \
        + mo(0.5*(1.0 + xi[2]), wedge[4])
    J1 = J1.T
    J2 = \
        - mo(0.5 * (1.0-xi[0]-xi[1]), wedge[0]) \
        - mo(0.5 * xi[0], wedge[1]) \
        - mo(0.5 * xi[1], wedge[2]) \
        + mo(0.5 * (1.0-xi[0]-xi[1]), wedge[3]) \
        + mo(0.5 * xi[0], wedge[4]) \
        + mo(0.5 * xi[1], wedge[5])
    J2 = J2.T
    det = J0[0]*J1[1]*J2[2] + J1[0]*J2[1]*J0[2] + J2[0]*J0[1]*J1[2] \
        - J0[2]*J1[1]*J2[0] - J1[2]*J2[1]*J0[0] - J2[2]*J0[1]*J1[0]

    return sum((scheme.weights * f(x)).T * abs(det.T))


class Felippa(object):
    '''
    Carlos Felippa,
    A compendium of FEM integration formulas for symbolic work,
    Engineering Computation,
    Volume 21, Number 8, 2004, pages 867-890.

    <https://people.sc.fsu.edu/~jburkardt/datasets/quadrature_rules_wedge/quadrature_rules_wedge.html>
    '''
    def __init__(self, index):
        if index == 1:
            self.weights = numpy.array([1.0])
            self.points = numpy.array([
                [1.0/3.0, 1.0/3.0, 0.0]
                ])
            self.degree = 1
        elif index == 2:
            self.weights = 1.0/6.0 * numpy.ones(6)
            self.points = self.s21_z(1.0/6.0, numpy.sqrt(1.0/3.0))
            self.degree = 2
        elif index == 3:
            self.weights = 1.0/6.0 * numpy.ones(6)
            self.points = self.s21_z(0.5, numpy.sqrt(1.0/3.0))
            self.degree = 2
        elif index == 4:
            self.weights = numpy.concatenate([
                0.6205044157722541E-01 * numpy.ones(6),
                0.3054215101536719E-01 * numpy.ones(6),
                0.9928070652356065E-01 * numpy.ones(3),
                0.4886744162458750E-01 * numpy.ones(3),
                ])
            self.points = numpy.concatenate([
                self.s21_z(0.4459484909159649, 0.7745966692414834),
                self.s21_z(0.9157621350977074E-01, 0.7745966692414834),
                self.s21(0.4459484909159649),
                self.s21(0.9157621350977074E-01),
                ])
            self.degree = 4
        elif index == 5:
            self.weights = numpy.concatenate([
                0.3498310570689643E-01 * numpy.ones(6),
                0.3677615355236283E-01 * numpy.ones(6),
                0.6250000000000000E-01 * numpy.ones(2),
                0.5597296913103428E-01 * numpy.ones(3),
                0.5884184568378053E-01 * numpy.ones(3),
                1.0000000000000000E-01 * numpy.ones(1),
                ])
            self.points = numpy.concatenate([
                self.s21_z(0.1012865073234563, 0.7745966692414834),
                self.s21_z(0.4701420641051151, 0.7745966692414834),
                self.s3_z(0.7745966692414834),
                self.s21(0.1012865073234563),
                self.s21(0.4701420641051151),
                self.s3(),
                ])
            self.degree = 5
        else:
            assert index == 6
            self.weights = numpy.concatenate([
                0.8843323515718317E-02 * numpy.ones(6),
                0.2031233592848984E-01 * numpy.ones(6),
                0.1441007403935041E-01 * numpy.ones(12),
                0.1657912966938509E-01 * numpy.ones(6),
                0.3808080193469984E-01 * numpy.ones(6),
                0.2701546376983638E-01 * numpy.ones(12),

                ])
            self.points = numpy.concatenate([
                self.s21_z(0.6308901449150223E-01, -0.8611363115940526),
                self.s21_z(0.2492867451709104, -0.8611363115940526),
                self.s111_z(
                    0.5314504984481695E-01,
                    0.3103524510337844,
                    0.8611363115940526
                    ),
                self.s21_z(0.6308901449150223E-01, 0.3399810435848563),
                self.s21_z(0.2492867451709104, 0.3399810435848563),
                self.s111_z(
                    0.5314504984481695E-01,
                    0.3103524510337844,
                    0.3399810435848563
                    ),
                ])
            self.degree = 6

        return

    def s3(self):
        return numpy.array([
            [1.0/3.0, 1.0/3.0, 0.0],
            ])

    def s3_z(self, z):
        return numpy.array([
            [1.0/3.0, 1.0/3.0, +z],
            [1.0/3.0, 1.0/3.0, -z],
            ])

    def s21(self, a):
        b = 1.0 - 2*a
        return numpy.array([
            [a, b, 0.0],
            [b, a, 0.0],
            [a, a, 0.0],
            ])

    def s21_z(self, a, z):
        b = 1.0 - 2*a
        return numpy.array([
            [a, b, +z],
            [b, a, +z],
            [a, a, +z],
            [a, b, -z],
            [b, a, -z],
            [a, a, -z],
            ])

    def s111_z(self, a, b, z):
        c = 1.0 - a - b
        return numpy.array([
            [b, c, +z],
            [a, b, +z],
            [c, a, +z],
            [c, b, +z],
            [a, c, +z],
            [b, a, +z],
            [b, c, -z],
            [a, b, -z],
            [c, a, -z],
            [c, b, -z],
            [a, c, -z],
            [b, a, -z],
            ])

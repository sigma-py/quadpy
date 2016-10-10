# -*- coding: utf-8 -*-
#
import math
import numpy


def area(radius):
    return 4*numpy.pi*radius**2


def show(scheme):
    from matplotlib import pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.set_aspect('equal')

    phi, theta = numpy.mgrid[0.0:numpy.pi:100j, 0.0:2.0*numpy.pi:100j]
    x = numpy.sin(phi) * numpy.cos(theta)
    y = numpy.sin(phi) * numpy.sin(theta)
    z = numpy.cos(phi)
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
        color='k',
        s=60
        )
    return


def integrate(f, midpoint, radius, rule):
    # w * f(x(xi)) * |det(J)|
    out = math.fsum([
        weight * radius**3 * f(radius*xi + midpoint)
        for xi, weight in zip(rule.points, rule.weights)
        ])
    return area(radius) * out


class Lebedev(object):
    '''
    Sphere integration schemes from

    Lebedev, V. I. (1976),
    Quadratures on a sphere,
    Zh. Vȳchisl. Mat. Mat. Fiz. 16 (2): 293–306.
    doi:10.1016/0041-5553(76)90100-2.

    <https://en.wikipedia.org/wiki/Lebedev_quadrature>
    <https://people.sc.fsu.edu/~jburkardt/datasets/sphere_lebedev_rule/sphere_lebedev_rule.html>
    '''
    def __init__(self, index):
        if index == 1:
            self.weights = numpy.concatenate([
                1.0/6.0 * numpy.ones(6),
                ])
            self.points = numpy.concatenate([
                self.a1(),
                ])
            self.degree = 3
        elif index == 2:
            self.weights = numpy.concatenate([
                1.0/15.0 * numpy.ones(6),
                0.075 * numpy.ones(8),
                ])
            self.points = numpy.concatenate([
                self.a1(),
                self.a3(),
                ])
            self.degree = 5
        elif index == 3:
            self.weights = numpy.concatenate([
                1.0/21.0 * numpy.ones(6),
                4.0/105.0 * numpy.ones(12),
                9.0/280.0 * numpy.ones(8),
                ])
            self.points = numpy.concatenate([
                self.a1(),
                self.a2(),
                self.a3(),
                ])
            self.degree = 7
        elif index == 4:
            self.weights = numpy.concatenate([
                1.0/105.0 * numpy.ones(6),
                9.0/280.0 * numpy.ones(8),
                1.0/35.0 * numpy.ones(24),
                ])
            self.points = numpy.concatenate([
                self.a1(),
                self.a3(),
                self.pq0(8.8807383397711525674e-01, 4.5970084338098304855e-01)
                ])
            self.degree = 9
        elif index == 5:
            self.weights = numpy.concatenate([
                4.0 / 315.0 * numpy.ones(6),
                64.0 / 2835.0 * numpy.ones(12),
                27.0 / 1280.0 * numpy.ones(8),
                14641.0 / 725760.0 * numpy.ones(24),
                ])
            self.points = numpy.concatenate([
                self.a1(),
                self.a2(),
                self.a3(),
                self.llm(3.0151134457776357367e-01, 9.0453403373329088755e-01)
                ])
            self.degree = 11
        elif index == 6:
            self.weights = numpy.concatenate([
                5.1306717973400001e-04 * numpy.ones(6),
                1.6604069565742001e-02 * numpy.ones(12),
                -2.9586038961039000e-02 * numpy.ones(8),
                1.6522170993716001e-02 * numpy.ones(24),
                2.6576207082158999e-02 * numpy.ones(24),
                ])
            self.points = numpy.concatenate([
                self.a1(),
                self.a2(),
                self.a3(),
                self.pq0(3.2077264898077640e-01, 9.4715622136258792e-01),
                self.llm(4.8038446141526142e-01, 7.3379938570534264e-01),
                ])
            self.degree = 13
        elif index == 7:
            self.weights = numpy.concatenate([
                1.1544011544012000e-02 * numpy.ones(6),
                1.1943909085856000e-02 * numpy.ones(8),
                1.1812303746904000e-02 * numpy.ones(24),
                1.1110555710602999e-02 * numpy.ones(24),
                1.1876501294537000e-02 * numpy.ones(24),
                ])
            self.points = numpy.concatenate([
                self.a1(),
                self.a3(),
                self.pq0(3.7424303909034118e-01, 9.2733065715117247e-01),
                self.llm(3.6960284645415015e-01, 8.5251831170126768e-01),
                self.llm(6.9435400660266633e-01, 1.8906355288539489e-01),
                ])
            self.degree = 15
        elif index == 8:
            self.weights = numpy.concatenate([
                3.8282704949370002e-03 * numpy.ones(6),
                9.7937375124880002e-03 * numpy.ones(8),
                9.6949963616630008e-03 * numpy.ones(24),
                8.2117372831910004e-03 * numpy.ones(24),
                9.9428148911780007e-03 * numpy.ones(24),
                9.5954713360710004e-03 * numpy.ones(24),
                ])
            self.points = numpy.concatenate([
                self.a1(),
                self.a3(),
                self.pq0(4.7836902881215021e-01, 8.7815891060406615e-01),
                self.llm(1.8511563534473630e-01, 9.6512403508659406e-01),
                self.llm(6.9042104838229212e-01, 2.1595729184584847e-01),
                self.llm(3.9568947305594182e-01, 8.2876998125259227e-01),
                ])
            self.degree = 17
        elif index == 9:
            self.weights = numpy.concatenate([
                5.9963136886199996e-04 * numpy.ones(6),
                7.3729997186210003e-03 * numpy.ones(12),
                7.2105153601440004e-03 * numpy.ones(8),
                7.1163554931180000e-03 * numpy.ones(24),
                6.7538294863139997e-03 * numpy.ones(24),
                7.5743941590539999e-03 * numpy.ones(24),
                6.9910873533029997e-03 * numpy.ones(48),
                ])
            self.points = numpy.concatenate([
                self.a1(),
                self.a2(),
                self.a3(),
                self.llm(6.7644104001142624e-01, 2.9129888220952688e-01),
                self.llm(4.1749612279654530e-01, 8.0708981835958249e-01),
                self.llm(1.5746766720390826e-01, 9.7488864367717321e-01),
                self.rsw(1.4035538117131835e-01, 4.4933283232695570e-01, 8.8227001126032267e-01),
                ])
            self.degree = 19
        else:
            raise ValueError('Illegal Lebedev index')

        return

    def a1(self):
        return numpy.array([
            [+1.0, 0.0, 0.0],
            [0.0, +1.0, 0.0],
            [0.0, 0.0, +1.0],
            [-1.0, 0.0, 0.0],
            [0.0, -1.0, 0.0],
            [0.0, 0.0, -1.0],
            ])

    def a2(self):
        return numpy.array([
            [+1.0, +1.0, 0.0],
            [+1.0, 0.0, +1.0],
            [0.0, +1.0, +1.0],
            [-1.0, +1.0, 0.0],
            [-1.0, 0.0, +1.0],
            [0.0, -1.0, +1.0],
            [+1.0, -1.0, 0.0],
            [+1.0, 0.0, -1.0],
            [0.0, +1.0, -1.0],
            [-1.0, -1.0, 0.0],
            [-1.0, 0.0, -1.0],
            [0.0, -1.0, -1.0],
            ]) / numpy.sqrt(2.0)

    def a3(self):
        return numpy.array([
            [+1.0, +1.0, +1.0],
            [+1.0, +1.0, -1.0],
            [+1.0, -1.0, +1.0],
            [+1.0, -1.0, -1.0],
            [-1.0, +1.0, +1.0],
            [-1.0, +1.0, -1.0],
            [-1.0, -1.0, +1.0],
            [-1.0, -1.0, -1.0],
            ]) / numpy.sqrt(3.0)

    def pq0(self, p, q):
        assert abs(p**2 + q**2 - 1.0) < 1.0e-12
        return numpy.array([
            [+p, +q, 0.0],
            [+p, 0.0, +q],
            [0.0, +p, +q],
            [+q, +p, 0.0],
            [+q, 0.0, +p],
            [0.0, +q, +p],
            [-p, +q, 0.0],
            [-p, 0.0, +q],
            [0.0, -p, +q],
            [+q, -p, 0.0],
            [+q, 0.0, -p],
            [0.0, +q, -p],
            [+p, -q, 0.0],
            [+p, 0.0, -q],
            [0.0, +p, -q],
            [-q, +p, 0.0],
            [-q, 0.0, +p],
            [0.0, -q, +p],
            [-p, -q, 0.0],
            [-p, 0.0, -q],
            [0.0, -p, -q],
            [-q, -p, 0.0],
            [-q, 0.0, -p],
            [0.0, -q, -p],
            ])

    def llm(self, l, m):
        assert abs(2*l**2 + m**2 - 1.0) < 1.0e-12
        return numpy.array([
            [+l, +l, +m],
            [+l, +m, +l],
            [+m, +l, +l],
            [-l, +l, +m],
            [-l, +m, +l],
            [+m, -l, +l],
            [+l, -l, +m],
            [+l, +m, -l],
            [+m, +l, -l],
            [+l, +l, -m],
            [+l, -m, +l],
            [-m, +l, +l],
            [-l, -l, +m],
            [-l, +m, -l],
            [-m, -l, -l],
            [-l, +l, -m],
            [-l, -m, +l],
            [-m, -l, +l],
            [+l, -l, -m],
            [+l, -m, -l],
            [-m, +l, -l],
            [-l, -l, -m],
            [-l, -m, -l],
            [-m, -l, -l],
            ])

    def rsw(self, r, s, w):
        assert abs(r**2 + s**2 + w**2 - 1.0) < 1.0e-12
        return numpy.array([
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
            ])

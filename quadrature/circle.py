# -*- coding: utf-8 -*-
#
import math
import numpy
import sympy
from . import line


def show(scheme):
    from matplotlib import pyplot as plt
    ax = plt.gca()
    # change default range so that new circles will work
    plt.axis('equal')
    ax.set_xlim((-1.5, 1.5))
    ax.set_ylim((-1.5, 1.5))

    circle1 = plt.Circle((0, 0), 1, color='k', fill=False)
    ax.add_artist(circle1)

    for tp, weight in zip(scheme.points, scheme.weights):
        color = 'b' if weight >= 0 else 'r'
        # highlight circle center
        plt.plot([tp[0]], [tp[1]], '.' + color)
        # plot circle
        # scale the circle volume according to the weight
        radius = numpy.sqrt(abs(weight)/sum(scheme.weights) / numpy.pi)
        circ = plt.Circle((tp[0], tp[1]), radius, color=color, alpha=0.5)
        ax.add_artist(circ)
    return


# def _get_det_J(quad, xi):
#         J0 = \
#             - numpy.outer(quad[0], 0.25*(1-xi[1])) \
#             + numpy.outer(quad[1], 0.25*(1-xi[1])) \
#             + numpy.outer(quad[2], 0.25*(1+xi[1])) \
#             - numpy.outer(quad[3], 0.25*(1+xi[1]))
#         J1 = \
#             - numpy.outer(quad[0], 0.25*(1-xi[0])) \
#             - numpy.outer(quad[1], 0.25*(1+xi[0])) \
#             + numpy.outer(quad[2], 0.25*(1+xi[0])) \
#             + numpy.outer(quad[3], 0.25*(1-xi[0]))
#         det = J0[0]*J1[1] - J1[0]*J0[1]
#         return det


def integrate(f, scheme):
    x = scheme.points.T
    return math.fsum(scheme.weights * f(x).T)


class Peirce(object):
    '''
    W.H. Peirce,
    Numerical integration over the planer annulus,
    J. Soc. Indust. Appl. Math.,
    Vol. 5, No. 2, June, 1957.
    '''
    def __init__(self, m):
        k = 4*m + 3
        self.degree = k
        theta = 2*numpy.pi * numpy.arange(1, k+2) / (k+1)
        p, _ = numpy.polynomial.legendre.leggauss(m+1)
        # scale points to [r0, r1] (where r0 = 0, r1 = 1 for now)
        p = numpy.sqrt(0.5*(p + 1.0))
        p_theta = numpy.dstack(numpy.meshgrid(p, theta)).reshape(-1, 2).T
        self.points = numpy.column_stack([
            p_theta[0] * numpy.cos(p_theta[1]),
            p_theta[0] * numpy.sin(p_theta[1]),
            ])

        # weights
        # C_{j} = 2*pi/(k+1) * B_j
        # where
        #   B_j = 1 / (2*P'_{m+1}(r_j^2) \
        #       * \int_0^1 P_{m+1}(r^2) / (r^2 - r_j^2) dr
        r2 = sympy.Symbol('r2')
        rj2 = sympy.Symbol('rj2')
        pm1 = sympy.prod([
            (r2 - p[i]**2) / p[i]**2 for i in range(m+1)
            ])
        pdiff = sympy.diff(pm1).subs({r2: rj2})

        b = numpy.empty(len(p))
        for j, rj2_val in enumerate(p):
            igrand = sympy.prod([
                (r2 - p[i]**2) / p[i]**2 for i in range(m+1) if i != j
                ])
            igral = sympy.integrate(igrand, (r2, 0, 1))
            b[j] = 1.0 / (2.0 * pdiff.subs({rj2: rj2_val})) * igral

        print(b)

        self.weights = numpy.tile(2*numpy.pi / (k+1) * b, k+1)
        self.weights = numpy.array(16 * [numpy.pi/16.0])

        return

# -*- coding: utf-8 -*-
#
import math
import numpy
import sympy

from . import helpers


def show(
        scheme,
        triangle=numpy.array([
            [-0.5, 0.0],
            [+0.5, 0.0],
            [0, 0.5 * (numpy.sqrt(3))],
            ]),
        show_axes=False
        ):
    '''Shows the quadrature points on a given triangle. The size of the circles
    around the points coincides with their weights.
    '''
    from matplotlib import pyplot as plt

    plt.plot(triangle[:, 0], triangle[:, 1], '-k')
    plt.plot(
        [triangle[-1, 0], triangle[0, 0]],
        [triangle[-1, 1], triangle[0, 1]],
        '-k'
        )

    if not show_axes:
        plt.gca().set_axis_off()

    transformed_pts = \
        + numpy.outer(
            (1.0 - scheme.points[:, 0] - scheme.points[:, 1]),
            triangle[0]
            ) \
        + numpy.outer(scheme.points[:, 0], triangle[1]) \
        + numpy.outer(scheme.points[:, 1], triangle[2])

    vol = integrate(lambda x: numpy.ones(1), triangle, Centroid())
    helpers.plot_disks(
        plt, transformed_pts, scheme.weights, vol
        )

    plt.axis('equal')
    plt.show()
    return


def integrate(f, triangle, scheme, sumfun=helpers.kahan_sum):
    xi = scheme.points.T
    x = \
        + numpy.multiply.outer(1.0 - xi[0] - xi[1], triangle[0]) \
        + numpy.multiply.outer(xi[0], triangle[1]) \
        + numpy.multiply.outer(xi[1], triangle[2])
    x = x.T

    # det is the signed volume of the triangle
    J0 = (triangle[1] - triangle[0]).T
    J1 = (triangle[2] - triangle[0]).T
    # The factor 0.5 is the volume of the reference triangle.
    det = 0.5 * (J0[0]*J1[1] - J1[0]*J0[1])

    return sumfun(scheme.weights * f(x) * abs(det), axis=-1)


def _s3():
    return numpy.array([
        [1.0/3.0, 1.0/3.0, 1.0/3.0]
        ])


def _s21(a):
    b = 1.0 - 2*a
    return numpy.array([
        [a, a, b],
        [a, b, a],
        [b, a, a],
        ])


def _s111(a, b):
    c = 1.0 - a - b
    return numpy.array([
        [a, b, c],
        [c, a, b],
        [b, c, a],
        [b, a, c],
        [c, b, a],
        [a, c, b],
        ])


class Centroid(object):
    def __init__(self):
        self.weights = numpy.array([1.0])
        bary = _s3()
        self.points = bary[:, [1, 2]]
        self.degree = 1
        self.name = 'centroid'
        return


class Vertex(object):
    def __init__(self):
        self.weights = 1.0/3.0 * numpy.ones(3)
        bary = _s21(0.0)
        self.points = bary[:, [1, 2]]
        self.degree = 1
        self.name = 'vertex'
        return


class SevenPoint(object):
    def __init__(self):
        self.weights = numpy.concatenate([
            0.45 * numpy.ones(1),
            0.05 * numpy.ones(3),
            2.0 / 15.0 * numpy.ones(3),
            ])
        bary = numpy.concatenate([
            _s3(),
            _s21(0.0),
            _s21(0.5),
            ])
        self.points = bary[:, [1, 2]]
        self.degree = 3
        self.name = 'seven-point'
        return


class HammerMarloweStroud(object):
    '''
    P.C. Hammer, O.J. Marlowe and A.H. Stroud,
    Numerical Integration Over Simplexes and Cones,
    Mathematical Tables and Other Aids to Computation,
    Vol. 10, No. 55, Jul. 1956, pp. 130-137,
    <https://doi.org/10.1090/S0025-5718-1956-0086389-6>.

    Abstract:
    In this paper we develop numerical integration formulas for simplexes and
    cones in n-space for n>=2. While several papers have been written on
    numerical integration in higher spaces, most of these have dealt with
    hyperrectangular regions. For certain exceptions see [3]. Hammer and Wymore
    [1] have given a first general type theory designed through systematic use
    of cartesian product regions and affine transformations to extend the
    possible usefulness of formulas for each region.

    Two of the schemes also appear in

    P.C. Hammer, Arthur H. Stroud,
    Numerical Evaluation of Multiple Integrals II,
    Mathematical Tables and Other Aids to Computation.
    Vol. 12, No. 64 (Oct., 1958), pp. 272-280,
    <http://www.jstor.org/stable/2002370>
    '''
    def __init__(self, index):
        self.name = 'HMS(%d)' % index
        if index == 1:
            self.weights = numpy.concatenate([
                1.0 * numpy.ones(1),
                ])
            bary = numpy.concatenate([
                _s3(),
                ])
            self.degree = 1
        elif index == 2:
            self.weights = numpy.concatenate([
                1.0/3.0 * numpy.ones(3),
                ])
            bary = numpy.concatenate([
                self._r(0.5),
                ])
            self.degree = 2
        elif index == 3:
            self.weights = numpy.concatenate([
                1.0/3.0 * numpy.ones(3),
                ])
            bary = numpy.concatenate([
                self._r(-0.5),
                ])
            self.degree = 2
        elif index == 4:
            self.weights = numpy.concatenate([
                -9.0/16.0 * numpy.ones(1),
                25.0/48.0 * numpy.ones(3),
                ])
            bary = numpy.concatenate([
                _s3(),
                self._r(0.4),
                ])
            self.degree = 3
        else:
            assert index == 5
            self.weights = numpy.concatenate([
                9.0/40.0 * numpy.ones(1),
                (155.0 - numpy.sqrt(15.0)) / 1200.0 * numpy.ones(3),
                (155.0 + numpy.sqrt(15.0)) / 1200.0 * numpy.ones(3),
                ])
            bary = numpy.concatenate([
                _s3(),
                self._r((1 + numpy.sqrt(15)) / 7.0),
                self._r((1 - numpy.sqrt(15)) / 7.0),
                ])
            self.degree = 5

        self.points = bary[:, 1:]
        return

    def _r(self, r):
        '''Given $r$ (as appearing in the article), it returns the barycentric
        coordinates of the three points.
        '''
        a = r + (1.0-r) / 3.0
        b = 0.5 * (1.0 - a)
        return numpy.array([
            [a, b, b],
            [b, a, b],
            [b, b, a],
            ])


def _newton_cotes(n, point_fun):
    '''
    Construction after

    P. Silvester,
    Symmetric quadrature formulae for simplexes
    Math. Comp., 24, 95-100 (1970),
    <https://doi.org/10.1090/S0025-5718-1970-0258283-6>.
    '''
    degree = n

    # points
    idx = numpy.array([
        [i, j, n-i-j]
        for i in range(n + 1)
        for j in range(n + 1 - i)
        ])
    bary = point_fun(idx, n)
    points = bary[:, [1, 2]]

    # weights
    if n == 0:
        weights = numpy.ones(1)
        return points, weights, degree

    def get_poly(t, m, n):
        return sympy.prod([
            sympy.poly(
                (t - point_fun(k, n)) / (point_fun(m, n) - point_fun(k, n))
                )
            for k in range(m)
            ])
    weights = numpy.empty(len(points))
    idx = 0
    for i in range(n + 1):
        for j in range(n + 1 - i):
            k = n - i - j
            # Define the polynomial which to integrate over the
            # tetrahedron.
            t = sympy.DeferredVector('t')
            g = get_poly(t[0], i, n) \
                * get_poly(t[1], j, n) \
                * get_poly(t[2], k, n)
            # The integral of monomials over a tetrahedron are well-known,
            # see Silvester.
            weights[idx] = numpy.sum([
                 c * numpy.prod([math.factorial(l) for l in m]) * 2
                 / math.factorial(numpy.sum(m) + 2)
                 for m, c in zip(g.monoms(), g.coeffs())
                 ])
            idx += 1
    return points, weights, degree


class NewtonCotesClosed(object):
    def __init__(self, n):
        self.points, self.weights, self.degree = \
            _newton_cotes(n, lambda k, n: k / float(n))
        self.name = 'NCC(%d)' % n
        return


class NewtonCotesOpen(object):
    def __init__(self, n):
        self.points, self.weights, self.degree = \
            _newton_cotes(n, lambda k, n: (k+1) / float(n+3))
        self.name = 'NCO(%d)' % n
        return


class Strang(object):
    '''
    See
    https://people.sc.fsu.edu/~jburkardt/datasets/quadrature_rules_tri/quadrature_rules_tri.html
    and

    Gilbert Strang, George Fix,
    An Analysis of the Finite Element Method,
    Cambridge, 1973,
    ISBN: 096140888X,
    LC: TA335.S77,
    <http://bookstore.siam.org/wc08/>.

    The same schemes are published as

    G.R. Cowper,
    Gaussian quadrature formulas for triangles,
    Numerical Methods in Engineering,
    Volume 7, Issue 3, 1973, Pages 405–408.
    DOI: 10.1002/nme.1620070316,
    <https://dx.doi.org/10.1002/nme.1620070316>.
    '''
    def __init__(self, index):
        self.name = 'Strang(%d)' % index
        if index == 1:
            self.weights = numpy.array([1.0/3.0, 1.0/3.0, 1.0/3.0])
            bary = _s21(1.0/6.0)
            self.degree = 2
        elif index == 2:
            self.weights = numpy.array([1.0/3.0, 1.0/3.0, 1.0/3.0])
            bary = _s21(0.5)
            self.degree = 2
        elif index == 3:
            self.weights = numpy.concatenate([
                -0.5625 * numpy.ones(1),
                25.0 / 48.0 * numpy.ones(3),
                ])
            bary = numpy.concatenate([
                _s3(),
                _s21(0.2),
                ])
            self.degree = 3
        elif index == 4:
            self.weights = 1.0/6.0 * numpy.ones(6)
            bary = _s111(0.659027622374092, 0.231933368553031)
            self.degree = 3
        elif index == 5:
            self.weights = numpy.concatenate([
                0.109951743655322 * numpy.ones(3),
                0.223381589678011 * numpy.ones(3),
                ])
            bary = numpy.concatenate([
                _s21(0.091576213509771),
                _s21(0.445948490915965),
                ])
            self.degree = 4
        elif index == 6:
            self.weights = numpy.concatenate([
                0.375 * numpy.ones(1),
                5.0 / 48.0 * numpy.ones(6),
                ])
            bary = numpy.concatenate([
                _s3(),
                _s111(0.736712498968435, 0.237932366472434),
                ])
            self.degree = 4
        elif index == 7:
            self.weights = numpy.concatenate([
                0.225 * numpy.ones(1),
                0.12593918054482717 * numpy.ones(3),
                0.13239415278850616 * numpy.ones(3),
                ])
            bary = numpy.concatenate([
                _s3(),
                _s21(0.10128650732345633),
                _s21(0.47014206410511505),
                ])
            self.degree = 5
        elif index == 8:
            self.weights = numpy.concatenate([
                0.205950504760887 * numpy.ones(3),
                0.063691414286223 * numpy.ones(6),
                ])
            bary = numpy.concatenate([
                _s21(0.437525248383384),
                _s111(0.797112651860071, 0.165409927389841),
                ])
            self.degree = 5
        elif index == 9:
            self.weights = numpy.concatenate([
                0.050844906370207 * numpy.ones(3),
                0.116786275726379 * numpy.ones(3),
                0.082851075618374 * numpy.ones(6),
                ])
            bary = numpy.concatenate([
                _s21(0.063089014491502),
                _s21(0.249286745170910),
                _s111(0.636502499121399, 0.310352451033785),
                ])
            self.degree = 6
        else:
            assert index == 10
            self.weights = numpy.concatenate([
                -0.149570044467670 * numpy.ones(1),
                0.175615257433204 * numpy.ones(3),
                0.053347235608839 * numpy.ones(3),
                0.077113760890257 * numpy.ones(6),
                ])
            bary = numpy.concatenate([
                _s3(),
                _s21(0.260345966079038),
                _s21(0.065130102902216),
                _s111(0.638444188569809, 0.312865496004875),
                ])
            self.degree = 7

        self.points = bary[:, [1, 2]]
        return


class LynessJespersen(object):
    '''
    J.N. Lyness, D. Jespersen,
    Moderate Degree Symmetric Quadrature Rules for the Triangle,
    J. Inst. Maths Applies (1975) 15, 19-32,
    doi: 10.1093/imamat/15.1.19,
    <https://dx.doi.org/10.1093/imamat/15.1.19>.

    Abstract:
    A variant formulation of the moment fitting equations for the construction
    of D3 (triangularly symmetric) quadrature rules for the triangle is
    derived. These equations are solved to produce weights and abscissas for
    quadrature rules of polynomial degree up to 11 for the triangle, some of
    which require fewer function evaluations than any presently available rule
    of the same polynomial degree. Cytolic rules of degrees up to 9 are also
    derived.
    '''
    def __init__(self, index):
        self.name = 'LJ(%d)' % index
        if index == 1:
            self.weights = numpy.concatenate([
                1.0/3.0 * numpy.ones(3),
                ])
            bary = numpy.concatenate([
                _s21(0.5),
                ])
            self.degree = 2
        elif index == 2:
            self.weights = numpy.concatenate([
                0.75 * numpy.ones(1),
                1.0/12.0 * numpy.ones(3),
                ])
            bary = numpy.concatenate([
                _s3(),
                _s21(0.0),
                ])
            self.degree = 2
        elif index == 3:
            self.weights = numpy.concatenate([
                -9.0/16.0 * numpy.ones(1),
                25.0/48.0 * numpy.ones(3),
                ])
            bary = numpy.concatenate([
                _s3(),
                _s21(0.2),
                ])
            self.degree = 3
        elif index == 4:
            self.weights = numpy.concatenate([
                9.0/20.0 * numpy.ones(1),
                1.0/20.0 * numpy.ones(3),
                2.0/15.0 * numpy.ones(3),
                ])
            bary = numpy.concatenate([
                _s3(),
                _s21(0.0),
                _s21(0.5),
                ])
            self.degree = 3
        elif index == 5:
            self.weights = numpy.concatenate([
                3.298552309659655E-01/3.0 * numpy.ones(3),
                6.701447690340345E-01/3.0 * numpy.ones(3),
                ])
            bary = numpy.concatenate([
                _s21(9.157621350977073E-02),
                _s21(4.459484909159649E-01),
                ])
            self.degree = 4
        elif index == 6:
            self.weights = numpy.concatenate([
                9.0/20.0 * numpy.ones(1),
                -1.0/60.0 * numpy.ones(3),
                1.0/10.0 * numpy.ones(6),
                ])
            bary = numpy.concatenate([
                _s3(),
                _s21(0.0),
                _s111(
                    (3.0 + numpy.sqrt(3.0)) / 6.0,
                    (3.0 - numpy.sqrt(3.0)) / 6.0
                    ),
                ])
            self.degree = 4
        elif index == 7:
            self.weights = numpy.concatenate([
                (11.0 - numpy.sqrt(13.0)) / 360.0 * numpy.ones(3),
                (10.0 - 2*numpy.sqrt(13.0)) / 45.0 * numpy.ones(3),
                (29.0 + 17*numpy.sqrt(13.0)) / 360.0 * numpy.ones(3),
                ])
            bary = numpy.concatenate([
                _s21(0.0),
                _s21(0.5),
                _s21((7.0 - numpy.sqrt(13.0)) / 18.0),
                ])
            self.degree = 4
        elif index == 8:
            self.weights = numpy.concatenate([
                9.0/40.0 * numpy.ones(1),
                (155.0 - numpy.sqrt(15.0)) / 1200.0 * numpy.ones(3),
                (155.0 + numpy.sqrt(15.0)) / 1200.0 * numpy.ones(3),
                ])
            bary = numpy.concatenate([
                _s3(),
                _s21((6.0 - numpy.sqrt(15.0))/21.0),
                _s21((6.0 + numpy.sqrt(15.0))/21.0),
                ])
            self.degree = 5
        elif index == 9:
            self.weights = numpy.concatenate([
                81.0/320.0 * numpy.ones(1),
                1.0/90.0 * numpy.ones(3),
                16.0/225.0 * numpy.ones(3),
                2401.0/14400.0 * numpy.ones(3),
                ])
            bary = numpy.concatenate([
                _s3(),
                _s21(0.0),
                _s21(0.5),
                _s21(1.0/7.0),
                ])
            self.degree = 5
        elif index == 10:
            self.weights = numpy.concatenate([
                3.503588271790222E-01 / 3.0 * numpy.ones(3),
                1.525347191106164E-01 / 3.0 * numpy.ones(3),
                4.971064537103375E-01 / 6.0 * numpy.ones(6),
                ])
            bary = numpy.concatenate([
                _s21(2.492867451709329E-01),
                _s21(6.308901449150177E-02),
                _s111(6.365024991213939E-01, 5.314504984483216E-02),
                ])
            self.degree = 6
        elif index == 11:
            self.weights = numpy.concatenate([
                -81.0/140.0 * numpy.ones(1),
                -5.0/252.0 * numpy.ones(3),
                17.0/315.0 * numpy.ones(3),
                128.0/315.0 * numpy.ones(3),
                9.0/210.0 * numpy.ones(6),
                ])
            bary = numpy.concatenate([
                _s3(),
                _s21(0.0),
                _s21(0.5),
                _s21(0.25),
                _s111(
                    (3.0 + numpy.sqrt(6.0)) / 6.0,
                    (3.0 - numpy.sqrt(6.0)) / 6.0,
                    ),
                ])
            self.degree = 6
        elif index == 12:
            self.weights = numpy.concatenate([
                 1.527089667883523E-01 * numpy.ones(1),
                 2.944076042366762E-01 / 3.0 * numpy.ones(3),
                 3.887052878418766E-01 / 3.0 * numpy.ones(3),
                 1.641781411330949E-01 / 6.0 * numpy.ones(6),
                ])
            bary = numpy.concatenate([
                _s3(),
                _s21(4.738308139536513E-01),
                _s21(1.721176696308175E-01),
                _s111(0.0, 8.653073540834571E-01),
                ])
            self.degree = 6
        elif index == 13:
            self.weights = numpy.concatenate([
                  -1.495700444677495E-01 * numpy.ones(1),
                  5.268457722996328E-01 / 3.0 * numpy.ones(3),
                  1.600417068265167E-01 / 3.0 * numpy.ones(3),
                  4.626825653415500E-01 / 6.0 * numpy.ones(6),
                ])
            bary = numpy.concatenate([
                _s3(),
                _s21(2.603459660790466E-01),
                _s21(6.513010290221623E-02),
                _s111(6.384441885698096E-01, 4.869031542531756E-02),
                ])
            self.degree = 7
        elif index == 14:
            self.weights = numpy.concatenate([
                1.763126156005252E-01 * numpy.ones(1),
                1.210901532763310E-02 / 3.0 * numpy.ones(3),
                3.499561757697094E-01 / 3.0 * numpy.ones(3),
                3.195119754425220E-01 / 3.0 * numpy.ones(3),
                1.421102178595603E-01 / 6.0 * numpy.ones(6),
                ])
            bary = numpy.concatenate([
                _s3(),
                _s21(0.0),
                _s21(1.549360602237604E-01),
                _s21(4.691507461438120E-01),
                _s111(0.0, 8.392991722729236E-01),
                ])
            self.degree = 7
        elif index == 15:
            self.weights = numpy.concatenate([
                1.443156076777862E-01 * numpy.ones(1),
                2.852749028018549E-01 / 3.0 * numpy.ones(3),
                9.737549286959440E-02 / 3.0 * numpy.ones(3),
                3.096521116041552E-01 / 3.0 * numpy.ones(3),
                1.633818850466092E-01 / 6.0 * numpy.ones(6),
                ])
            bary = numpy.concatenate([
                _s3(),
                _s21(4.592925882927229E-01),
                _s21(5.054722831703103E-02),
                _s21(1.705693077517601E-01),
                _s111(8.394777409957211E-03, 7.284923929554041E-01),
                ])
            self.degree = 8
        elif index == 16:
            self.weights = numpy.concatenate([
                +1.207273935292775E-02 / 3.0 * numpy.ones(3),
                -8.491579879151455E-01 / 3.0 * numpy.ones(3),
                +1.042367468891334E+00 / 3.0 * numpy.ones(3),
                +1.947229791412260E-01 / 3.0 * numpy.ones(3),
                +4.511852767201322E-01 / 3.0 * numpy.ones(3),
                +1.488095238055238E-01 / 6.0 * numpy.ones(6),
                ])
            bary = numpy.concatenate([
                _s21(0.0),
                _s21(0.5),
                _s21(4.956813941755582E-01),
                _s21(9.032775751426533E-02),
                _s21(2.341547497073052E-01),
                _s111(0.0, 7.236067977499750E-01),
                ])
            self.degree = 8
        elif index == 17:
            self.weights = numpy.concatenate([
                -2.834183851113958E-01 * numpy.ones(1),
                2.097208857979572E-01 / 3.0 * numpy.ones(3),
                5.127273801480265E-02 / 3.0 * numpy.ones(3),
                6.564896469913508E-01 / 3.0 * numpy.ones(3),
                3.659351143072855E-01 / 6.0 * numpy.ones(6),
                ])
            bary = numpy.concatenate([
                _s3(),
                _s21(4.766654393821525E-01),
                _s21(3.377184405448033E-02),
                _s21(2.703478891654040E-01),
                _s111(5.146433548666149E-02, 7.458294907672514E-01),
                ])
            self.degree = 8
        elif index == 18:
            self.weights = numpy.concatenate([
                9.713579628279610E-02 * numpy.ones(1),
                9.400410068141950E-02 / 3.0 * numpy.ones(3),
                2.334826230143263E-01 / 3.0 * numpy.ones(3),
                2.389432167816271E-01 / 3.0 * numpy.ones(3),
                7.673302697609430E-02 / 3.0 * numpy.ones(3),
                2.597012362637364E-01 / 6.0 * numpy.ones(6),
                ])
            bary = numpy.concatenate([
                _s3(),
                _s21(4.896825191987370E-01),
                _s21(4.370895914929355E-01),
                _s21(1.882035356190322E-01),
                _s21(4.472951339445297E-02),
                _s111(3.683841205473626E-02, 7.411985987844980E-01),
                ])
            self.degree = 9
        elif index == 19:
            self.weights = numpy.concatenate([
                1.133624844599192E-01 * numpy.ones(1),
                1.062573789846330E-03 / 3.0 * numpy.ones(3),
                4.803411513859279E-02 / 3.0 * numpy.ones(3),
                2.524243006337300E-01 / 3.0 * numpy.ones(3),
                7.819254371487040E-02 / 3.0 * numpy.ones(3),
                2.472227459993048E-01 / 3.0 * numpy.ones(3),
                2.597012362637364E-01 / 6.0 * numpy.ones(6),
                ])
            bary = numpy.concatenate([
                _s3(),
                _s21(0.0),
                _s21(0.5),
                _s21(4.497793381870162E-01),
                _s21(4.694744319909033E-02),
                _s21(1.918719127374489E-01),
                _s111(3.683841205473626E-02, 7.411985987844980E-01),
                ])
            self.degree = 9
        elif index == 20:
            self.weights = numpy.concatenate([
                4.097919300803106E-02 / 3.0 * numpy.ones(3),
                1.085536215102866E-01 / 3.0 * numpy.ones(3),
                2.781018986881812E-03 / 3.0 * numpy.ones(3),
                1.779689321422668E-01 / 3.0 * numpy.ones(3),
                2.314486047444677E-01 / 3.0 * numpy.ones(3),
                3.140226717732234E-01 / 6.0 * numpy.ones(6),
                1.242459578348437E-01 / 6.0 * numpy.ones(6),
                ])
            bary = numpy.concatenate([
                _s21(3.236494811127173E-02),
                _s21(1.193509122825931E-01),
                _s21(5.346110482707572E-01),
                _s21(2.033099004312816E-01),
                _s21(3.989693029658558E-01),
                _s111(5.017813831049474E-02, 5.932012134282132E-01),
                _s111(2.102201653616613E-02, 8.074890031597923E-01),
                ])
            self.degree = 11
        else:
            assert index == 21
            self.weights = numpy.concatenate([
                 8.797730116222190E-02 * numpy.ones(1),
                 2.623293466120857E-02 / 3.0 * numpy.ones(3),
                 1.142447159818060E-01 / 3.0 * numpy.ones(3),
                 5.656634416839376E-02 / 3.0 * numpy.ones(3),
                 2.164790926342230E-01 / 3.0 * numpy.ones(3),
                 2.079874161166116E-01 / 3.0 * numpy.ones(3),
                 4.417430269980344E-02 / 6.0 * numpy.ones(6),
                 2.463378925757316E-01 / 6.0 * numpy.ones(6),
                ])
            bary = numpy.concatenate([
                _s3(),
                _s21(2.598914092828833E-02),
                _s21(9.428750264792270E-02),
                _s21(4.946367750172147E-01),
                _s21(2.073433826145142E-01),
                _s21(4.389078057004907E-01),
                _s111(0.0, 8.588702812826364E-01),
                _s111(4.484167758913055E-02, 6.779376548825902E-01),
                ])
            self.degree = 11

        self.points = bary[:, [1, 2]]
        return


class Hillion(object):
    '''
    P. Hillion,
    Numerical Integration on a Triangle,
    International Journal for Numerical Methods in Engineering,
    Vol. 11, 797-815 (1977).
    DOI:10.1002/nme.1620110504,
    <https://dx.doi.org/10.1002/nme.1620110504>.

    Note that the schemes here are not fully symmetric. Also note that in the
    article, the quadrature constants are specified with low precision such
    that the tests are failing. What is needed here is a reimplementation of
    Hillion's method to retrieve more digits.
    '''
    def __init__(self, index):
        self.name = 'Hillion(%d)' % index
        if index == 1:
            self.weights = numpy.concatenate([
                1.0 * numpy.ones(1),
                ])
            self.points = numpy.array([
                [1.0/3.0, 1.0/3.0]
                ])
            self.degree = 1
        elif index == 2:
            self.weights = 2.0 * numpy.concatenate([
                1.0/6.0 * numpy.ones(2),
                1.0/6.0 * numpy.ones(1),
                ])
            self.points = numpy.concatenate([
                self._symm(0.0, 0.5),
                numpy.array([[0.5, 0.5]]),
                ])
            self.degree = 2
        else:
            assert index == 3
            self.weights = 2.0 * numpy.concatenate([
                1.0/6.0 * numpy.ones(2),
                1.0/6.0 * numpy.ones(1),
                ])
            self.points = numpy.concatenate([
                2.0/3.0 - self._symm(0.0, 0.5),
                2.0/3.0 - numpy.array([[0.5, 0.5]]),
                ])
            self.degree = 2
        # elif index == 4:
        #     self.weights = 2.0 * numpy.concatenate([
        #         1.0/18.0 * numpy.ones(1),
        #         2.0/9.0 * numpy.ones(2),
        #         ])
        #     self.points = numpy.concatenate([
        #         numpy.array([[0.0, 0.0]]),
        #         self._symm(0.591506351, 0.158493649),
        #         ])
        #     self.degree = 2
        # elif index == 5:
        #     self.weights = 2.0 * numpy.concatenate([
        #         1.0/18.0 * numpy.ones(1),
        #         2.0/9.0 * numpy.ones(2),
        #         ])
        #     self.points = numpy.concatenate([
        #         2.0/3.0 - numpy.array([[0.0, 0.0]]),
        #         2.0/3.0 - self._symm(0.591506351, 0.158493649),
        #         ])
        #     self.degree = 2
        # elif index == 6:
        #     self.weights = 2.0 * numpy.concatenate([
        #         1.0/8.0 * numpy.ones(4),
        #         ])
        #     lmbda = 0.655308609
        #     mu = 0.247060398
        #     self.points = numpy.concatenate([
        #         self._symm(lmbda, mu),
        #         2.0/3.0 - self._symm(lmbda, mu),
        #         ])
        #     self.degree = 2
        # elif index == 7:
        #     self.weights = 2.0 * numpy.concatenate([
        #         0.159020691 * numpy.ones(2),
        #         0.090979309 * numpy.ones(2),
        #         ])
        #     self.points = numpy.concatenate([
        #         self._symm(0.666390246, 0.280019915),
        #         self._symm(0.178558728, 0.075031109),
        #         ])
        #     self.degree = 3
        # elif index == 8:
        #     self.weights = 2.0 * numpy.concatenate([
        #         0.065104166 * numpy.ones(2),
        #         0.192191138 * numpy.ones(1),
        #         0.177600528 * numpy.ones(1),
        #         ])
        #     lambda2 = 0.433949142
        #     lambda3 = 0.175574667
        #     self.points = numpy.concatenate([
        #         self._symm(0.0, 0.8),
        #         numpy.array([[lambda2, lambda2]]),
        #         numpy.array([[lambda3, lambda3]]),
        #         ])
        #     self.degree = 3
        # elif index == 9:
        #     self.weights = 2.0 * numpy.concatenate([
        #         9.0/32.0 * numpy.ones(1),
        #         25.0/96.0 * numpy.ones(3),
        #         ])
        #     self.points = numpy.concatenate([
        #         numpy.array([[1.0/3.0, 1.0/3.0]]),
        #         self._symm(0.2, 0.6),
        #         numpy.array([[0.2, 0.2]]),
        #         ])
        #     self.degree = 3
        # elif index == 10:
        #     self.weights = 2.0 * numpy.concatenate([
        #         0.036232077 * numpy.ones(2),
        #         0.083559589 * numpy.ones(2),
        #         25.0/96.0 * numpy.ones(1),
        #         ])
        #     self.points = numpy.concatenate([
        #         self._symm(0.939332590, 0.0),
        #         self._symm(0.0, 0.340667409),
        #         numpy.array([[0.4, 0.4]]),
        #         ])
        #     self.degree = 3

        return

    def _symm(self, a, b):
        return numpy.array([
            [a, b],
            [b, a],
            ])


class LaursenGellert(object):
    '''
    M.E. Laursen, M. Gellert,
    Some criteria for numerically integrated matrices and quadrature formulas
    for triangles,
    International Journal for Numerical Methods in Engineering,
    Volume 12, Issue 1, 1978, Pages 67–76.
    DOI: 10.1002/nme.1620120107,
    <https://dx.doi.org/10.1002/nme.1620120107>.

    Abstract:
    For a wide class of finite element matrices integrated numerically rather
    than exactly, a definable number of sampling points is found to be
    sufficient for keeping their theoretical properties unchanged. A systematic
    criterion limiting the number of possible point configurations for
    numerical quadrature formulas on triangles is established. Some new high
    order formulas are presented. Tables containing optimal formulas with
    respect to minimum number of sampling points and required degrees of
    accuracy are given. They are arranged so as to assist with selection of
    suitable quadrature formulas for finite element computer programming.
    '''
    def __init__(self, index):
        self.name = 'LG(%s)' % index
        if index == '1':
            self.weights = numpy.concatenate([
                1.0 * numpy.ones(1),
                ])
            bary = numpy.concatenate([
                _s3(),
                ])
            self.degree = 1
        elif index == '2a':
            self.weights = numpy.concatenate([
                1.0/3.0 * numpy.ones(3),
                ])
            bary = numpy.concatenate([
                _s21(1.0/6.0),
                ])
            self.degree = 2
        elif index == '2b':
            self.weights = numpy.concatenate([
                1.0/3.0 * numpy.ones(3),
                ])
            bary = numpy.concatenate([
                _s21(0.5),
                ])
            self.degree = 2
        elif index == '3':
            self.weights = numpy.concatenate([
                -0.5625 * numpy.ones(1),
                25.0/48.0 * numpy.ones(3),
                ])
            bary = numpy.concatenate([
                _s3(),
                _s21(0.2),
                ])
            self.degree = 3
        elif index == '4':
            self.weights = numpy.concatenate([
                1.0/6.0 * numpy.ones(6),
                ])
            bary = numpy.concatenate([
                _s111(0.659027622374092, 0.231933368553031),
                ])
            self.degree = 3
        elif index == '5':
            self.weights = numpy.concatenate([
                0.109951743655322 * numpy.ones(3),
                0.223381589678011 * numpy.ones(3),
                ])
            bary = numpy.concatenate([
                _s21(0.091576213509771),
                _s21(0.445948490915965),
                ])
            self.degree = 4
        elif index == '6':
            self.weights = numpy.concatenate([
                0.375 * numpy.ones(1),
                5.0/48.0 * numpy.ones(6),
                ])
            bary = numpy.concatenate([
                _s3(),
                _s111(0.736712498968435, 0.237932366472434),
                ])
            self.degree = 4
        elif index == '7':
            self.weights = numpy.concatenate([
                0.225 * numpy.ones(1),
                0.125939180544827 * numpy.ones(3),
                0.132394152788506 * numpy.ones(3),
                ])
            bary = numpy.concatenate([
                _s3(),
                _s21(0.101286507323456),
                _s21(0.470142064105115),
                ])
            self.degree = 5
        elif index == '8':
            self.weights = numpy.concatenate([
                0.205950504760887 * numpy.ones(3),
                0.063691414286223 * numpy.ones(6),
                ])
            bary = numpy.concatenate([
                _s21(0.437525248383384),
                _s111(0.797112651860071, 0.165409927389841),
                ])
            self.degree = 5
        elif index == '9':
            self.weights = numpy.concatenate([
                0.050844906370207 * numpy.ones(3),
                0.116786275726379 * numpy.ones(3),
                0.082851075618374 * numpy.ones(6),
                ])
            bary = numpy.concatenate([
                _s21(0.063089014491502),
                _s21(0.249286745170910),
                _s111(0.636502499121399, 0.310352451033785),
                ])
            self.degree = 6
        elif index == '10':
            self.weights = numpy.concatenate([
                -0.149570044467670 * numpy.ones(1),
                +0.175615257433204 * numpy.ones(3),
                +0.053347235608839 * numpy.ones(3),
                +0.077113760890257 * numpy.ones(6),
                ])
            bary = numpy.concatenate([
                _s3(),
                _s21(0.260345966079038),
                _s21(0.065130102902216),
                _s111(0.638444188569809, 0.312865496004875),
                ])
            self.degree = 7
        elif index == '11':
            self.weights = numpy.concatenate([
                0.053077801790233 * numpy.ones(3),
                0.070853083692136 * numpy.ones(6),
                0.069274682079415 * numpy.ones(6),
                ])
            bary = numpy.concatenate([
                _s21(0.064930513159165),
                _s111(0.284575584249173, 0.517039939069325),
                _s111(0.313559184384932, 0.043863471792371),
                ])
            self.degree = 7
        elif index == '12':
            self.weights = numpy.concatenate([
                0.144315607677787 * numpy.ones(1),
                0.103217370534718 * numpy.ones(3),
                0.032458497623198 * numpy.ones(3),
                0.095091634267284 * numpy.ones(3),
                0.027230314174435 * numpy.ones(6),
                ])
            bary = numpy.concatenate([
                _s3(),
                _s21(0.170569307751761),
                _s21(0.050547228317031),
                _s21(0.459292588292723),
                _s111(0.008394777409958, 0.263112829634638),
                ])
            self.degree = 8
        elif index == '13':
            self.weights = numpy.concatenate([
                0.097135796282799 * numpy.ones(1),
                0.031334700227139 * numpy.ones(3),
                0.077827541004774 * numpy.ones(3),
                0.079647738927210 * numpy.ones(3),
                0.025577675658698 * numpy.ones(3),
                0.043283539377289 * numpy.ones(6),
                ])
            bary = numpy.concatenate([
                _s3(),
                _s21(0.489682519198738),
                _s21(0.437089591492937),
                _s21(0.188203535619033),
                _s21(0.044729513394453),
                _s111(0.036838412054736, 0.221962989160766),
                ])
            self.degree = 9
        elif index == '14':
            self.weights = numpy.concatenate([
                0.051617202569021 * numpy.ones(3),
                0.094080073458356 * numpy.ones(3),
                0.025993571032320 * numpy.ones(3),
                0.045469538047619 * numpy.ones(6),
                0.035351705089199 * numpy.ones(6),
                ])
            bary = numpy.concatenate([
                _s21(0.481519834783311),
                _s21(0.403603979817940),
                _s21(0.045189009784377),
                _s111(0.136991201264904, 0.218290070971381),
                _s111(0.030424361728820, 0.222063165537318),
                ])
            self.degree = 9
        elif index == '15a':
            self.weights = numpy.concatenate([
                0.079894504741240 * numpy.ones(1),
                0.071123802232377 * numpy.ones(3),
                0.008223818690464 * numpy.ones(3),
                0.045430592296170 * numpy.ones(6),
                0.037359856234305 * numpy.ones(6),
                0.030886656884564 * numpy.ones(6),
                ])
            bary = numpy.concatenate([
                _s3(),
                _s21(0.425086210602091),
                _s21(0.023308867510000),
                _s111(0.147925626209534, 0.223766973576973),
                _s111(0.029946031954171, 0.358740141864431),
                _s111(0.035632559587504, 0.143295370426867),
                ])
            self.degree = 10
        else:
            assert index == '15b'
            self.weights = numpy.concatenate([
                0.081743329146286 * numpy.ones(1),
                0.045957963604745 * numpy.ones(3),
                0.013352968813150 * numpy.ones(3),
                0.063904906396424 * numpy.ones(6),
                0.034184648162959 * numpy.ones(6),
                0.025297757707288 * numpy.ones(6),
                ])
            bary = numpy.concatenate([
                _s3(),
                _s21(0.142161101056564),
                _s21(0.032055373216944),
                _s111(0.148132885783821, 0.321812995288835),
                _s111(0.029619889488730, 0.369146781827811),
                _s111(0.028367665339938, 0.163701733737182),
                ])
            self.degree = 10

        self.points = bary[:, 1:]
        return


class Cubtri(object):
    '''
    See
    https://people.sc.fsu.edu/~jburkardt/datasets/quadrature_rules_tri/quadrature_rules_tri.html
    and

    Laurie, D. P.,
    Algorithm 584: CUBTRI: Automatic Cubature over a Triangle,
    ACM Trans. Math. Softw.,
    June 1982,
    <http://dl.acm.org/citation.cfm?id=356001>.
    '''
    def __init__(self):
        self.name = 'CUBTRI'
        self.weights = numpy.concatenate([
            0.0378610912003147 * numpy.ones(1),
            0.0376204254131829 * numpy.ones(3),
            0.0783573522441174 * numpy.ones(3),
            0.1162714796569659 * numpy.ones(3),
            0.0134442673751655 * numpy.ones(3),
            0.0375097224552317 * numpy.ones(6),
            ])

        bary = numpy.concatenate([
            _s3(),
            _s21(0.1012865073234563),
            _s21(0.4701420641051151),
            _s21(0.2321023267750504),
            _s21(0.0294808608844396),
            _s111(0.7384168123405100, 0.2321023267750504),
            ])
        self.points = bary[:, [1, 2]]
        self.degree = 8
        return


class Triex(object):
    '''
    See
    https://people.sc.fsu.edu/~jburkardt/datasets/quadrature_rules_tri/quadrature_rules_tri.html
    and

    E. de Doncker and I. Robinson,
    Algorithm 612: TRIEX: Integration Over a TRIangle Using Nonlinear
    EXtrapolation,
    ACM Trans. Math. Softw.,
    March 1984,
    <http://dl.acm.org/citation.cfm?id=356070&CFID=836775288&CFTOKEN=89206835>.
    '''
    def __init__(self, index):
        self.name = 'TRIEX(%d)' % index
        if index == 19:
            self.weights = numpy.concatenate([
                9.71357962827961025E-002 * numpy.ones(1),
                3.13347002271398278E-002 * numpy.ones(3),
                7.78275410047754301E-002 * numpy.ones(3),
                7.96477389272090969E-002 * numpy.ones(3),
                2.55776756586981006E-002 * numpy.ones(3),
                4.32835393772893970E-002 * numpy.ones(6),
                ])
            bary = numpy.concatenate([
                _s3(),
                _s21(0.48968251919873701),
                _s21(0.43708959149293553),
                _s21(0.18820353561903219),
                _s21(4.47295133944529688E-002),
                _s111(0.74119859878449801, 3.68384120547362581E-002),
                ])
            self.points = bary[:, [1, 2]]
            self.degree = 9
        else:
            assert index == 28
            self.weights = numpy.concatenate([
                0.08797730116222190 * numpy.ones(1),
                0.008744311553736190 * numpy.ones(3),
                0.03808157199393533 * numpy.ones(3),
                0.01885544805613125 * numpy.ones(3),
                0.07215969754474100 * numpy.ones(3),
                0.06932913870553720 * numpy.ones(3),
                0.04105631542928860 * numpy.ones(6),
                0.007362383783300573 * numpy.ones(6),
                ])
            bary = numpy.concatenate([
                _s3(),
                _s21(0.02598914092828833),
                _s21(0.09428750264792270),
                _s21(0.4946367750172147),
                _s21(0.2073433826145142),
                _s21(0.4389078057004907),
                _s111(0.6779376548825902, 0.04484167758913055),
                _s111(0.8588702812826364, 0.0),
                ])
            self.points = bary[:, [1, 2]]
            self.degree = 11

        return


class Dunavant(object):
    '''
    Triangle integration schemes from

    D. A. Dunavant,
    High Degree Efficient Symmetrical Gaussian Quadrature Rules for the
    Triangle,
    Article in International Journal for Numerical Methods in Engineering,
    21(6):1129-1148, June 1985,
    10.1002/nme.1620210612,
    <https://dx.doi.org/10.1002/nme.1620210612>.
    '''
    def __init__(self, index):
        self.name = 'Dunavant(%d)' % index
        if index == 1:
            self.weights = numpy.array([1.0])
            bary = _s3()
            self.degree = 1
        elif index == 2:
            self.weights = 1.0/3.0 * numpy.ones(3)
            bary = _s21(1.0/6.0)
            self.degree = 2
        elif index == 3:
            self.weights = numpy.concatenate([
                -0.5625 * numpy.ones(1),
                25.0 / 48.0 * numpy.ones(3),
                ])
            bary = numpy.concatenate([
                _s3(),
                _s21(0.2),
                ])
            self.degree = 3
        elif index == 4:
            self.weights = numpy.concatenate([
                0.223381589678011 * numpy.ones(3),
                0.109951743655322 * numpy.ones(3),
                ])
            bary = numpy.concatenate([
                _s21(0.445948490915965),
                _s21(0.091576213509771),
                ])
            self.degree = 4
        elif index == 5:
            self.weights = numpy.concatenate([
                0.225 * numpy.ones(1),
                0.132394152788506 * numpy.ones(3),
                0.125939180544827 * numpy.ones(3),
                ])
            bary = numpy.concatenate([
                _s3(),
                _s21(0.4701420641051),
                _s21(0.101286507323456),
                ])
            self.degree = 5
        elif index == 6:
            self.weights = numpy.concatenate([
                0.116786275726379 * numpy.ones(3),
                0.050844906370207 * numpy.ones(3),
                0.082851075618374 * numpy.ones(6),
                ])
            bary = numpy.concatenate([
                _s21(0.249286745170910),
                _s21(0.063089014491502),
                _s111(0.053145049844817, 0.310352451033784),
                ])
            self.degree = 6
        elif index == 7:
            self.weights = numpy.concatenate([
                -0.149570044467682 * numpy.ones(1),
                0.175615257433208 * numpy.ones(3),
                0.053347235608838 * numpy.ones(3),
                0.077113760890257 * numpy.ones(6),
                ])
            bary = numpy.concatenate([
                _s3(),
                _s21(0.260345966079040),
                _s21(0.065130102902216),
                _s111(0.048690315425316, 0.312865496004874),
                ])
            self.degree = 7
        elif index == 8:
            self.weights = numpy.concatenate([
                0.144315607677787 * numpy.ones(1),
                0.095091634267285 * numpy.ones(3),
                0.103217370534718 * numpy.ones(3),
                0.032458497623198 * numpy.ones(3),
                0.027230314174435 * numpy.ones(6),
                ])
            bary = numpy.concatenate([
                _s3(),
                _s21(0.459292588292723),
                _s21(0.170569307751760),
                _s21(0.050547228317031),
                _s111(0.008394777409958, 0.263112829634638),
                ])
            self.degree = 8
        elif index == 9:
            self.weights = numpy.concatenate([
                0.097135796282799 * numpy.ones(1),
                0.031334700227139 * numpy.ones(3),
                0.077827541004774 * numpy.ones(3),
                0.079647738927210 * numpy.ones(3),
                0.025577675658698 * numpy.ones(3),
                0.043283539377289 * numpy.ones(6),
                ])
            bary = numpy.concatenate([
                _s3(),
                _s21(0.489682519198738),
                _s21(0.437089591492937),
                _s21(0.188203535619033),
                _s21(0.044729513394453),
                _s111(0.036838412054736, 0.221962989160766),
                ])
            self.degree = 9
        elif index == 10:
            self.weights = numpy.concatenate([
                0.090817990382754 * numpy.ones(1),
                0.036725957756467 * numpy.ones(3),
                0.045321059435528 * numpy.ones(3),
                0.072757916845420 * numpy.ones(6),
                0.028327242531057 * numpy.ones(6),
                0.009421666963733 * numpy.ones(6),
                ])
            bary = numpy.concatenate([
                _s3(),
                _s21(0.485577633383657),
                _s21(0.109481575485037),
                _s111(0.141707219414880, 0.307939838764121),
                _s111(0.025003534762686, 0.246672560639903),
                _s111(0.009540815400299, 0.066803251012200),
                ])
            self.degree = 10
        elif index == 11:
            self.weights = numpy.concatenate([
                0.000927006328961 * numpy.ones(3),
                0.077149534914813 * numpy.ones(3),
                0.059322977380774 * numpy.ones(3),
                0.036184540503418 * numpy.ones(3),
                0.013659731002678 * numpy.ones(3),
                0.052337111962204 * numpy.ones(6),
                0.020707659639141 * numpy.ones(6),
                ])
            bary = numpy.concatenate([
                _s21(0.534611048270758),
                _s21(0.398969302965855),
                _s21(0.203309900431282),
                _s21(0.119350912282581),
                _s21(0.032364948111276),
                _s111(0.050178138310495, 0.356620648261293),
                _s111(0.021022016536166, 0.171488980304042),
                ])
            self.degree = 11
        elif index == 12:
            self.weights = numpy.concatenate([
                0.025731066440455 * numpy.ones(3),
                0.043692544538038 * numpy.ones(3),
                0.062858224217885 * numpy.ones(3),
                0.034796112930709 * numpy.ones(3),
                0.006166261051559 * numpy.ones(3),
                0.040371557766381 * numpy.ones(6),
                0.022356773202303 * numpy.ones(6),
                0.017316231108659 * numpy.ones(6),
                ])
            bary = numpy.concatenate([
                _s21(0.488217389773805),
                _s21(0.439724392294460),
                _s21(0.271210385012116),
                _s21(0.127576145541586),
                _s21(0.021317350453210),
                _s111(0.115343494534698, 0.275713269685514),
                _s111(0.022838332222257, 0.281325580989940),
                _s111(0.025734050548330, 0.116251915907597),
                ])
            self.degree = 12
        elif index == 13:
            self.weights = numpy.concatenate([
                0.052520923400802 * numpy.ones(1),
                0.011280145209330 * numpy.ones(3),
                0.031423518362454 * numpy.ones(3),
                0.047072502504194 * numpy.ones(3),
                0.047363586536355 * numpy.ones(3),
                0.031167529045794 * numpy.ones(3),
                0.007975771465074 * numpy.ones(3),
                0.036848402728732 * numpy.ones(6),
                0.017401463303822 * numpy.ones(6),
                0.015521786839045 * numpy.ones(6),
                ])
            bary = numpy.concatenate([
                _s3(),
                _s21(0.495048184939705),
                _s21(0.468716635109574),
                _s21(0.414521336801277),
                _s21(0.229399572042831),
                _s21(0.114424495196330),
                _s21(0.024811391363459),
                _s111(0.094853828379579, 0.268794997058761),
                _s111(0.018100773278807, 0.291730066734288),
                _s111(0.022233076674090, 0.126357385491669),
                ])
            self.degree = 13
        elif index == 14:
            self.weights = numpy.concatenate([
                0.021883581369429 * numpy.ones(3),
                0.032788353544125 * numpy.ones(3),
                0.051774104507292 * numpy.ones(3),
                0.042162588736993 * numpy.ones(3),
                0.014433699669777 * numpy.ones(3),
                0.004923403602400 * numpy.ones(3),
                0.024665753212564 * numpy.ones(6),
                0.038571510787061 * numpy.ones(6),
                0.014436308113534 * numpy.ones(6),
                0.005010228838501 * numpy.ones(6),
                ])
            bary = numpy.concatenate([
                _s21(0.488963910362179),
                _s21(0.417644719340454),
                _s21(0.273477528308839),
                _s21(0.177205532412543),
                _s21(0.061799883090873),
                _s21(0.019390961248701),
                _s111(0.057124757403648, 0.172266687821356),
                _s111(0.092916249356972, 0.336861459796345),
                _s111(0.014646950055654, 0.298372882136258),
                _s111(0.001268330932872, 0.118974497696957),
                ])
            self.degree = 14
        elif index == 15:
            self.weights = numpy.concatenate([
                0.001916875642849 * numpy.ones(3),
                0.044249027271145 * numpy.ones(3),
                0.051186548718852 * numpy.ones(3),
                0.023687735870688 * numpy.ones(3),
                0.013289775690021 * numpy.ones(3),
                0.004748916608192 * numpy.ones(3),
                0.038550072599593 * numpy.ones(6),
                0.027215814320624 * numpy.ones(6),
                0.002182077366797 * numpy.ones(6),
                0.021505319847731 * numpy.ones(6),
                0.007673942631049 * numpy.ones(6),
                ])
            bary = numpy.concatenate([
                _s21(0.506972916858243),
                _s21(0.431406354283023),
                _s21(0.277693644847144),
                _s21(0.126464891041254),
                _s21(0.070808385974686),
                _s21(0.018965170241073),
                _s111(0.133734161966621, 0.261311371140087),
                _s111(0.036366677396917, 0.575586555512814),
                _s111(-0.010174883126571, 0.285712220049916),
                _s111(0.036843869875878, 0.215599664072284),
                _s111(0.012459809331199, 0.103575616576386),
                ])
            self.degree = 15
        elif index == 16:
            self.weights = numpy.concatenate([
                0.046875697427642 * numpy.ones(1),
                0.006405878578585 * numpy.ones(3),
                0.041710296739387 * numpy.ones(3),
                0.026891484250064 * numpy.ones(3),
                0.042132522761650 * numpy.ones(3),
                0.030000266842773 * numpy.ones(3),
                0.014200098925024 * numpy.ones(3),
                0.003582462351273 * numpy.ones(3),
                0.032773147460627 * numpy.ones(6),
                0.015298306248441 * numpy.ones(6),
                0.002386244192839 * numpy.ones(6),
                0.019084792755899 * numpy.ones(6),
                0.006850054546542 * numpy.ones(6),
                ])
            bary = numpy.concatenate([
                _s3(),
                _s21(0.497380541948438),
                _s21(0.413469438549352),
                _s21(0.470458599066991),
                _s21(0.240553749969521),
                _s21(0.147965794222573),
                _s21(0.075465187657474),
                _s21(0.016596402623025),
                _s111(0.103575692245252, 0.296555596579887),
                _s111(0.020083411655416, 0.337723063403079),
                _s111(-0.004341002614139, 0.204748281642812),
                _s111(0.041941786468010, 0.189358492130623),
                _s111(0.014317320230681, 0.085283615682657),
                ])
            self.degree = 16
        elif index == 17:
            self.weights = numpy.concatenate([
                0.033437199290803 * numpy.ones(1),
                0.005093415440507 * numpy.ones(3),
                0.014670864527638 * numpy.ones(3),
                0.024350878353672 * numpy.ones(3),
                0.031107550868969 * numpy.ones(3),
                0.031257111218620 * numpy.ones(3),
                0.024815654339665 * numpy.ones(3),
                0.014056073070557 * numpy.ones(3),
                0.003194676173779 * numpy.ones(3),
                0.008119655318993 * numpy.ones(6),
                0.026805742283163 * numpy.ones(6),
                0.018459993210822 * numpy.ones(6),
                0.008476868534328 * numpy.ones(6),
                0.018292796770025 * numpy.ones(6),
                0.006665632004165 * numpy.ones(6),
                ])
            bary = numpy.concatenate([
                _s3(),
                _s21(0.497170540556774),
                _s21(0.482176322624625),
                _s21(0.450239969020782),
                _s21(0.400266239377397),
                _s21(0.252141267970953),
                _s21(0.162047004658461),
                _s21(0.075875882260746),
                _s21(0.015654726967822),
                _s111(0.334319867363658, 0.655493203809423),
                _s111(0.292221537796944, 0.572337590532020),
                _s111(0.319574885423190, 0.626001190286228),
                _s111(0.190704224192292, 0.796427214974071),
                _s111(0.180483211648746, 0.752351005937729),
                _s111(0.080711313679564, 0.904625504095608),
                ])
            self.degree = 17
        elif index == 18:
            self.weights = numpy.concatenate([
                0.030809939937647 * numpy.ones(1),
                0.009072436679404 * numpy.ones(3),
                0.018761316939594 * numpy.ones(3),
                0.019441097985477 * numpy.ones(3),
                0.027753948610810 * numpy.ones(3),
                0.032256225351457 * numpy.ones(3),
                0.025074032616922 * numpy.ones(3),
                0.015271927971832 * numpy.ones(3),
                0.006793922022963 * numpy.ones(3),
                -0.002223098729920 * numpy.ones(3),
                0.006331914076406 * numpy.ones(6),
                0.027257538049138 * numpy.ones(6),
                0.017676785649465 * numpy.ones(6),
                0.018379484638070 * numpy.ones(6),
                0.008104732808192 * numpy.ones(6),
                0.007634129070725 * numpy.ones(6),
                0.000046187660794 * numpy.ones(6),
                ])
            bary = numpy.concatenate([
                _s3(),
                _s21(0.493344808630921),
                _s21(0.469210594241957),
                _s21(0.436281395887006),
                _s21(0.394846170673416),
                _s21(0.249794568803157),
                _s21(0.161432193743843),
                _s21(0.076598227485371),
                _s21(0.024252439353450),
                _s21(0.043146367216965),
                _s111(0.358911494940944, 0.632657968856636),
                _s111(0.294402476751957, 0.574410971510855),
                _s111(0.325017801641814, 0.624779046792512),
                _s111(0.184737559666046, 0.748933176523037),
                _s111(0.218796800013321, 0.769207005420443),
                _s111(0.101179597136408, 0.883962302273467),
                _s111(0.020874755282586, 1.014347260005363),
                ])

            self.degree = 18
        elif index == 19:
            self.weights = numpy.concatenate([
                0.032906331388919 * numpy.ones(1),
                #
                0.010330731891272 * numpy.ones(3),
                0.022387247263016 * numpy.ones(3),
                0.030266125869468 * numpy.ones(3),
                0.030490967802198 * numpy.ones(3),
                0.024159212741641 * numpy.ones(3),
                0.016050803586801 * numpy.ones(3),
                0.008084580261784 * numpy.ones(3),
                0.002079362027485 * numpy.ones(3),
                #
                0.003884876904981 * numpy.ones(6),
                0.025574160612022 * numpy.ones(6),
                0.008880903573338 * numpy.ones(6),
                0.016124546761731 * numpy.ones(6),
                0.002491941817491 * numpy.ones(6),
                0.018242840118951 * numpy.ones(6),
                0.010258563736199 * numpy.ones(6),
                0.003799928855302 * numpy.ones(6),
                ])
            bary = numpy.concatenate([
                _s3(),
                _s21(0.489609987073006),
                _s21(0.454536892697893),
                _s21(0.401416680649431),
                _s21(0.255551654403098),
                _s21(0.177077942152130),
                _s21(0.110061053227952),
                _s21(0.055528624251840),
                _s21(0.012621863777229),
                _s111(0.395754787356943, 0.600633794794645),
                _s111(0.307929983880436, 0.557603261588784),
                _s111(0.264566948406520, 0.720987025817365),
                _s111(0.358539352205951, 0.594527068955871),
                _s111(0.157807405968595, 0.839331473680839),
                _s111(0.075050596975911, 0.701087978926173),
                _s111(0.142421601113383, 0.822931324069857),
                _s111(0.065494628082938, 0.924344252620784),
                ])
            self.degree = 19
        else:
            assert index == 20
            self.weights = numpy.concatenate([
                0.033057055541624 * numpy.ones(1),
                #
                0.000867019185663 * numpy.ones(3),
                0.011660052716448 * numpy.ones(3),
                0.022876936356421 * numpy.ones(3),
                0.030448982673938 * numpy.ones(3),
                0.030624891725355 * numpy.ones(3),
                0.024368057676800 * numpy.ones(3),
                0.015997432032024 * numpy.ones(3),
                0.007698301815602 * numpy.ones(3),
                -0.000632060497488 * numpy.ones(3),
                0.001751134301193 * numpy.ones(3),
                #
                0.016465839189576 * numpy.ones(6),
                0.004839033540485 * numpy.ones(6),
                0.025804906534650 * numpy.ones(6),
                0.008471091054441 * numpy.ones(6),
                0.018354914106280 * numpy.ones(6),
                0.000704404677908 * numpy.ones(6),
                0.010112684927462 * numpy.ones(6),
                0.003573909385950 * numpy.ones(6),
                ])
            bary = numpy.concatenate([
                _s3(),
                _s21(0.500950464352200),
                _s21(0.488212957934729),
                _s21(0.455136681950283),
                _s21(0.401996259318289),
                _s21(0.255892909759421),
                _s21(0.176488255995106),
                _s21(0.104170855336758),
                _s21(0.053068963840930),
                _s21(0.041618715196029),
                _s21(0.011581921406822),
                _s111(0.344855770229001, 0.606402646106160),
                _s111(0.377843269594854, 0.615842614456541),
                _s111(0.306635479062357, 0.559048000390295),
                _s111(0.249419362774742, 0.736606743262866),
                _s111(0.212775724802802, 0.711675142287434),
                _s111(0.146965436053239, 0.861402717154987),
                _s111(0.137726978828923, 0.835586957912363),
                _s111(0.059696109149007, 0.929756171556853),
                ])
            self.degree = 20

        # convert barycentric coordinates to reference triangle
        self.points = bary[:, [1, 2]]
        return


class CoolsHaegemans(object):
    '''
    R. Cools, A. Haegemans,
    Construction of minimal cubature formulae for the square and the triangle
    using invariant theory,
    Department of Computer Science, K.U.Leuven,
    TW Reports vol:TW96, Sept. 1987,
    <https://lirias.kuleuven.be/handle/123456789/131869>.
    '''
    def __init__(self, index):
        self.name = 'CH(%d)' % index
        assert index == 1
        self.weights = 2.0 * numpy.concatenate([
            0.16058343856681218798E-09 * numpy.ones(3),
            0.26530624434780379347E-01 * numpy.ones(3),
            0.29285717640155892159E-01 * numpy.ones(3),
            0.43909556791220782402E-01 * numpy.ones(3),
            0.66940767639916174192E-01 * numpy.ones(3),
            ])
        bary = numpy.concatenate([
            self._r3(
                0.34579201116826902882E+00,
                0.36231682215692616667E+01
                ),
            self._r3(
                0.65101993458939166328E-01,
                0.87016510156356306078E+00
                ),
            self._r3(
                0.65177530364879570754E+00,
                0.31347788752373300717E+00
                ),
            self._r3(
                0.31325121067172530696E+00,
                0.63062143431895614010E+00
                ),
            self._r3(
                0.51334692063945414949E+00,
                0.28104124731511039057E+00
                ),
            ])
        self.degree = 8
        # elif index == 2:
        #     self.weights = 2.0 * numpy.concatenate([
        #         0.15319130036758557631E-06 * numpy.ones(3),
        #         0.13260526227928785221E-01 * numpy.ones(3),
        #         0.15646439344539042136E-01 * numpy.ones(3),
        #         0.21704258224807323311E-01 * numpy.ones(3),
        #         0.21797613600129922367E-01 * numpy.ones(3),
        #         0.38587913508193459468E-01 * numpy.ones(3),
        #         0.39699584282594413022E-01 * numpy.ones(3),
        #         0.47910534861520060665E-01 * numpy.ones(1),
        #         ])
        #     bary = numpy.concatenate([
        #         self._r3(
        #             +0.58469201683584513031E-01,
        #             -0.54887778772527519316E+00
        #             ),
        #         self._r3(
        #             0.50849285064031410705E-01,
        #             0.90799059794957813439E+00
        #             ),
        #         self._r3(
        #             0.51586732419949574487E+00,
        #             0.46312452842927062902E+00
        #             ),
        #         self._r3(
        #             0.24311033191739048230E+00,
        #             0.72180595182371959467E-00
        #             ),
        #         self._r3(
        #             0.75397765920922660134E-00,
        #             0.20647569839132397633E+00
        #             ),
        #         self._r3(
        #             0.42209207910846960294E-00,
        #             0.12689533413411127327E+00
        #             ),
        #         self._r3(
        #             0.19823878346663354068E+00,
        #             0.62124412566393319745E+00
        #             ),
        #         numpy.array([[1.0/3.0, 1.0/3.0, 1.0/3.0]])
        #         ])
        #     self.degree = 10

        self.points = bary[:, 1:]
        return

    def _r3(self, a, b):
        # rotation group R_3
        c = 1.0 - a - b
        return numpy.array([
            [a, b, c],
            [c, a, b],
            [b, c, a],
            ])


class BerntsenEspelid(object):
    '''
    J. Berntsen, T.O. Espelid,
    Degree 13 symmetric quadrature rules for the triangle,
    Reports in Informatics, Dept. of Informatics, University of Bergen,
    (1990).

    Abstract:
    In this paper we develop some tools based on the theory of moments to be
    used to construct triangur symmetric quadrature rules for the triangle. We
    use this technique to construct rules of polynomial degree 13 focusing on
    rules that have all evaluation points inside the triangle, and all weights
    positive.
    '''
    def __init__(self, index):
        self.name = 'BE(%d)' % index
        if index == 1:
            # This first schemes was separately published as
            #
            # Berntsen and Espelid,
            # Algorithm 706: DCUTRI: An Algorithm for Adaptive Cubature over a
            # Collection of Triangles,
            # ACM Trans. Math. Softw.,
            # Sept. 1992,
            # 10.1145/131766.131772,
            # <http://dl.acm.org/citation.cfm?id=131772>.
            self.weights = numpy.concatenate([
                0.051739766065744133555179145422 * numpy.ones(1),
                0.008007799555564801597804123460 * numpy.ones(3),
                0.046868898981821644823226732071 * numpy.ones(3),
                0.046590940183976487960361770070 * numpy.ones(3),
                0.031016943313796381407646220131 * numpy.ones(3),
                0.010791612736631273623178240136 * numpy.ones(3),
                0.032195534242431618819414482205 * numpy.ones(3),
                0.015445834210701583817692900053 * numpy.ones(6),
                0.017822989923178661888748319485 * numpy.ones(6),
                0.037038683681384627918546472190 * numpy.ones(6),
                ])
            bary = numpy.concatenate([
                _s3(),
                _s21(0.024862168537947217274823955239),
                _s21(0.414192542538082326221847602214),
                _s21(0.230293878161404779868453507244),
                _s21(0.113919981661733719124857214943),
                _s21(0.495457300025082323058213517632),
                _s21(0.468861354847056503251458179727),
                _s111(
                    0.022076289653624405142446876931,
                    0.851306504174348550389457672223
                    ),
                _s111(
                    0.018620522802520968955913511549,
                    0.689441970728591295496647976487
                    ),
                _s111(
                    0.096506481292159228736516560903,
                    0.635867859433372768286976979827
                    ),
                ])
        elif index == 2:
            self.weights = numpy.concatenate([
                0.058696079612719031799193912788 * numpy.ones(1),
                0.007850768296100080327451819370 * numpy.ones(3),
                0.050668953175886963421095258589 * numpy.ones(3),
                0.050080326090509066160747699962 * numpy.ones(3),
                0.031647114592298319035326893473 * numpy.ones(3),
                0.005356903791090860889118181848 * numpy.ones(3),
                0.031492563075968795690055730726 * numpy.ones(3),
                0.015802532215260751359123743555 * numpy.ones(6),
                0.015981637780928405322919308674 * numpy.ones(6),
                0.036551502224097295256193503655 * numpy.ones(6),
                ])
            bary = numpy.concatenate([
                _s3(),
                _s21(0.024607188643230218187849951620),
                _s21(0.420308753101194683716920517937),
                _s21(0.227900255506160619646298949779),
                _s21(0.116213058883517905247155308064),
                _s21(0.5),
                _s21(0.476602980049079152951254192421),
                _s111(
                    0.022797894538248612547720754462,
                    0.851775587145410469734660000132
                    ),
                _s111(
                    0.016275770991088540943703616092,
                    0.692797317566660854594116271938
                    ),
                _s111(
                    0.089733060451605359079629076100,
                    0.637955883864209538412552781228
                    ),
                ])
        elif index == 3:
            self.weights = numpy.concatenate([
                -4.438917939249711e-15 * numpy.ones(3),
                0.023875084055169335843543623613 * numpy.ones(3),
                0.063189783598782833129430995388 * numpy.ones(3),
                0.008045069816524589599830031859 * numpy.ones(3),
                0.027856097113552204952023523591 * numpy.ones(3),
                0.050685061067025767745642589150 * numpy.ones(3),
                0.014867088321983380610493967543 * numpy.ones(6),
                0.021575699816275772518477875728 * numpy.ones(6),
                0.043398330702882367361429063273 * numpy.ones(6),
                ])
            bary = numpy.concatenate([
                _s21(-1.097321247106281159287766916114),
                _s21(0.488287850733405315708960134736),
                _s21(0.271000295524474716503595027679),
                _s21(0.024788431033661361058352074973),
                _s21(0.107120353118147709346761786284),
                _s21(0.440323874478061332339068546065),
                _s111(
                    0.020821520846631616958730687380,
                    0.850459062644356742678494398953
                    ),
                _s111(
                    0.022919482804812809947480096117,
                    0.683758575887968213394629723103
                    ),
                _s111(
                    0.115458022821994138042223116054,
                    0.631364930935447484201224031403
                    ),
                ])
        else:
            assert index == 4
            self.weights = numpy.concatenate([
                0.055141401445961668095892272765 * numpy.ones(1),
                0.000011142520455322162070507537 * numpy.ones(3),
                0.008019330681470505488363363198 * numpy.ones(3),
                0.033429216779221783453803543232 * numpy.ones(3),
                0.046966588930899169431852266167 * numpy.ones(3),
                0.031079169485602998741093672276 * numpy.ones(3),
                0.048947942555161210000640851464 * numpy.ones(3),
                0.005884459601338707440236321752 * numpy.ones(3),
                0.015445834210701583817692900053 * numpy.ones(6),
                0.017822989923178661888748319485 * numpy.ones(6),
                0.037038683681384627918546472190 * numpy.ones(6),
                ])
            bary = numpy.concatenate([
                _s3(),
                _s21(0.0),
                _s21(0.024978640633391274114293084881),
                _s21(0.474489920436516855163277733910),
                _s21(0.230836272600280459320993940175),
                _s21(0.114080598593243463483923394518),
                _s21(0.417965185286509715766771174230),
                _s21(0.5),
                _s111(
                    0.022076289653624405142446876931,
                    0.851306504174348550389457672223
                    ),
                _s111(
                    0.018620522802520968955913511549,
                    0.689441970728591295496647976487
                    ),
                _s111(
                    0.096506481292159228736516560903,
                    0.635867859433872768286976979827
                    ),
                ])

        self.degree = 13
        self.points = bary[:, [1, 2]]
        return


class LiuVinokur(object):
    '''
    Y. Liu and M. Vinokur,
    Exact Integrations of Polynomials and Symmetric Quadrature Formulas over
    Arbitrary Polyhedral Grids,
    Journal of Computational Physics, 140, 122–147 (1998).
    DOI: 10.1006/jcph.1998.5884,
    <https://dx.doi.org/10.1006/jcph.1998.5884>.
    '''
    def __init__(self, index):
        self.name = 'LV(%d)' % index
        if index == 1:
            self.weights = numpy.concatenate([
                1.0 * numpy.ones(1),
                ])
            bary = numpy.concatenate([
                _s3(),
                ])
            self.degree = 1
        elif index == 2:
            self.weights = numpy.concatenate([
                1.0/3.0 * numpy.ones(3),
                ])
            bary = numpy.concatenate([
                self._r_alpha(1.0),
                ])
            self.degree = 1
        elif index == 3:
            self.weights = numpy.concatenate([
                1.0/3.0 * numpy.ones(3),
                ])
            bary = numpy.concatenate([
                self._r_alpha(-0.5),
                ])
            self.degree = 2
        elif index == 4:
            self.weights = numpy.concatenate([
                0.75 * numpy.ones(1),
                1.0/12.0 * numpy.ones(3),
                ])
            bary = numpy.concatenate([
                _s3(),
                self._r_alpha(1.0),
                ])
            self.degree = 2
        elif index == 5:
            self.weights = numpy.concatenate([
                -9.0/16.0 * numpy.ones(1),
                25.0/48.0 * numpy.ones(3),
                ])
            bary = numpy.concatenate([
                _s3(),
                # Wrongly specified in the article as 25 (instead of 2/5).
                self._r_alpha(0.4),
                ])
            self.degree = 3
        elif index == 6:
            self.weights = numpy.concatenate([
                (1.0 + numpy.sqrt(21.0)) / 120.0 * numpy.ones(3),
                (39.0 - numpy.sqrt(21.0)) / 120.0 * numpy.ones(3),
                ])
            bary = numpy.concatenate([
                self._r_alpha(1.0),
                self._r_alpha((1.0 - numpy.sqrt(21.0)) / 10.0),
                ])
            self.degree = 3
        elif index == 7:
            self.weights = numpy.concatenate([
                9.0/20.0 * numpy.ones(1),
                1.0/20.0 * numpy.ones(3),
                2.0/15.0 * numpy.ones(3),
                ])
            bary = numpy.concatenate([
                _s3(),
                self._r_alpha(1.0),
                self._r_alpha(-0.5),
                ])
            self.degree = 3
        elif index == 8:
            sqrt10 = numpy.sqrt(10)
            alpha1 = (-10 + 5*sqrt10 + numpy.sqrt(950.0 - 220*sqrt10)) / 30.0
            alpha2 = (-10 + 5*sqrt10 - numpy.sqrt(950.0 - 220*sqrt10)) / 30.0
            self.weights = numpy.concatenate([
                (5*alpha2-2) / (60*alpha1**2 * (alpha2 - alpha1))
                * numpy.ones(3),
                (5*alpha1-2) / (60*alpha2**2 * (alpha1 - alpha2))
                * numpy.ones(3),
                ])
            bary = numpy.concatenate([
                self._r_alpha(alpha1),
                self._r_alpha(alpha2),
                ])
            self.degree = 4
        elif index == 9:
            self.weights = numpy.concatenate([
                27.0/80.0 * numpy.ones(1),
                8.0/105.0 * numpy.ones(3),
                81.0/560.0 * numpy.ones(3),
                ])
            bary = numpy.concatenate([
                _s3(),
                self._r_alpha(-0.5),
                self._r_alpha(2.0/3.0),
                ])
            self.degree = 4
        elif index == 10:
            self.weights = numpy.concatenate([
                (11.0 - numpy.sqrt(13)) / 360.0 * numpy.ones(3),
                (80.0 - 16*numpy.sqrt(13)) / 360.0 * numpy.ones(3),
                (29.0 + 17*numpy.sqrt(13)) / 360.0 * numpy.ones(3),
                ])
            bary = numpy.concatenate([
                self._r_alpha(1.0),
                self._r_alpha(-0.5),
                self._r_alpha((-1.0 + numpy.sqrt(13.0)) / 6.0),
                ])
            self.degree = 4
        elif index == 11:
            self.weights = numpy.concatenate([
                0.45 * numpy.ones(1),
                -1.0 / 60.0 * numpy.ones(3),
                0.1 * numpy.ones(6),
                ])
            bary = numpy.concatenate([
                _s3(),
                self._r_alpha(1.0),
                self._r_gamma_delta(
                    (3.0 + numpy.sqrt(3.0)) / 6.0,
                    (3.0 - numpy.sqrt(3.0)) / 6.0
                    ),
                ])
            self.degree = 4
        elif index == 12:
            self.weights = numpy.concatenate([
                9.0/40.0 * numpy.ones(1),
                (155.0 - numpy.sqrt(15.0))/1200.0 * numpy.ones(3),
                (155.0 + numpy.sqrt(15.0))/1200.0 * numpy.ones(3),
                ])
            bary = numpy.concatenate([
                _s3(),
                self._r_alpha((1.0 + numpy.sqrt(15.0)) / 7.0),
                self._r_alpha((1.0 - numpy.sqrt(15.0)) / 7.0),
                ])
            self.degree = 5
        else:
            assert index == 13
            self.weights = numpy.concatenate([
                81.0/320.0 * numpy.ones(1),
                1.0/90.0 * numpy.ones(3),
                16.0/225.0 * numpy.ones(3),
                2401.0/14400.0 * numpy.ones(3),
                ])
            bary = numpy.concatenate([
                _s3(),
                self._r_alpha(1.0),
                self._r_alpha(-0.5),
                self._r_alpha(4.0/7.0),
                ])
            self.degree = 5

        self.points = bary[:, 1:]
        return

    def _r_alpha(self, alpha):
        '''From the article:

        mu_i = (1 + (n-1) alpha) / n,
        mu_j = (1 - alpha) / n    for j!=i,

        where n is the number of vertices
        '''
        a = (1.0 + 2*alpha) / 3.0
        b = (1.0 - alpha) / 3.0
        return numpy.array([
            [a, b, b],
            [b, a, b],
            [b, b, a],
            ])

    def _r_gamma_delta(self, gamma, delta):
        '''From the article:

        mu_i = (1 + (n-1) gamma - delta) / n,
        mu_j = (1 + (n-1) delta - gamma) / n,
        mu_k = (1 - gamma - delta) / n    for k!=i, k!=j,

        where n is the number of vertices
        '''
        a = (1.0 + 2*gamma - delta) / 3.0
        b = (1.0 + 2*delta - gamma) / 3.0
        c = (1.0 - gamma - delta) / 3.0
        return numpy.array([
            [a, b, c],
            [c, a, b],
            [b, c, a],
            [a, c, b],
            [b, a, c],
            [c, b, a],
            ])


class WandzuraXiao(object):
    '''
    S. Wandzurat (sic!), H. Xiao,
    Symmetric quadrature rules on a triangle,
    Computers & Mathematics with Applications
    Volume 45, Issue 12, June 2003, Pages 1829-1840,
    doi:10.1016/S0898-1221(03)90004-6,
    <https://dx.doi.org/10.1016/S0898-1221(03)90004-6>.

    Abstract:
    We present a class of quadrature rules on triangles in R2 which, somewhat
    similar to Gaussian rules on intervals in R1, have rapid convergence,
    positive weights, and symmetry. By a scheme combining simple group theory
    and numerical optimization, we obtain quadrature rules of this kind up to
    the order 30 on triangles. This scheme, essentially a formalization and
    generalization of the approach used by Lyness and Jespersen over 25 years
    ago, can be easily extended to other regions in R2 and surfaces in higher
    dimensions, such as squares, spheres. We present example formulae and
    relevant numerical results.

    Note that in the above article, the authors present the coordinates in the
    symmetric triangle [[-0.5, -sqrt(3)/2], [-0.5, +sqrt(3)/2], [1, 0]]. These
    have been transformed to barycentric coordinates here.
    '''
    def __init__(self, index):
        self.name = 'WX(%d)' % index
        if index == 1:
            self.weights = numpy.concatenate([
                0.2250000000000000E+00 * numpy.ones(1),
                0.13239415278850623+00 * numpy.ones(3),
                0.12593918054482713+00 * numpy.ones(3),
                ])
            bary = numpy.concatenate([
                _s3(),
                _s21(4.70142064105115109474e-01),
                _s21(1.01286507323456370644e-01),
                ])
            self.degree = 5
        elif index == 2:
            self.weights = numpy.concatenate([
                0.8352339980519638E-01 * numpy.ones(1),
                0.7229850592056743E-02 * numpy.ones(3),
                0.7449217792098051E-01 * numpy.ones(3),
                0.7864647340310853E-01 * numpy.ones(3),
                0.6928323087107504E-02 * numpy.ones(3),
                0.2951832033477940E-01 * numpy.ones(6),
                0.3957936719606124E-01 * numpy.ones(6),
                ])
            bary = numpy.concatenate([
                _s3(),
                _s21(4.978654329544749e-01),
                _s21(4.280124497290562e-01),
                _s21(1.847564127432246e-01),
                _s21(2.048121857167755e-02),
                _s111(1.365735762560334e-01, 3.500298989727198e-02),
                _s111(3.327436005886387e-01, 3.754907025844267e-02),
                ])
            self.degree = 10
        elif index == 3:
            self.weights = numpy.concatenate([
                0.3266181884880529E-01 * numpy.ones(3),
                0.2741281803136436E-01 * numpy.ones(3),
                0.2651003659870330E-01 * numpy.ones(3),
                0.2921596213648611E-01 * numpy.ones(3),
                0.1058460806624399E-01 * numpy.ones(3),
                0.3614643064092035E-02 * numpy.ones(3),
                0.8527748101709436E-02 * numpy.ones(6),
                0.1391617651669193E-01 * numpy.ones(6),
                0.4291932940734835E-02 * numpy.ones(6),
                0.1623532928177489E-01 * numpy.ones(6),
                0.2560734092126239E-01 * numpy.ones(6),
                0.3308819553164567E-01 * numpy.ones(6),
                ])
            bary = numpy.concatenate([
                _s21(4.582807963691251e-01),
                _s21(4.036104645791306e-01),
                _s21(2.931971679130254e-01),
                _s21(1.464677869427729e-01),
                _s21(5.636286766560340e-02),
                _s21(1.657512685837032e-02),
                _s111(2.395345541547944e-01, 9.912203309224784e-03),
                _s111(4.048788073183400e-01, 1.580377063022801e-02),
                _s111(9.500211311304484e-02, 5.143608816970661e-03),
                _s111(1.497531073222740e-01, 4.892232575298879e-02),
                _s111(2.869196124413350e-01, 6.876874863251921e-02),
                _s111(2.818356680990846e-01, 1.684044181246992e-01),
                ])
            self.degree = 15
        elif index == 4:
            self.weights = numpy.concatenate([
                0.2761042699769952E-01 * numpy.ones(1),
                0.1779029547326740E-02 * numpy.ones(3),
                0.2011239811396117E-01 * numpy.ones(3),
                0.2681784725933157E-01 * numpy.ones(3),
                0.2452313380150201E-01 * numpy.ones(3),
                0.1639457841069539E-01 * numpy.ones(3),
                0.1479590739864960E-01 * numpy.ones(3),
                0.4579282277704251E-02 * numpy.ones(3),
                0.1651826515576217E-02 * numpy.ones(3),
                0.2349170908575584E-02 * numpy.ones(6),
                0.4465925754181793E-02 * numpy.ones(6),
                0.6099566807907972E-02 * numpy.ones(6),
                0.6891081327188203E-02 * numpy.ones(6),
                0.7997475072478163E-02 * numpy.ones(6),
                0.7386134285336024E-02 * numpy.ones(6),
                0.1279933187864826E-01 * numpy.ones(6),
                0.1725807117569655E-01 * numpy.ones(6),
                0.1867294590293547E-01 * numpy.ones(6),
                0.2281822405839526E-01 * numpy.ones(6),
                ])
            bary = numpy.concatenate([
                _s3(),
                _s21(4.992496753377855e-01),
                _s21(4.529301240305246e-01),
                _s21(3.977639379552368e-01),
                _s21(2.645002025327873e-01),
                _s21(2.110189640920768e-01),
                _s21(1.077356071712713e-01),
                _s21(3.906908783780266e-02),
                _s21(1.117437972932958e-02),
                _s111(6.354966590835220e-02, 5.349618187337257e-03),
                _s111(1.571069189407070e-01, 7.954817066198938e-03),
                _s111(3.956421143643741e-01, 1.042239828126379e-02),
                _s111(2.731675707129105e-01, 1.096441479612335e-02),
                _s111(1.017853824850171e-01, 3.856671208546233e-02),
                _s111(4.466585491764138e-01, 3.558050781721819e-02),
                _s111(1.990107941495031e-01, 4.967081636276412e-02),
                _s111(3.242611836922827e-01, 5.851972508433175e-02),
                _s111(2.085313632101329e-01, 1.214977870043942e-01),
                _s111(3.231705665362575e-01, 1.407108449439387e-01),
                ])
            self.degree = 20
        elif index == 5:
            self.weights = numpy.concatenate([
                0.8005581880020417E-02 * numpy.ones(3),
                0.1594707683239050E-01 * numpy.ones(3),
                0.1310914123079553E-01 * numpy.ones(3),
                0.1958300096563562E-01 * numpy.ones(3),
                0.1647088544153727E-01 * numpy.ones(3),
                0.8547279074092100E-02 * numpy.ones(3),
                0.8161885857226492E-02 * numpy.ones(3),
                0.6121146539983779E-02 * numpy.ones(3),
                0.2908498264936665E-02 * numpy.ones(3),
                0.6922752456619963E-03 * numpy.ones(3),
                0.1248289199277397E-02 * numpy.ones(6),
                0.3404752908803022E-02 * numpy.ones(6),
                0.3359654326064051E-02 * numpy.ones(6),
                0.1716156539496754E-02 * numpy.ones(6),
                0.1480856316715606E-02 * numpy.ones(6),
                0.3511312610728685E-02 * numpy.ones(6),
                0.7393550149706484E-02 * numpy.ones(6),
                0.7983087477376558E-02 * numpy.ones(6),
                0.4355962613158041E-02 * numpy.ones(6),
                0.7365056701417832E-02 * numpy.ones(6),
                0.1096357284641955E-01 * numpy.ones(6),
                0.1174996174354112E-01 * numpy.ones(6),
                0.1001560071379857E-01 * numpy.ones(6),
                0.1330964078762868E-01 * numpy.ones(6),
                0.1415444650522614E-01 * numpy.ones(6),
                0.1488137956116801E-01 * numpy.ones(6),
                ])
            bary = numpy.concatenate([
                _s21(4.860267584634129e-01),
                _s21(4.344106993361743e-01),
                _s21(3.898891352439638e-01),
                _s21(2.984432340198044e-01),
                _s21(2.340441723373718e-01),
                _s21(1.514683346090176e-01),
                _s21(1.127338935459936e-01),
                _s21(7.771569209152626e-02),
                _s21(3.489309361429704e-02),
                _s21(7.258184620932340e-03),
                _s111(2.272144521533641e-01, 1.292352704442205e-03),
                _s111(4.350105548535717e-01, 5.399701272116197e-03),
                _s111(3.203095992722045e-01, 6.384003033975016e-03),
                _s111(9.175032228000519e-02, 5.028211501993063e-03),
                _s111(3.801083585872434e-02, 6.826758621781874e-03),
                _s111(1.574252184853117e-01, 1.001619963992954e-02),
                _s111(2.398896597785332e-01, 2.575781317338999e-02),
                _s111(3.619431181260605e-01, 3.022789811991581e-02),
                _s111(8.355196095482852e-02, 3.050499010716208e-02),
                _s111(1.484432207324181e-01, 4.595654736256934e-02),
                _s111(2.837397087275350e-01, 6.744280054027756e-02),
                _s111(4.068993751187875e-01, 7.004509141591060e-02),
                _s111(1.941139870248925e-01, 8.391152464011660e-02),
                _s111(3.241343470007031e-01, 1.203755356771527e-01),
                _s111(2.292774835559810e-01, 1.480668991573667e-01),
                _s111(3.256181225959837e-01, 1.917718658673251e-01),
                ])
            self.degree = 25
        else:
            assert index == 6
            self.weights = numpy.concatenate([
                0.1557996020289920E-01 * numpy.ones(1),
                0.3177233700534134E-02 * numpy.ones(3),
                0.1048342663573077E-01 * numpy.ones(3),
                0.1320945957774363E-01 * numpy.ones(3),
                0.1497500696627150E-01 * numpy.ones(3),
                0.1498790444338419E-01 * numpy.ones(3),
                0.1333886474102166E-01 * numpy.ones(3),
                0.1088917111390201E-01 * numpy.ones(3),
                0.8189440660893461E-02 * numpy.ones(3),
                0.5575387588607785E-02 * numpy.ones(3),
                0.3191216473411976E-02 * numpy.ones(3),
                0.1296715144327045E-02 * numpy.ones(3),
                0.2982628261349172E-03 * numpy.ones(3),
                0.9989056850788964E-03 * numpy.ones(6),
                0.4628508491732533E-03 * numpy.ones(6),
                0.1234451336382413E-02 * numpy.ones(6),
                0.5707198522432062E-03 * numpy.ones(6),
                0.1126946125877624E-02 * numpy.ones(6),
                0.1747866949407337E-02 * numpy.ones(6),
                0.1182818815031657E-02 * numpy.ones(6),
                0.1990839294675034E-02 * numpy.ones(6),
                0.1900412795035980E-02 * numpy.ones(6),
                0.4498365808817451E-02 * numpy.ones(6),
                0.3478719460274719E-02 * numpy.ones(6),
                0.4102399036723953E-02 * numpy.ones(6),
                0.4021761549744162E-02 * numpy.ones(6),
                0.6033164660795066E-02 * numpy.ones(6),
                0.3946290302129598E-02 * numpy.ones(6),
                0.6644044537680268E-02 * numpy.ones(6),
                0.8254305856078458E-02 * numpy.ones(6),
                0.6496056633406411E-02 * numpy.ones(6),
                0.9252778144146602E-02 * numpy.ones(6),
                0.9164920726294280E-02 * numpy.ones(6),
                0.1156952462809767E-01 * numpy.ones(6),
                0.1176111646760917E-01 * numpy.ones(6),
                0.1382470218216540E-01 * numpy.ones(6),
                ])
            bary = numpy.concatenate([
                _s3(),
                _s21(4.963349417836173e-01),
                _s21(4.585021620985177e-01),
                _s21(4.245095219372948e-01),
                _s21(3.820470700539167e-01),
                _s21(2.809878457960759e-01),
                _s21(2.273489758540344e-01),
                _s21(1.745591115087298e-01),
                _s21(1.232584272014366e-01),
                _s21(8.008422889219684e-02),
                _s21(4.777446740789878e-02),
                _s21(2.172051468014151e-02),
                _s21(4.764677615437014e-03),
                _s111(4.152952709133117e-01, 9.253711933464581e-04),
                _s111(6.118990978534904e-02, 1.385925855563926e-03),
                _s111(1.649086901369066e-01, 3.682415455910748e-03),
                _s111(2.503506223200247e-02, 3.903223424159349e-03),
                _s111(3.060644651510957e-01, 3.233248155010545e-03),
                _s111(1.070732837302181e-01, 6.467432112236475e-03),
                _s111(2.299575493455843e-01, 3.247475491332678e-03),
                _s111(3.370366333057829e-01, 8.675090806753808e-03),
                _s111(5.625657618206073e-02, 1.559702646731387e-02),
                _s111(4.024513752124010e-01, 1.797672125368521e-02),
                _s111(2.436547020108285e-01, 1.712424535388934e-02),
                _s111(1.653895856145327e-01, 2.288340534658188e-02),
                _s111(9.930187449584685e-02, 3.273759728776667e-02),
                _s111(3.084783330690551e-01, 3.382101234234092e-02),
                _s111(4.606683185921130e-01, 3.554761446001527e-02),
                _s111(2.188152994539296e-01, 5.053979030686654e-02),
                _s111(3.792095515602741e-01, 5.701471491573221e-02),
                _s111(1.429608194181854e-01, 6.415280642120340e-02),
                _s111(2.837312821059250e-01, 8.050114828762560e-02),
                _s111(1.967374410044408e-01, 1.043670681345305e-01),
                _s111(3.558891412116621e-01, 1.138448944287513e-01),
                _s111(2.598186853519115e-01, 1.453634877155237e-01),
                _s111(3.219231812312984e-01, 1.899456528219787e-01),
                ])
            self.degree = 30

        self.points = bary[:, [1, 2]]
        return


class TaylorWingateBos(object):
    '''
    Mark A. Taylor, Beth A. Wingate, Len P. Bos,
    Several new quadrature formulas for polynomial integration in the triangle,
    arXiv,
    Submitted on 27 Jan 2005 (v1), last revised 8 Feb 2007 (this version, v2).

    Abstract:
    We present several new quadrature formulas in the triangle for exact
    integration of polynomials. The points were computed numerically with a
    cardinal function algorithm which imposes that the number of quadrature
    points N be equal to the dimension of a lower dimensional polynomial space.
    Quadrature forumulas are presented for up to degree d=25, all which have
    positive weights and contain no points outside the triangle. Seven of these
    quadrature formulas improve on previously known results.
    '''
    def __init__(self, index):
        self.name = 'TWB(%d)' % index
        if index == 1:
            self.weights = 1.0/3.0 * numpy.ones(3)
            bary = _s21(1.0/6.0)
            self.degree = 2
        elif index == 2:
            self.weights = 0.5 * numpy.concatenate([
                0.2199034873106 * numpy.ones(3),
                0.4467631793560 * numpy.ones(3),
                ])
            bary = numpy.concatenate([
                _s21(0.0915762135098),
                _s21(0.4459484909160),
                ])
            self.degree = 4
        # elif index == 3:
            # not symmetric?
            # self.degree = 5
        elif index == 4:
            self.weights = 0.5 * numpy.concatenate([
                0.0102558174092 * numpy.ones(3),
                0.1116047046647 * numpy.ones(6),
                0.1679775595335 * numpy.ones(3),
                0.2652238803946 * numpy.ones(3),
                ])
            bary = numpy.concatenate([
                _s21(0.0),
                _s111(0.7839656651012, 0.0421382841642),
                _s21(0.4743880861752),
                _s21(0.2385615300181),
                ])
            self.degree = 7
        elif index == 5:
            self.weights = 0.5 * numpy.concatenate([
                 0.0519871420646 * numpy.ones(3),
                 0.0707034101784 * numpy.ones(6),
                 0.0909390760952 * numpy.ones(6),
                 0.1032344051380 * numpy.ones(3),
                 0.1881601469167 * numpy.ones(3),
                ])
            bary = numpy.concatenate([
                _s21(0.0451890097844),
                _s111(0.7475124727339, 0.0304243617288),
                _s111(0.1369912012649, 0.2182900709714),
                _s21(0.4815198347833),
                _s21(0.4036039798179),
                ])
            self.degree = 9
        # elif index == 6:
            # not symmetric?
            # self.degree = 11
        # elif index == 7:
            # not symmetric?
            # self.degree = 13
        else:
            assert index == 8
            self.weights = 0.5 * numpy.concatenate([
                 0.0010616711990 * numpy.ones(3),
                 0.0131460236101 * numpy.ones(6),
                 0.0242881926949 * numpy.ones(6),
                 0.0316799866332 * numpy.ones(6),
                 0.0349317947036 * numpy.ones(3),
                 0.0383664533945 * numpy.ones(3),
                 0.0578369491210 * numpy.ones(6),
                 0.0725821687394 * numpy.ones(6),
                 0.0897856524107 * numpy.ones(3),
                 0.1034544533617 * numpy.ones(3),
                ])
            bary = numpy.concatenate([
                _s21(0.0),
                _s111(0.0573330873026, 0.0151382269814),
                _s111(0.8159625040711, 0.1659719969565),
                _s111(0.3165475556378, 0.0186886898773),
                _s21(0.4903668903754),
                _s21(0.0875134669581),
                _s111(0.0935526036219, 0.2079865423167),
                _s111(0.0974892983467, 0.5380088595149),
                _s21(0.2217145894873),
                _s21(0.3860471669296),
                ])
            self.degree = 14
        # elif index == 9:
            # not symmetric?
            # self.degree = 16
        # elif index == 10:
            # not symmetric?
            # self.degree = 18
        # elif index == 11:
            # not symmetric?
            # self.degree = 20
        # elif index == 12:
        # Not working?
        #     self.weights = 0.5 * numpy.concatenate([
        #         0.0006704436439 * numpy.ones(3),
        #         0.0045472608074 * numpy.ones(6),
        #         0.0052077585320 * numpy.ones(6),
        #         0.0065435432887 * numpy.ones(3),
        #         0.0092737841533 * numpy.ones(6),
        #         0.0095937782623 * numpy.ones(6),
        #         0.0114247809167 * numpy.ones(6),
        #         0.0117216964174 * numpy.ones(6),
        #         0.0188197155232 * numpy.ones(6),
        #         0.0235260980271 * numpy.ones(3),
        #         0.0235571466151 * numpy.ones(3),
        #         0.0268246207430 * numpy.ones(6),
        #         0.0314289776779 * numpy.ones(6),
        #         0.0337196192159 * numpy.ones(6),
        #         0.0427745294213 * numpy.ones(6),
        #         0.0441138932737 * numpy.ones(3),
        #         0.0461469594684 * numpy.ones(6),
        #         0.0469152468624 * numpy.ones(3),
        #         0.0551199980347 * numpy.ones(1),
        #         ])
        #     bary = numpy.concatenate([
        #         _s21(0.0035524391922),
        #         _s111(0.9553548273730, 0.0087898929093),
        #         _s111(0.8865264879047, 0.1082329745017),
        #         _s21(0.0466397432150),
        #         _s111(0.2075720456946, 0.0082759241284),
        #         _s111(0.0858119489725, 0.0314836947701),
        #         _s111(0.6688778233826, 0.0095150760625),
        #         _s111(0.4379999543113, 0.0099859785681),
        #         _s111(0.7974931072148, 0.0405093994119),
        #         _s21(0.3864215551955),
        #         _s21(0.0954935310336),
        #         _s111(0.2745425238718, 0.0479840480721),
        #         _s111(0.4053472446667, 0.5429849622344),
        #         _s111(0.5429849622344, 0.4053472446667),
        #         _s111(0.1195059712009, 0.3057122990643),
        #         _s21(0.2009377128319),
        #         _s111(0.2160775200005, 0.3121360256673),
        #         _s21(0.4376579903849),
        #         _s3(),
        #         ])
        #     self.degree = 21
        # elif index == 13:
            # not symmetric?
            # self.degree = 23
        # elif index == 14:
            # not symmetric?
            # self.degree = 25

        self.points = bary[:, [1, 2]]
        return


class ZhangCuiLiu(object):
    '''
    Linbo Zhang, Tao Cui and Hui Liu,
    A set of symmetric quadrature rules on triangles and tetrahedra,
    Journal of Computational Mathematics
    Vol. 27, No. 1 (January 2009), pp. 89-96,
    <http://www.jstor.org/stable/43693493>.

    Abstract:
    We present a program for computing symmetric quadrature rules on triangles
    and tetrahedra. A set of rules are obtained by using this program.
    Quadrature rules up to order 21 on triangles and up to order 14 on
    tetrahedra have been obtained which are useful for use in finite element
    computations. All rules presented here have positive weights with points
    lying within the integration domain.
    '''
    def __init__(self, index):
        self.name = 'ZCL(%d)' % index
        if index == 1:
            self.weights = numpy.concatenate([
                0.1443156076777871682510911104890646 * numpy.ones(1),
                0.1032173705347182502817915502921290 * numpy.ones(3),
                0.0324584976231980803109259283417806 * numpy.ones(3),
                0.0950916342672846247938961043885843 * numpy.ones(3),
                0.0272303141744349942648446900739089 * numpy.ones(6),
                ])
            bary = numpy.concatenate([
                _s3(),
                _s21(0.1705693077517602066222935014914645),
                _s21(0.0505472283170309754584235505965989),
                _s21(0.4592925882927231560288155144941693),
                _s111(
                    0.2631128296346381134217857862846436,
                    0.0083947774099576053372138345392944
                    ),
                ])
            self.degree = 8
        elif index == 2:
            self.weights = numpy.concatenate([
                0.0585962852260285941278938063477560 * numpy.ones(1),
                0.0017351512297252675680618638808094 * numpy.ones(3),
                0.0261637825586145217778288591819783 * numpy.ones(3),
                0.0039197292424018290965208275701454 * numpy.ones(3),
                0.0122473597569408660972869899262505 * numpy.ones(3),
                0.0281996285032579601073663071515657 * numpy.ones(3),
                0.0508870871859594852960348275454540 * numpy.ones(3),
                0.0504534399016035991910208971341189 * numpy.ones(3),
                0.0170636442122334512900253993849472 * numpy.ones(6),
                0.0096834664255066004075209630934194 * numpy.ones(6),
                0.0363857559284850056220113277642717 * numpy.ones(6),
                0.0069646633735184124253997225042413 * numpy.ones(6),
                ])
            bary = numpy.concatenate([
                _s3(),
                _s21(0.0099797608064584324152935295820524),
                _s21(0.4799778935211883898105528650883899),
                _s21(0.1538119591769669000000000000000000),
                _s21(0.0740234771169878100000000000000000),
                _s21(0.1303546825033300000000000000000000),
                _s21(0.2306172260266531342996053700983831),
                _s21(0.4223320834191478241144087137913939),
                _s111(
                    0.7862373859346610033296221140330900,
                    0.1906163600319009042461432828653034
                    ),
                _s111(
                    0.6305521436606074416224090755688129,
                    0.3623231377435471446183267343597729
                    ),
                _s111(
                    0.6265773298563063142335123137534265,
                    0.2907712058836674150248168174816732
                    ),
                _s111(
                    0.9142099849296254122399670993850469,
                    0.0711657108777507625475924502924336
                    ),
                ])
            self.degree = 14
        else:
            assert index == 3
            self.weights = numpy.concatenate([
                 0.0125376079944966565735856367723948 * numpy.ones(1),
                 0.0274718698764242137484535496073598 * numpy.ones(3),
                 0.0097652722770514230413646914294237 * numpy.ones(3),
                 0.0013984195353918235239233631597867 * numpy.ones(3),
                 0.0092921026251851826304282034030330 * numpy.ones(3),
                 0.0165778760323669253260236250351840 * numpy.ones(3),
                 0.0206677623486650769614219700129729 * numpy.ones(6),
                 0.0208222355211545073068785561993297 * numpy.ones(6),
                 0.0095686384198490606888758450458320 * numpy.ones(6),
                 0.0244527709689724638856439207024089 * numpy.ones(6),
                 0.0031557306306305340038264003207296 * numpy.ones(6),
                 0.0121367963653212969370133090807574 * numpy.ones(6),
                 0.0149664801438864490365249118515707 * numpy.ones(6),
                 0.0063275933217777395693240327504398 * numpy.ones(6),
                 0.0013425603120636958849798512981433 * numpy.ones(6),
                 0.0027760769163475540677293561558015 * numpy.ones(6),
                 0.0107398444741849415551734474479517 * numpy.ones(6),
                 0.0053678057381874532052474100212697 * numpy.ones(6),
                 ])
            bary = numpy.concatenate([
                _s3(),
                _s21(0.2158743059329919731902545438401828),
                _s21(0.0753767665297472780972854309459163),
                _s21(0.0103008281372217921136862160096969),
                _s21(0.4936022112987001655119208321450536),
                _s21(0.4615509381069252967410487102915180),
                _s111(
                    .3286214064242369933034974609509133,
                    .4293405702582103752139588004663984
                    ),
                _s111(
                    .2604803617865687564195930170811535,
                    .1015775342809694461687550061961797
                    ),
                _s111(
                    .1370742358464553000000000000000000,
                    .7100659730011301599879040745464079
                    ),
                _s111(
                    .1467269458722997843041609884874530,
                    .4985454776784148493896226967076119
                    ),
                _s111(
                    .0269989777425532900000000000000000,
                    .0491867226725820016197037125775872
                    ),
                _s111(
                    .0618717859336170268417124700122339,
                    .7796601465405693953603506190768108
                    ),
                _s111(
                    .0477243674276219962083526801042934,
                    .3704915391495476369201496202567388
                    ),
                _s111(
                    .1206005151863643799672337870400794,
                    .8633469487547526484979879960925217
                    ),
                _s111(
                    .0026971477967097876716489145012827,
                    .0561949381877455029878923019865887
                    ),
                _s111(
                    .0030156332779423626572762598234710,
                    .2086750067484213509575944630613577
                    ),
                _s111(
                    .0299053757884570188069287738643386,
                    .7211512409120340910281041502050941
                    ),
                _s111(
                    .0067566542224609885399458175192278,
                    .6400554419405418899040536682721647
                    ),
                ])
            self.degree = 20

        self.points = bary[:, [1, 2]]
        return


class XiaoGimbutas(object):
    '''
    Hong Xiao, Zydrunas Gimbutas,
    A numerical algorithm for the construction of efficient quadrature rules in
    two and higher dimensions,
    Computers & Mathematics with Applications,
    Volume 59, Issue 2, January 2010, Pages 663–676,
    <http://dx.doi.org/10.1016/j.camwa.2009.10.027>.

    Abstract:
    We present a numerical algorithm for the construction of efficient,
    high-order quadratures in two and higher dimensions. Quadrature rules
    constructed via this algorithm possess positive weights and interior nodes,
    resembling the Gaussian quadratures in one dimension. In addition, rules
    can be generated with varying degrees of symmetry, adaptable to individual
    domains. We illustrate the performance of our method with numerical
    examples, and report quadrature rules for polynomials on triangles,
    squares, and cubes, up to degree 50. These formulae are near optimal in the
    number of nodes used, and many of them appear to be new.

    Data adapted from
    <https://people.sc.fsu.edu/~jburkardt/f_src/triangle_symq_rule/triangle_symq_rule.f90>.
    '''
    def __init__(self, index):
        self.name = 'XG(%d)' % index
        if index == 1:
            bary = numpy.concatenate([
                _s3(),
                ])
            self.weights = numpy.concatenate([
                1.0 * numpy.ones(1),
                ])
        elif index == 2:
            bary = numpy.concatenate([
                _s21(0.16666666666666666666666666666667),
                ])
            self.weights = numpy.concatenate([
                0.33333333333333333333333333333333 * numpy.ones(3),
                ])
        elif index == 3:
            bary = numpy.concatenate([
                _s21(0.44594849091596488631832925388305),
                _s21(0.091576213509770743459571463402202),
                ])
            self.weights = numpy.concatenate([
                0.22338158967801146569500700843312 * numpy.ones(3),
                0.10995174365532186763832632490021 * numpy.ones(3),
                ])
        elif index == 4:
            bary = numpy.concatenate([
                _s21(0.44594849091596488631832925388305),
                _s21(0.091576213509770743459571463402202),
                ])
            self.weights = numpy.concatenate([
                0.22338158967801146569500700843312 * numpy.ones(3),
                0.10995174365532186763832632490021 * numpy.ones(3),
                ])
        elif index == 5:
            bary = numpy.concatenate([
                _s21(0.10128650732345633880098736191512),
                _s21(0.47014206410511508977044120951345),
                _s3(),
                ])
            self.weights = numpy.concatenate([
                0.12593918054482715259568394550018 * numpy.ones(3),
                0.13239415278850618073764938783315 * numpy.ones(3),
                0.225 * numpy.ones(1),
                ])
        elif index == 6:
            bary = numpy.concatenate([
                _s21(0.21942998254978296000012662922294),
                _s21(0.4801379641122150440289415322955),
                _s111(
                    0.14161901592396815841177387638623,
                    0.019371724361240788379858517913303
                    ),
                ])
            self.weights = numpy.concatenate([
                0.17133312415298103012021658880582 * numpy.ones(3),
                0.08073108959303097830959370021371 * numpy.ones(3),
                0.040634559793660662451761522156903 * numpy.ones(6),
                ])
        elif index == 7:
            bary = numpy.concatenate([
                _s21(0.47319565368925105124001630452474),
                _s21(0.057797640054506434529631658383363),
                _s111(
                    0.25933901186578568205854125686548,
                    0.046971206130085487780904240865465
                    ),
                _s21(0.24166360639724739069040686292902),
                ])
            self.weights = numpy.concatenate([
                0.053180833296760454560427647904751 * numpy.ones(3),
                0.040918170394056866787122268013906 * numpy.ones(3),
                0.055754540540691087840162442176342 * numpy.ones(6),
                0.127725248561133836305458533062 * numpy.ones(3),
                ])
        elif index == 8:
            bary = numpy.concatenate([
                _s21(0.17056930775176020662229350149146),
                _s21(0.45929258829272315602881551449417),
                _s3(),
                _s21(0.050547228317030975458423550596599),
                _s111(
                    0.26311282963463811342178578628464,
                    0.0083947774099576053372138345392966
                    ),
                ])
            self.weights = numpy.concatenate([
                0.10321737053471825028179155029213 * numpy.ones(3),
                0.095091634267284624793896104388585 * numpy.ones(3),
                0.14431560767778716825109111048907 * numpy.ones(1),
                0.032458497623198080310925928341781 * numpy.ones(3),
                0.02723031417443499426484469007391 * numpy.ones(6),
                ])
        elif index == 9:
            bary = numpy.concatenate([
                _s21(0.48968251919873762778370692483619),
                _s3(),
                _s21(0.18820353561903273024096128046734),
                _s111(
                    0.22196298916076569567510252769319,
                    0.036838412054736283634817598783386
                    ),
                _s21(0.43708959149293663726993036443535),
                _s21(0.044729513394452709865106589966278),
                ])
            self.weights = numpy.concatenate([
                0.03133470022713907053685483128721 * numpy.ones(3),
                0.097135796282798833819241982507292 * numpy.ones(1),
                0.079647738927210253032891774264047 * numpy.ones(3),
                0.04328353937728937728937728937729 * numpy.ones(6),
                0.077827541004774279316739356299405 * numpy.ones(3),
                0.025577675658698031261678798559 * numpy.ones(3),
                ])
        elif index == 10:
            bary = numpy.concatenate([
                _s21(0.49517345980117049293764079093536),
                _s21(0.019139415242841232124217710336254),
                _s111(
                    0.13373475510086914136425983641841,
                    0.034723620482327441447814775265889
                    ),
                _s3(),
                _s111(
                    0.32669313628133689449055887255581,
                    0.037582727341191656547456899984084
                    ),
                _s21(0.184485012685246491320535450702),
                _s21(0.4282348209437188701355815327906),
                ])
            self.weights = numpy.concatenate([
                0.0097925904984183028163296675923584 * numpy.ones(3),
                0.0063853592301186533643020285096736 * numpy.ones(3),
                0.028962281463256338983864796880846 * numpy.ones(6),
                0.083614874373973922670067589784662 * numpy.ones(1),
                0.038739049086018899324602317319021 * numpy.ones(6),
                0.078633769746377276933363606180449 * numpy.ones(3),
                0.075247327968543982712381272722906 * numpy.ones(3),
                ])
        elif index == 11:
            bary = numpy.concatenate([
                _s21(0.03084689563558803473172144500853),
                _s21(0.49878016517846077698471569037401),
                _s111(
                    0.159303619837693496248052920664,
                    0.014366662569555593251324747847985
                    ),
                _s3(),
                _s21(0.11320782728669394677971888159308),
                _s21(0.43665501639317610910521443993827),
                _s21(0.21448345861926932473774314809805),
                _s111(
                    0.31063121631346312612263084490447,
                    0.047664066972150744009649231386716
                    ),
                ])
            self.weights = numpy.concatenate([
                0.012249296950707962508685782317507 * numpy.ones(3),
                0.012465491873881379968438177032217 * numpy.ones(3),
                0.014557623337809246861550792470953 * numpy.ones(6),
                0.081445134709351286204184323120947 * numpy.ones(1),
                0.040129242381308313537086313491967 * numpy.ones(3),
                0.063094872159898675064938725063184 * numpy.ones(3),
                0.067845107743695139643324822340212 * numpy.ones(3),
                0.040642848655886470076514910219681 * numpy.ones(6),
                ])
        elif index == 12:
            bary = numpy.concatenate([
                _s21(0.27146250701492608487799594802255),
                _s21(0.10925782765935429058374731200431),
                _s21(0.44011164865859311101345370307291),
                _s111(
                    0.25545422863851734653135572748307,
                    0.1162960196779265866307481341497
                    ),
                _s111(
                    0.12727971723358936878784187474738,
                    0.021382490256170589594154514820557
                    ),
                _s111(
                    0.29165567973834096053391381193394,
                    0.023034156355267139481641501928932
                    ),
                _s21(0.48820375094554155177800812386881),
                _s21(0.024646363436335594766731197214324),
                ])
            self.weights = numpy.concatenate([
                0.06254121319590276046932378144552 * numpy.ones(3),
                0.028486052068877544999694473907515 * numpy.ones(3),
                0.049918334928060942119051722060995 * numpy.ones(3),
                0.04322736365941421054913327091686 * numpy.ones(6),
                0.015083677576511438585850590462767 * numpy.ones(6),
                0.021783585038607557932631436339292 * numpy.ones(6),
                0.024266838081452033150717972575608 * numpy.ones(3),
                0.0079316425099736384593147879058598 * numpy.ones(3),
                ])
        elif index == 13:
            bary = numpy.concatenate([
                _s21(0.49613589474104609045695539592253),
                _s21(0.46960868965349191238533212311234),
                _s21(0.23111028494908224276146574912332),
                _s111(
                    0.29207868857663640526462753629005,
                    0.018988004383759018504130968446279
                    ),
                _s3(),
                _s111(
                    0.26674525331035119231416830685932,
                    0.097736031066016514304542996851621
                    ),
                _s21(0.41447757027905457442041850715054),
                _s21(0.11355991257213317272595867060872),
                _s111(
                    0.12679977578383735601052059313016,
                    0.021966344206529192279566209528841
                    ),
                _s21(0.024895931491216385389981473879154),
                ])
            self.weights = numpy.concatenate([
                0.009941476361072588041789045547843 * numpy.ones(3),
                0.032781241603722970240012998128976 * numpy.ones(3),
                0.046062409592778250776496057583897 * numpy.ones(3),
                0.01812549864620088072889653109465 * numpy.ones(6),
                0.051622646664290813070987066955006 * numpy.ones(1),
                0.037211960457261537984544519052895 * numpy.ones(6),
                0.046947095542155186331202975178876 * numpy.ones(3),
                0.03090309797575979114247508630053 * numpy.ones(3),
                0.015393072683782175152765652752479 * numpy.ones(6),
                0.0080293997952584213786147424748364 * numpy.ones(3),
                ])
        elif index == 14:
            bary = numpy.concatenate([
                _s21(0.41764471934045392250944082218564),
                _s111(
                    0.29837288213625775297083151805961,
                    0.014646950055654409670541327920073
                    ),
                _s21(0.061799883090872601267478828436936),
                _s111(
                    0.33686145979634500174405519708893,
                    0.09291624935697182475824858954872
                    ),
                _s21(0.2734775283088386597549442832627),
                _s21(0.17720553241254343695661069046506),
                _s21(0.01939096124870104817825009505454),
                _s21(0.48896391036217863867737602045239),
                _s111(
                    0.17226668782135557837528960161366,
                    0.057124757403647939035677124218913
                    ),
                _s111(
                    0.1189744976969568453981819619299,
                    0.0012683309328720250872464010954932
                    ),
                ])
            self.weights = numpy.concatenate([
                0.032788353544125350641310978738626 * numpy.ones(3),
                0.014436308113533840496088691999016 * numpy.ones(6),
                0.014433699669776667601709921480653 * numpy.ones(3),
                0.038571510787060683228489027810411 * numpy.ones(6),
                0.051774104507291586314784910166398 * numpy.ones(3),
                0.042162588736993017538230437324186 * numpy.ones(3),
                0.0049234036024000816818260235090422 * numpy.ones(3),
                0.021883581369428890640844945963326 * numpy.ones(3),
                0.024665753212563673962875245183637 * numpy.ones(6),
                0.0050102288385006717698600930824892 * numpy.ones(6),
                ])
        elif index == 15:
            bary = numpy.concatenate([
                _s21(0.12997822993307786699561917269956),
                _s3(),
                _s21(0.4600769492970597257053156775539),
                _s111(
                    0.18232178340719133268537229161241,
                    0.084594221482191775263412202029957
                    ),
                _s111(
                    0.15020038406523876215738569668549,
                    0.016027089786345402681672035469101
                    ),
                _s111(
                    0.3231113151637126779928254051729,
                    0.097650442430242313818837437676215
                    ),
                _s21(0.49168581663029724013927558778923),
                _s21(0.22153234079514198906474838336725),
                _s21(0.39693373740906058790144135569373),
                _s111(
                    0.30794768148367287186227817077653,
                    0.018454251904633121970689825962834
                    ),
                _s21(0.056341917696100096257056132266807),
                _s111(
                    0.038035229301109283846095501988108,
                    0.0011135352740136955272965228742432
                    ),
                ])
            self.weights = numpy.concatenate([
                0.0073975040670460990270507910097438 * numpy.ones(3),
                0.029730419748071321921339990936424 * numpy.ones(1),
                0.021594087936438450651427293397615 * numpy.ones(3),
                0.024230008783125606042212940063436 * numpy.ones(6),
                0.011228504298878057474685258866626 * numpy.ones(6),
                0.031075220470510947420754385438807 * numpy.ones(6),
                0.015832276350021797046061171800196 * numpy.ones(3),
                0.046287286105198075256693754502225 * numpy.ones(3),
                0.046336041391207233100957981029885 * numpy.ones(3),
                0.016436762092827891766792893356416 * numpy.ones(6),
                0.015084474247597067156477285721111 * numpy.ones(3),
                0.0024752660145579158559970517215904 * numpy.ones(6),
                ])
        elif index == 16:
            bary = numpy.concatenate([
                _s111(
                    0.41376948582708518207145189607931,
                    0.0096649544036601953182082922101098
                    ),
                _s111(
                    0.30417944822947972367002083274056,
                    0.030305943355186323776981456256415
                    ),
                _s21(0.066674472240238279376313692073284),
                _s111(
                    0.089609089022705868291418742903643,
                    0.010812972776103686948420403353724
                    ),
                _s111(
                    0.296615372400382997741575070111,
                    0.10665316053614843395449857111644
                    ),
                _s21(0.24132168070137835085171818387765),
                _s21(0.41279809595522368530012767799294),
                _s111(
                    0.16976335515028977415819553167341,
                    0.051354315344013092580421269283027
                    ),
                _s21(0.15006373658703512444564473227114),
                _s111(
                    0.21404877992584728424070996916893,
                    0.0036969427073555421358471671923928
                    ),
                _s21(0.46954803099668497462765090400142),
                _s3(),
                _s21(0.017041629405718390626196889609321),
                ])
            self.weights = numpy.concatenate([
                0.0081822105532221371452165419336272 * numpy.ones(6),
                0.013983607124653564572845464305998 * numpy.ones(6),
                0.012425425595561008336910473694903 * numpy.ones(3),
                0.0057518699704971587666798025660687 * numpy.ones(6),
                0.031646061681983245867023842159643 * numpy.ones(6),
                0.041184041069792549230560837167821 * numpy.ones(3),
                0.040985219786815363638137069316155 * numpy.ones(3),
                0.017653081047103284436546741351738 * numpy.ones(6),
                0.028783496702748910891928652496484 * numpy.ones(3),
                0.0046146906397291343640099082519866 * numpy.ones(6),
                0.027093669467710448397907318813214 * numpy.ones(3),
                0.046227910314191341490827042152998 * numpy.ones(1),
                0.0037891352382642220363020333223085 * numpy.ones(3),
                ])
        elif index == 17:
            bary = numpy.concatenate([
                _s21(0.417103444361599201605122916372),
                _s111(
                    0.072505470799002426274525642459242,
                    0.011575175903180615348042118401617
                    ),
                _s21(0.18035811626637062018584307222938),
                _s111(
                    0.41547545929522905484233417495494,
                    0.013229672760086893307016713349105
                    ),
                _s111(
                    0.27179187005535484351747174019367,
                    0.013135870834002694966603639429968
                    ),
                _s111(
                    0.29921894247697032086774068478606,
                    0.1575054779268699050276109155237
                    ),
                _s21(0.28570650243658662800183787874469),
                _s111(
                    0.30628159174618654153530949385487,
                    0.067349377867361196395968203347312
                    ),
                _s111(
                    0.16872251349525945693224173036633,
                    0.078042340568282424171330870519312
                    ),
                _s21(0.066654063479596929756766748398501),
                _s111(
                    0.15919228747279267872330207968591,
                    0.016017642362119296842407312133267
                    ),
                _s21(0.014755491660753953403154179300809),
                _s21(0.46559787161889030189589159084694),
                ])
            self.weights = numpy.concatenate([
                0.027310926528102107523778847554031 * numpy.ones(3),
                0.0045843484017358668429496073707797 * numpy.ones(6),
                0.026312630588017984956034415173581 * numpy.ones(3),
                0.010398439955839536463628565378088 * numpy.ones(6),
                0.0086922145010011915687739254404959 * numpy.ones(6),
                0.026171625935336987257122893805802 * numpy.ones(6),
                0.037716237152795280016428805967956 * numpy.ones(3),
                0.022487772546691066442322155780169 * numpy.ones(6),
                0.020557898320454517496005829581422 * numpy.ones(6),
                0.012459000802305442095250042865414 * numpy.ones(3),
                0.0079783002059295933242826221450654 * numpy.ones(6),
                0.0027738875776376421546124925923297 * numpy.ones(3),
                0.025019450950497357797057530176382 * numpy.ones(3),
                ])
        elif index == 18:
            bary = numpy.concatenate([
                _s111(
                    0.38504403441316367334400254247437,
                    0.090427040354340612427383113261328
                    ),
                _s21(0.47491821132404573588789755091754),
                _s21(0.15163850697260486492387353795772),
                _s111(
                    0.047276141832651782522284038985054,
                    0.012498932483495440128048193579532
                    ),
                _s3(),
                _s111(
                    0.30206195771287080772484323648551,
                    0.054011735339024234680444362470851
                    ),
                _s111(
                    0.25650615977424154068897765977749,
                    0.010505018819241935598686033442105
                    ),
                _s21(0.41106710187591949855469549486746),
                _s111(
                    0.17847912556588763355267204638677,
                    0.066122458028403387700539471853985
                    ),
                _s21(0.26561460990537421478430796115175),
                _s21(0.0037589443410683458570246273328608),
                _s111(
                    0.26857330639601384733212028806857,
                    0.1490669101257738392001911394479
                    ),
                _s111(
                    0.41106566867461836291309677848251,
                    0.011691824674667085270423426497857
                    ),
                _s111(
                    0.13277883027138932992144407050471,
                    0.014331524778941953568448671295635
                    ),
                _s21(0.072438705567332870474262063744802),
                ])
            self.weights = numpy.concatenate([
                0.015328258194553140867046286819207 * numpy.ones(6),
                0.013107027491738755678601531003486 * numpy.ones(3),
                0.020318338845458397305216768560988 * numpy.ones(3),
                0.0042175167747444429098438771600713 * numpy.ones(6),
                0.0307485212391158553993533382016 * numpy.ones(1),
                0.016365908413986565958152216113746 * numpy.ones(6),
                0.007729835280006227008092796341026 * numpy.ones(6),
                0.033471994059847898118769734621442 * numpy.ones(3),
                0.016911653917480078794565533238269 * numpy.ones(6),
                0.031116396602006131196893892501586 * numpy.ones(3),
                0.00053200561694778056109294261721747 * numpy.ones(3),
                0.027592886488579478020095933346207 * numpy.ones(6),
                0.00958612447436150376044024017261 * numpy.ones(6),
                0.0076417049727196359508471137212573 * numpy.ones(6),
                0.0137902866047669388014726908033 * numpy.ones(3),
                ])
        elif index == 19:
            bary = numpy.concatenate([
                _s111(
                    0.14242228257112693843245883040253,
                    0.0050051423523504110778304506166263
                    ),
                _s21(0.052526279854103524808062704663702),
                _s111(
                    0.060083899962702342608937343673389,
                    0.0097770614386768393836590377412627
                    ),
                _s111(
                    0.13070066996053455915725076932197,
                    0.039142449434608843634133960160205
                    ),
                _s111(
                    0.3113183832239868349383174472771,
                    0.12931280976797897536964522732728
                    ),
                _s21(0.11144805571699863784701750730831),
                _s21(0.011639027327922553950714476224703),
                _s21(0.25516213315312480457667696158388),
                _s111(
                    0.22143394188911346418996437516731,
                    0.074561189304355101555510307213326
                    ),
                _s21(0.40396971796638608940856849621176),
                _s111(
                    0.35402592699971188743478682850316,
                    0.040888314464978101885612505318782
                    ),
                _s21(0.17817100607962746695652260484561),
                _s21(0.45919438895682762351548822881514),
                _s3(),
                _s21(0.4925124498658742341616474622229),
                _s111(
                    0.24189410400689262254316444034362,
                    0.014923638907438453942330405213637
                    ),
                _s111(
                    0.36462041433871005975346220121911,
                    0.0020691038491023166754900870236664
                    ),
                ])
            self.weights = numpy.concatenate([
                0.0029256924878800713098553783949684 * numpy.ones(6),
                0.0071093936227949467800810622029638 * numpy.ones(3),
                0.0033273888405939044557062050203018 * numpy.ones(6),
                0.0096955190816242013097359494290231 * numpy.ones(6),
                0.02634626470744536227647499553445 * numpy.ones(6),
                0.015234956517004834387172491255861 * numpy.ones(3),
                0.0017651924183085399617460072817281 * numpy.ones(3),
                0.031752854587529979653570928084442 * numpy.ones(3),
                0.018108074590430505045143705048976 * numpy.ones(6),
                0.031537358645239615180042918165454 * numpy.ones(3),
                0.016102209460939429742159784802699 * numpy.ones(6),
                0.024651981053584831413480786435791 * numpy.ones(3),
                0.022983570977123249830659946994427 * numpy.ones(3),
                0.034469160850905275376758916895335 * numpy.ones(1),
                0.010321882182418862490961700458231 * numpy.ones(3),
                0.0084559248390934794657511081039159 * numpy.ones(6),
                0.0032821375148397373168551337436636 * numpy.ones(6),
                ])
        elif index == 20:
            bary = numpy.concatenate([
                _s21(0.18629499774454094277173250805763),
                _s21(0.037310880598884693993503716287207),
                _s21(0.47624561154049901321035187223259),
                _s111(
                    0.064090585608434060049865667617829,
                    0.0048549376076237532883651515789124
                    ),
                _s111(
                    0.21560705739009440302944796795049,
                    0.1062272047202700423625543020539
                    ),
                _s21(0.44555105695592481518722227843421),
                _s111(
                    0.15913370765706722456461535142835,
                    0.007570780504696528263233458528561
                    ),
                _s111(
                    0.31786012383577201597596767163026,
                    0.13980807199179989859415148396842
                    ),
                _s21(0.25457926767333911238485264812738),
                _s3(),
                _s111(
                    0.19851813222878817838002170300785,
                    0.046560364907664316513769554605267
                    ),
                _s21(0.39342534781709985867896374350328),
                _s111(
                    0.099952296288138655382907823347114,
                    0.038363684775374595372611352447152
                    ),
                _s111(
                    0.42002375881622407964644291694727,
                    0.0098315482928025605785112357286863
                    ),
                _s111(
                    0.33313481730958748942762509766254,
                    0.054987479142986810634227248515814
                    ),
                _s111(
                    0.28058141142366523483748653054406,
                    0.010737212856011087330524561007898
                    ),
                _s21(0.010976141028397764039648011858869),
                _s21(0.10938359671171459701687247489983),
                ])
            self.weights = numpy.concatenate([
                0.018346925948505828959626427059069 * numpy.ones(3),
                0.0043225508213311550523463610789034 * numpy.ones(3),
                0.01420365060681688098384023296224 * numpy.ones(3),
                0.0022597392042517310589617575291965 * numpy.ones(6),
                0.015445215644198459688740825325427 * numpy.ones(6),
                0.0189047998664648953736570813204 * numpy.ones(3),
                0.0044057948371169951280656354078927 * numpy.ones(6),
                0.023383491463655473865092155468892 * numpy.ones(6),
                0.028166402615040495064320892062365 * numpy.ones(3),
                0.027820221402906231815783394126949 * numpy.ones(1),
                0.0119727971579093800402926280725 * numpy.ones(6),
                0.027576101258140918026775445431639 * numpy.ones(3),
                0.0082914230552277154255247253178153 * numpy.ones(6),
                0.0073913630005105959564384372263853 * numpy.ones(6),
                0.017334451134438666124365313029018 * numpy.ones(6),
                0.007156400476915370729862673215586 * numpy.ones(6),
                0.0015976815821332397239472923665329 * numpy.ones(3),
                0.015660461552149066842203501824447 * numpy.ones(3),
                ])
        elif index == 21:
            bary = numpy.concatenate([
                _s21(0.29893623531498257064359772308441),
                _s111(
                    0.28918949607859476265178213045511,
                    0.2052955593351615050362430442403
                    ),
                _s21(0.49700787546868558053252834048897),
                _s111(
                    0.23787338259799400007987661964201,
                    0.0069318090314680546433631007229396
                    ),
                _s21(0.40361758654638511638909196596921),
                _s111(
                    0.31886531079482824349012448174201,
                    0.12377940040549275351294332498264
                    ),
                _s111(
                    0.23187362537040097344930256078915,
                    0.038991362623220328367410435043323
                    ),
                _s111(
                    0.13316712294137028349253207861886,
                    0.009536247529710560208892415762726
                    ),
                _s111(
                    0.34680797980991109291701749655779,
                    0.053052191701216797924664686360099
                    ),
                _s21(0.11898857762271940719491988769835),
                _s21(0.19028871809127849772894603455961),
                _s111(
                    0.21659962318998254434769089503197,
                    0.10045802007411443233285291232308
                    ),
                _s21(0.48159786865321661952255739159726),
                _s21(0.44981279177536239217046290397329),
                _s111(
                    0.12882980796205153822433172696037,
                    0.049451065568540524908252963750219
                    ),
                _s21(0.053627575546144939337769469477774),
                _s111(
                    0.36095340801892220726094318830771,
                    0.010254635872924476042110655609965
                    ),
                _s21(0.010742456432828395187390974191145),
                _s111(
                    0.055719565072371984276290560539928,
                    0.010301903643423833683094036414693
                    ),
                ])
            self.weights = numpy.concatenate([
                0.021451121929132340492889362227647 * numpy.ones(3),
                0.017495416155763123021303138559512 * numpy.ones(6),
                0.0044378296970658790383403090296877 * numpy.ones(3),
                0.0042061202881497290911101965620686 * numpy.ones(6),
                0.023000704653283862600844227306397 * numpy.ones(3),
                0.018447484847932833533759924384536 * numpy.ones(6),
                0.010469904185324844538190538000319 * numpy.ones(6),
                0.0044808131219014752414327802853715 * numpy.ones(6),
                0.014500305918971019948688199549555 * numpy.ones(6),
                0.013656032452230197729482481889348 * numpy.ones(3),
                0.01945524186075071049855868631961 * numpy.ones(3),
                0.01590403670542797286779859639871 * numpy.ones(6),
                0.012214410163384381933819953520748 * numpy.ones(3),
                0.019614475227824020274460093128843 * numpy.ones(3),
                0.0098119718225504107741390354414909 * numpy.ones(6),
                0.0071520851012836513519050029718477 * numpy.ones(3),
                0.0068398848579343058904332242076064 * numpy.ones(6),
                0.0015086992723786891139071205207421 * numpy.ones(3),
                0.0032654285840440852427074148200636 * numpy.ones(6),
                ])
        elif index == 22:
            bary = numpy.concatenate([
                _s21(0.38518452462730214123319240996523),
                _s21(0.45776941136767207640932772763919),
                _s111(
                    0.069842169467443620335145946540006,
                    0.0078762822215823575258769613065786
                    ),
                _s21(0.29455825902995012612850997981184),
                _s21(0.1885105236302838921183610560885),
                _s21(0.42198188879353495701319985757606),
                _s111(
                    0.090398831166407777876305279663081,
                    0.044752284348335851078030442383568
                    ),
                _s111(
                    0.41134176402055874902748172759845,
                    0.038275234700863778355102261426328
                    ),
                _s111(
                    0.33210610500744639162553114408612,
                    0.10274707598693136100019092523082
                    ),
                _s111(
                    0.36257628043246727776894407843762,
                    0.0074002412347106898705588429129972
                    ),
                _s111(
                    0.29006682411666881865087514557371,
                    0.19108129796672010632565718043561
                    ),
                _s21(0.49616117840970864625343538363109),
                _s111(
                    0.28793180282417184703524300028589,
                    0.043991645393455821384417518876891
                    ),
                _s111(
                    0.21678693336494116733654375214589,
                    0.10868994186267196803882615499311
                    ),
                _s111(
                    0.14587371987352517588843360150048,
                    0.0091447113749640379493310892692486
                    ),
                _s111(
                    0.17629743482450005895373068221636,
                    0.048254924114641353001610975187971
                    ),
                _s111(
                    0.24399064603949309015835974948164,
                    0.0091639092481851589089713652940034
                    ),
                _s21(0.029108470670807505348851241236437),
                _s21(0.11543153821920495401933025088712),
                _s111(
                    0.01793432105293900640642624485353,
                    0.0017984649889483413284386278162501
                    ),
                ])
            self.weights = numpy.concatenate([
                0.013493083883610660009406780553498 * numpy.ones(3),
                0.013861399524234191485624638578506 * numpy.ones(3),
                0.0025954384742312777262088758558782 * numpy.ones(6),
                0.021075763957452182357605355581786 * numpy.ones(3),
                0.016021299125148891598848867667707 * numpy.ones(3),
                0.018853092553841287696668830967479 * numpy.ones(3),
                0.0075175778177883762121860713931544 * numpy.ones(6),
                0.0111973134719627695701315471314 * numpy.ones(6),
                0.017719093489510217389489043676352 * numpy.ones(6),
                0.0049042603975569638755759701153663 * numpy.ones(6),
                0.02170641955550895569886443536384 * numpy.ones(6),
                0.0052893396659844186994995227530383 * numpy.ones(3),
                0.011662222867343002859586400182817 * numpy.ones(6),
                0.015710162622570316548498743419387 * numpy.ones(6),
                0.0041066870715755557400711362137619 * numpy.ones(6),
                0.010563584967746898041783996760287 * numpy.ones(6),
                0.0050540768975846013973328165987688 * numpy.ones(6),
                0.0035691091658563765156663553713728 * numpy.ones(3),
                0.014415713128104603332904719583625 * numpy.ones(3),
                0.00064042853117142575882509442715151 * numpy.ones(6),
                ])
        elif index == 23:
            bary = numpy.concatenate([
                _s111(
                    0.15950379892475727415245345439626,
                    0.023870253654353568145989397165732
                    ),
                _s111(
                    0.11410136032236457465857255918642,
                    0.0051898217608444771157820321250864
                    ),
                _s21(0.039007268757032081338961413036489),
                _s111(
                    0.095539878171734936251830053893987,
                    0.032741029188706355465001046356455
                    ),
                _s111(
                    0.31116226805170197720486448816666,
                    0.0024475998559663065672699201626314
                    ),
                _s111(
                    0.20561723205805207173864427144621,
                    0.0087252895853085094127617966521407
                    ),
                _s111(
                    0.047261629449725311018960000702985,
                    0.0071625399102444347353662323320398
                    ),
                _s111(
                    0.35850959356962509915759656771972,
                    0.068526954187212964349612688916122
                    ),
                _s21(0.48032887733730851585672483583372),
                _s111(
                    0.24048277203501271124730468088991,
                    0.10172832932728420582705117968801
                    ),
                _s21(0.086841048207633144209733518872567),
                _s21(0.39432350601154154945222196625564),
                _s21(0.26625131787724729725073948861588),
                _s111(
                    0.17293230312922398454015135446875,
                    0.058351575237515396068776691486705
                    ),
                _s111(
                    0.31630430765383809669876232906715,
                    0.15483015540551621605503727622331
                    ),
                _s21(0.1371293873116476212950221464774),
                _s21(0.49895943120958633412436522792224),
                _s111(
                    0.39775857680300767100249353600535,
                    0.014758969729945113484256970638379
                    ),
                _s21(0.44469244212772752248677495233215),
                _s21(0.19874980639653623072964999268399),
                _s111(
                    0.27879416981410228425076371995414,
                    0.032993708192532783107825481478511
                    ),
                _s21(0.0090164402055983186982591212039678),
                _s3(),
                ])
            self.weights = numpy.concatenate([
                0.0025281660553822629672074087883227 * numpy.ones(6),
                0.0022250197297245146982342815356651 * numpy.ones(6),
                0.0039157402590329353410701464318236 * numpy.ones(3),
                0.0053280304311947847437074609178299 * numpy.ones(6),
                0.0022811036762558344816005716937097 * numpy.ones(6),
                0.00411475034441609163469364874074 * numpy.ones(6),
                0.0019525913278907262766784861684304 * numpy.ones(6),
                0.014981113393199167384539739220532 * numpy.ones(6),
                0.011397889267800758912968413920949 * numpy.ones(3),
                0.016121241637017152849030375731228 * numpy.ones(6),
                0.0089599170255135424579356891752058 * numpy.ones(3),
                0.023674608463128020582350901652009 * numpy.ones(3),
                0.02380786288749976010679067189595 * numpy.ones(3),
                0.01047025649313006499243731020335 * numpy.ones(6),
                0.020844395858968808472718831254232 * numpy.ones(6),
                0.01455944939274174623428865499317 * numpy.ones(3),
                0.0024075446041814094752379167408305 * numpy.ones(3),
                0.0070977788345218238282808725566654 * numpy.ones(6),
                0.018951950669338883871080760285167 * numpy.ones(3),
                0.019935277880105022211225727688395 * numpy.ones(3),
                0.010175574656707035965842040604179 * numpy.ones(6),
                0.0010653612328293149864178021030952 * numpy.ones(3),
                0.025253060323036207692073780850926 * numpy.ones(1),
                ])
        elif index == 24:
            bary = numpy.concatenate([
                _s3(),
                _s21(0.41889097491060276780569831520914),
                _s21(0.16236063371692634521241934570799),
                _s111(
                    0.24147976007359402415226382105443,
                    0.17036728246244368117207889225294
                    ),
                _s111(
                    0.32897580892422636485452350651907,
                    0.16975979586073598843335396088583
                    ),
                _s21(0.040985629001117093392719996156396),
                _s21(0.0067312708878883654472690622726455),
                _s111(
                    0.093167409779881173436528549918676,
                    0.038318225821019367604389912994659
                    ),
                _s111(
                    0.39452027980019438558519271303259,
                    0.092656481520757486730406950326856
                    ),
                _s111(
                    0.16267741639447744436113261749936,
                    0.041188714248475323525533078597563
                    ),
                _s111(
                    0.25358901421887948684883593032112,
                    0.039570904970158003249540321664373
                    ),
                _s111(
                    0.36225224131779131247753949238273,
                    0.038592700174896100170570335774836
                    ),
                _s111(
                    0.28162257770616086095582688169845,
                    0.094534961736598957233915521157653
                    ),
                _s111(
                    0.38327266499265926647924168442454,
                    0.007387994632294168743732664906526
                    ),
                _s111(
                    0.27375035251626049961263706789435,
                    0.0075460031623127730020318602704793
                    ),
                _s21(0.49625527767573514001389657545647),
                _s111(
                    0.094121342797366048030078297547108,
                    0.0072345584577821073842913917176328
                    ),
                _s21(0.26423131543827256642185930872271),
                _s111(
                    0.18039615188676572950824015046539,
                    0.095566269527365202849336707332193
                    ),
                _s21(0.48061256179250326980246841259339),
                _s111(
                    0.17473734628280572255989948921258,
                    0.0079879218808479043737869646278553
                    ),
                _s21(0.096328495599215201173201535935881),
                _s111(
                    0.037291472051291242304216084831107,
                    0.0080749108702087576310869003329958
                    ),
                _s21(0.37535292670208629082280160094505),
                ])
            self.weights = numpy.concatenate([
                0.012545689845600317538365061859693 * numpy.ones(1),
                0.013110532701885239618009366987256 * numpy.ones(3),
                0.010379016056400192481851923942355 * numpy.ones(3),
                0.014145045806484844874765314309247 * numpy.ones(6),
                0.015274442601324624824852790264274 * numpy.ones(6),
                0.0038336997309291829397721127368518 * numpy.ones(3),
                0.00061725450549664320419212747548882 * numpy.ones(3),
                0.005366271454167767620352544823192 * numpy.ones(6),
                0.015031854349741338124379725976935 * numpy.ones(6),
                0.0072041341747974274505684332803814 * numpy.ones(6),
                0.0089048769281635644125166070048955 * numpy.ones(6),
                0.0099472518756824178795286531005395 * numpy.ones(6),
                0.014352351578157452140914597865316 * numpy.ones(6),
                0.0042421492668037687246552423274008 * numpy.ones(6),
                0.004081275077116450441131016718494 * numpy.ones(6),
                0.0043432467221706974057545473509849 * numpy.ones(3),
                0.0025892123823979845679883365231832 * numpy.ones(6),
                0.020520008671509842678908172943379 * numpy.ones(3),
                0.01184356214254310744956332163135 * numpy.ones(6),
                0.010352494770852601438941460699465 * numpy.ones(3),
                0.0037072267642463084761216462460982 * numpy.ones(6),
                0.010027393067388905272557992426155 * numpy.ones(3),
                0.0017969475854465761889633765854462 * numpy.ones(6),
                0.01899458651735265609462072817133 * numpy.ones(3),
                ])
        elif index == 25:
            bary = numpy.concatenate([
                _s21(0.38764203040456340905403082721048),
                _s21(0.21100450806149663445070386542963),
                _s111(
                    0.44041692747934335009593371624103,
                    0.0018188666342743112766696723682803
                    ),
                _s21(0.29949231580450848407125892240039),
                _s111(
                    0.15900790619732789070330324480018,
                    0.036960141579671466022337535655029
                    ),
                _s21(0.037222925992440749253778180782458),
                _s111(
                    0.17735379675725289478107376926366,
                    0.078858068005635243914956961270277
                    ),
                _s21(0.14510924357450031655763093628465),
                _s111(
                    0.27006673582095941033489077509943,
                    0.068847529431497900296745799538768
                    ),
                _s111(
                    0.34139103302114986743291331307247,
                    0.11599980764096015542626603561511
                    ),
                _s111(
                    0.37393797971958445012298596356906,
                    0.048317434287376906056298575156586
                    ),
                _s111(
                    0.099133063341682212333239560395697,
                    0.0071283145012573502514738370963329
                    ),
                _s111(
                    0.29950641862967453779250339452183,
                    0.2036929105842509641504530263135
                    ),
                _s21(0.42475930454057472879850455682139),
                _s111(
                    0.17862984860361624224222586297493,
                    0.0072361617479481281736797158331883
                    ),
                _s111(
                    0.36206880189597202166404995230001,
                    0.012913883250032515811837236638836
                    ),
                _s111(
                    0.088792915489366582036185478491327,
                    0.037687949784259056157934844905562
                    ),
                _s111(
                    0.23362281014171527436003646558132,
                    0.13700669408707093392067534012605
                    ),
                _s21(0.46220870874870612934990684931817),
                _s111(
                    0.25659540970901987231619859240365,
                    0.024540060247524317710385555000791
                    ),
                _s21(0.092949701700769890257855236623887),
                _s111(
                    0.041068819111784664276338540242807,
                    0.0071888282616930215474865464613266
                    ),
                _s21(0.007835344282603702842444951610214),
                _s21(0.4890393696603954899865458352552),
                _s111(
                    0.27941618864926070287458309718122,
                    0.00089146431749809121705039411647924
                    ),
                ])
            self.weights = numpy.concatenate([
                0.013689851548272245320429203276484 * numpy.ones(3),
                0.011587263236010591535133261201506 * numpy.ones(3),
                0.0016748178319347049952071121902043 * numpy.ones(6),
                0.01801764070170147315121202886233 * numpy.ones(3),
                0.0063114780247592747702748960908622 * numpy.ones(6),
                0.0033972977219047357678698140289552 * numpy.ones(3),
                0.0095150215674557713456831852235084 * numpy.ones(6),
                0.011491525862564797113493652010033 * numpy.ones(3),
                0.0108843936124369191658703405362 * numpy.ones(6),
                0.01584035228789843632915119649126 * numpy.ones(6),
                0.010640170695508785039236842558269 * numpy.ones(6),
                0.0025452716253490144217593778231769 * numpy.ones(6),
                0.017913820892276058503532532883596 * numpy.ones(6),
                0.015911310137458419121728140818281 * numpy.ones(3),
                0.0032637396820492422392932473232668 * numpy.ones(6),
                0.0054546383679744283888000521419753 * numpy.ones(6),
                0.0052725619214294191574035660969367 * numpy.ones(6),
                0.013740082592022550833626226582498 * numpy.ones(6),
                0.013654275187528014462542890132225 * numpy.ones(3),
                0.0073143409079328459530336271323825 * numpy.ones(6),
                0.0091828212598200358239500658732806 * numpy.ones(3),
                0.0016929836341273321280349280552852 * numpy.ones(6),
                0.00080651028832461670020354395896671 * numpy.ones(3),
                0.0084440859465210771191602047950713 * numpy.ones(3),
                0.0015117020784588803378981330586815 * numpy.ones(6),
                ])
        elif index == 26:
            bary = numpy.concatenate([
                _s111(
                    0.080071654940316590838670809367462,
                    0.0047946609754366080496785851291657
                    ),
                _s111(
                    0.031643611571530775880568694708158,
                    0.029155196206835808295179789331484
                    ),
                _s21(0.066737122576466204004912920243451),
                _s111(
                    0.075380047515398673674403319207208,
                    0.026209364022498616142328431058363
                    ),
                _s111(
                    0.033100034336032280699576132203046,
                    0.0056981179168751913974406376049972
                    ),
                _s21(0.0063401164920767662570293326034631),
                _s111(
                    0.13248618961456728907976402318946,
                    0.041724722742120915302855032431712
                    ),
                _s111(
                    0.10868713291440214809032135393378,
                    0.10004565910652748810723470715608
                    ),
                _s21(0.49375303289638477722397143285906),
                _s111(
                    0.25027231329052646476620244210613,
                    0.12061440220524896877101964029377
                    ),
                _s111(
                    0.38902206204276180468775743835644,
                    0.029537942516907761430810798332271
                    ),
                _s111(
                    0.35850929642766155649301394664776,
                    0.087378465163844456157009687830804
                    ),
                _s21(0.38878749710759403945319188490627),
                _s111(
                    0.18686917947622156628946218281807,
                    0.076311901512959348062553797782758
                    ),
                _s21(0.273147100929078727685059556205),
                _s111(
                    0.41470590959030631628242487475306,
                    0.002057530965370791445095127630898
                    ),
                _s111(
                    0.31941530538343875268651555354955,
                    0.17047872849724888036507699302919
                    ),
                _s111(
                    0.14373762619976407168790501086106,
                    0.0079996080914842375618315615950868
                    ),
                _s21(0.47182856332116602399781116458116),
                _s21(0.15420143036454427894089009106215),
                _s111(
                    0.28378813885947043930105587200173,
                    0.051165873685137733569942363619641
                    ),
                _s21(0.21204316330220561881487732137648),
                _s21(0.43598541938438321511425602403781),
                _s111(
                    0.21654666647347715986673330645389,
                    0.022784599250895603686015228976256
                    ),
                _s3(),
                _s111(
                    0.31289850307487997053616803429133,
                    0.0094732979122135237896856652699842
                    ),
                _s111(
                    0.22643479740771754107085813522835,
                    0.00046400773217562334693505498682472
                    ),
                ])
            self.weights = numpy.concatenate([
                0.0013985264481602720573708664332252 * numpy.ones(6),
                0.0012055647737168855619845177402375 * numpy.ones(6),
                0.004913825302966016981007142607323 * numpy.ones(3),
                0.0033055447129676699445092700984397 * numpy.ones(6),
                0.0010857073429967548036368156160328 * numpy.ones(6),
                0.00052695311668187181624262721967129 * numpy.ones(3),
                0.0064035978997128186079669204264305 * numpy.ones(6),
                0.0046142110763783177116764180157482 * numpy.ones(6),
                0.0053021591818673456358156535624305 * numpy.ones(3),
                0.014379473227598739028032494046501 * numpy.ones(6),
                0.0082597672170868395541711408505236 * numpy.ones(6),
                0.013727958216085702692626792164197 * numpy.ones(6),
                0.019468067837182877110795160590948 * numpy.ones(3),
                0.010397645528174323211119160889642 * numpy.ones(6),
                0.019535646923247540479483246571919 * numpy.ones(3),
                0.0018571474709980839928543177592289 * numpy.ones(6),
                0.017599167180695208001678034213183 * numpy.ones(6),
                0.0029667616626565056104558924660712 * numpy.ones(6),
                0.01152850363465689191468910016604 * numpy.ones(3),
                0.013255259448545270082976377764583 * numpy.ones(3),
                0.010107124432088683924275287655253 * numpy.ones(6),
                0.01694434507852809070881172141097 * numpy.ones(3),
                0.016412400602587902609478997606524 * numpy.ones(3),
                0.0062693378460805675891755560958426 * numpy.ones(6),
                0.020486662589223242648548008188463 * numpy.ones(1),
                0.0045915583873986369459740069540588 * numpy.ones(6),
                0.0011395489158682133180844934604399 * numpy.ones(6),
                ])
        elif index == 27:
            bary = numpy.concatenate([
                _s21(0.38071402118118718796829089353041),
                _s21(0.44666780370386460487281769646911),
                _s21(0.41614137880541216728641979593873),
                _s21(0.080304647788438414426593495392939),
                _s111(
                    0.28704219659349661588432818174586,
                    0.030730604727272810845711493729038
                    ),
                _s111(
                    0.34508784171556842030052014331628,
                    0.12915264006344966410950076008599
                    ),
                _s21(0.23340040666987110652849456037186),
                _s111(
                    0.37593015704866182275994503752674,
                    0.028033486095249930819305696289918
                    ),
                _s21(0.30116546516650913524224010144214),
                _s21(0.1747799663549000372916494747766),
                _s111(
                    0.31694558893313200191103822152509,
                    0.20913092113766868968994557526703
                    ),
                _s111(
                    0.40722839304271988463728864988376,
                    0.066038912849738622384089279411697
                    ),
                _s111(
                    0.21355359845782394893368907639196,
                    0.041030576819181776106109047276815
                    ),
                _s111(
                    0.32885287806889269108666611335766,
                    0.0052996403717989955028946951693256
                    ),
                _s21(0.48556505418516278581248566197757),
                _s111(
                    0.13929530614214873916039472808702,
                    0.063073995414950836834842369776647
                    ),
                _s21(0.032571520180181552088550070594271),
                _s111(
                    0.25524625469697804130922508895929,
                    0.14896285093824008952482094132574
                    ),
                _s111(
                    0.20837601560037408588971921829953,
                    0.094697082433130656936636331033664
                    ),
                _s111(
                    0.44001055194621543988288199437566,
                    0.0055807170152600880295439706314057
                    ),
                _s21(0.12757090190467754450607207557053),
                _s111(
                    0.30222094122782113648292368051704,
                    0.075076902433196199210983219190141
                    ),
                _s111(
                    0.081946802583533692436184081313867,
                    0.006982529324458954022873857211846
                    ),
                _s111(
                    0.034364969912142005015294073277258,
                    0.0060935694037647949524948145245339
                    ),
                _s111(
                    0.080112073847101132529494201513186,
                    0.035034422527697321017424139691128
                    ),
                _s21(0.0066392191809587244576727440727105),
                _s111(
                    0.14721343189892252152803619049473,
                    0.01935200131803894324014422917084
                    ),
                _s111(
                    0.22971965325784322947520568074313,
                    0.0073324725490404178237779817972998
                    ),
                _s111(
                    0.14765552111986981014776075370802,
                    0.00049032844346289946370399533189533
                    ),
                ])
            self.weights = numpy.concatenate([
                0.0095600849674599199772576169200097 * numpy.ones(3),
                0.009410159809454225705038652025013 * numpy.ones(3),
                0.012050227024150432419406686253757 * numpy.ones(3),
                0.005212621872801876185064498019083 * numpy.ones(3),
                0.0055317948337667311730831625299371 * numpy.ones(6),
                0.012557436204036536566796674175289 * numpy.ones(6),
                0.01347131539804937431696711666857 * numpy.ones(3),
                0.006395152699454402968162347285102 * numpy.ones(6),
                0.015747965781362654080346754630112 * numpy.ones(3),
                0.011282442544698378679654541532637 * numpy.ones(3),
                0.013715393230550838092393500376707 * numpy.ones(6),
                0.0098622701189895686063015375505419 * numpy.ones(6),
                0.0064553729049296897504941446202009 * numpy.ones(6),
                0.0029278263617991022689234036214218 * numpy.ones(6),
                0.0071172374128746417554186634291104 * numpy.ones(3),
                0.0071306353104870229945258138606624 * numpy.ones(6),
                0.0027773395289541807812794466476234 * numpy.ones(3),
                0.012347663130861363158087059453057 * numpy.ones(6),
                0.01069370058961626192472955794923 * numpy.ones(6),
                0.0032424675976393404835269067602781 * numpy.ones(6),
                0.0097432449228177301892647340422733 * numpy.ones(3),
                0.010930611092913285709683376898382 * numpy.ones(6),
                0.0020151231272897023219345279609102 * numpy.ones(6),
                0.0011967736084731628118931094838474 * numpy.ones(6),
                0.0043271580353607124803304518608826 * numpy.ones(6),
                0.00057544240567050243730277095982355 * numpy.ones(3),
                0.0046223871117811095525012414534135 * numpy.ones(6),
                0.0033942537388070269670006887510759 * numpy.ones(6),
                0.00084660613576385057279842151172511 * numpy.ones(6),
                ])
        elif index == 28:
            bary = numpy.concatenate([
                _s21(0.30398292251648411800214900044008),
                _s111(
                    0.045505400558346431727807991781216,
                    0.021524385369456063093113322309787
                    ),
                _s111(
                    0.21339445476708733127892273592128,
                    0.049069669357559449497290014316515
                    ),
                _s21(0.0048041261966580130434937831717468),
                _s21(0.45827990424041196248566014616575),
                _s111(
                    0.24210251191931968798493919542294,
                    0.17765845029637026549346157665718
                    ),
                _s21(0.38626797357004206715306041066432),
                _s111(
                    0.32719073201917001654062780931358,
                    0.18981235629273680638945450911907
                    ),
                _s111(
                    0.14199816693317424158549365870721,
                    0.0044583820232892991106870594130978
                    ),
                _s111(
                    0.17539639319146174366217488397583,
                    0.087677976484352009621555620911487
                    ),
                _s111(
                    0.39213961333441454326181776757535,
                    0.063180327634410631288716858638044
                    ),
                _s21(0.25826407215046216036176320139954),
                _s111(
                    0.33414561503592132590885785374857,
                    0.0041494641339236350830453966612408
                    ),
                _s111(
                    0.17414619605118213887800128063667,
                    0.022794804925916208200139716412466
                    ),
                _s111(
                    0.27583900807182419432492586989451,
                    0.022700844371796940319777755342365
                    ),
                _s111(
                    0.026222667164652288463562817190334,
                    0.0061496485426639406096403132305216
                    ),
                _s111(
                    0.24304720236592617066340524620428,
                    0.11203362934227092558880481787411
                    ),
                _s111(
                    0.22999298405790715858777861386225,
                    0.0047814897729871129243110975886578
                    ),
                _s21(0.10589584417862763192509688251439),
                _s21(0.42955220211889932894027225677447),
                _s111(
                    0.29560828087240168665661714946321,
                    0.062448742179632832778367362185847
                    ),
                _s111(
                    0.12139345075409118696645762862512,
                    0.050211185913428064801792249229359
                    ),
                _s111(
                    0.37382380031020968435675418592448,
                    0.025727998742878701946738615612541
                    ),
                _s111(
                    0.44278340652024356022983523781947,
                    0.0056465659934661067092750280877533
                    ),
                _s111(
                    0.33779866320058232199335327638265,
                    0.11808906971509501899057676106634
                    ),
                _s111(
                    0.092229189195282251917034530595928,
                    0.018242291012294672032147019482753
                    ),
                _s21(0.48484113256258945491645544181893),
                _s21(0.15863768886305965531677723825212),
                _s21(0.060839192392758723217217197090234),
                _s111(
                    0.070290740478132760242433078634059,
                    0.0012002556014871294389479856551818
                    ),
                ])
            self.weights = numpy.concatenate([
                0.014362466300646133348027811264005 * numpy.ones(3),
                0.0021175395576808984412355289862133 * numpy.ones(6),
                0.0058956722345142461623934541375487 * numpy.ones(6),
                0.00031113520868149634857368224076693 * numpy.ones(3),
                0.0088517050108931599307105507663108 * numpy.ones(3),
                0.012557256216676878683396712136353 * numpy.ones(6),
                0.014210586390448186336263636315239 * numpy.ones(3),
                0.01310611412409938527984042280243 * numpy.ones(6),
                0.0017788954192555133655591802186126 * numpy.ones(6),
                0.0080119292866949631249745711753644 * numpy.ones(6),
                0.0084646957342766040111574135901113 * numpy.ones(6),
                0.013092748589088967023817607375257 * numpy.ones(3),
                0.0023767642676877834854026483470515 * numpy.ones(6),
                0.004276942335017944414567123499132 * numpy.ones(6),
                0.0051312773820372967076614848891642 * numpy.ones(6),
                0.0009485404258894004150457932126171 * numpy.ones(6),
                0.010245986561704200664923983034398 * numpy.ones(6),
                0.0023245453644640500783644556963718 * numpy.ones(6),
                0.0078099433710863412309597005727567 * numpy.ones(3),
                0.013547595603325509375043201955697 * numpy.ones(3),
                0.0089251161402765502463289795801276 * numpy.ones(6),
                0.0059178837723538902653686919548426 * numpy.ones(6),
                0.0062605346429685520732367493629901 * numpy.ones(6),
                0.0031401776434880361487834620924543 * numpy.ones(6),
                0.012787138422955800232021009842014 * numpy.ones(6),
                0.0032049161907930344576513517309852 * numpy.ones(6),
                0.0076821105785950227360524351900793 * numpy.ones(3),
                0.011942669640248226695693570371398 * numpy.ones(3),
                0.0051282520046780679426911116937038 * numpy.ones(3),
                0.00072513459498608292483699650528074 * numpy.ones(6),
                ])
        elif index == 29:
            bary = numpy.concatenate([
                _s111(
                    0.058942108840229254909367118374802,
                    0.0027287432479209962438792745425867
                    ),
                _s111(
                    0.34978801000933185820097161634175,
                    0.15717769986719341947263402850903
                    ),
                _s111(
                    0.32300182354355017923536854359506,
                    0.0021009666448274998140089350274275
                    ),
                _s21(0.49891482463768621240406300243043),
                _s111(
                    0.15814585424951614707804832797578,
                    0.068165808813746372963658642301914
                    ),
                _s111(
                    0.029549468261353386167561625955184,
                    0.010830958603609277837545055226642
                    ),
                _s111(
                    0.29181917342653711427349171295998,
                    0.21893234198017248474276829577936
                    ),
                _s111(
                    0.075522178512995823960412864572712,
                    0.021286896240733180131524820543607
                    ),
                _s111(
                    0.11711666950889887867191977240736,
                    0.040847216576102416159398844906716
                    ),
                _s111(
                    0.011166218108169205892804217668527,
                    0.0016034964960437314456385894414862
                    ),
                _s21(0.43438042676173062781327364063836),
                _s21(0.041097335627118032252800529751715),
                _s111(
                    0.20804564927908717362195007963557,
                    0.10154598522683395980648078606482
                    ),
                _s111(
                    0.39221011498043485569761702233906,
                    0.041527061268822623183170728076752
                    ),
                _s111(
                    0.3597112755099759108283990366837,
                    0.09384904114513237175837095732163
                    ),
                _s111(
                    0.2458746994828754694858404833332,
                    0.0086860298043840857771487176445225
                    ),
                _s111(
                    0.16700273817492316672365811386467,
                    0.017589124044045597694563605245405
                    ),
                _s21(0.20840530513240084162740931094792),
                _s111(
                    0.11500859863194643021826019106213,
                    0.005523524512212532521621970687098
                    ),
                _s21(0.16074588443196363301272685071313),
                _s111(
                    0.31539539811731915817097947090314,
                    0.023858926942655562219750175473007
                    ),
                _s111(
                    0.2232226502248206089186413426559,
                    0.040295334544771788772348237499174
                    ),
                _s111(
                    0.28839585991873242258544416691439,
                    0.067878404311447032986904882533314
                    ),
                _s21(0.48840160293260276689430859540177),
                _s111(
                    0.26736635027277558219825934618575,
                    0.13953560718108261149967220769318
                    ),
                _s21(0.30238641121512847136905799060117),
                _s111(
                    0.40576539529889156272216195395399,
                    0.0080665857041665784150860629765953
                    ),
                _s21(0.11442681299442553798565205312521),
                _s21(0.46476243108073894476491473417097),
                _s21(0.073721881390099808650445025675325),
                _s21(0.39061917878326375927003022057437),
                _s111(
                    0.18632072767535956156184471835257,
                    0.00012344681228738954130308260746232
                    ),
                ])
            self.weights = numpy.concatenate([
                0.00076925297147624881541618592874533 * numpy.ones(6),
                0.010452214696922434116475781603701 * numpy.ones(6),
                0.0013528239492524084865187580268467 * numpy.ones(6),
                0.0015164621031569022804849652219559 * numpy.ones(3),
                0.0064391367075338057816430532666317 * numpy.ones(6),
                0.0012721944371716311876992357567245 * numpy.ones(6),
                0.0135452465720075725102823727014 * numpy.ones(6),
                0.0027615307495904436564342003737674 * numpy.ones(6),
                0.0045126868454873290954419124247585 * numpy.ones(6),
                0.00026949223185870118305706579242832 * numpy.ones(6),
                0.011171100295167859073077023799404 * numpy.ones(3),
                0.0028604915061568867333824976118823 * numpy.ones(3),
                0.0093863306762836635630770545413452 * numpy.ones(6),
                0.0078229184583064779212229331699618 * numpy.ones(6),
                0.010630396363196612234552811744833 * numpy.ones(6),
                0.0030512388233451906013434613602514 * numpy.ones(6),
                0.0038253666717300118697667524740192 * numpy.ones(6),
                0.012539203442623228232605538347198 * numpy.ones(3),
                0.0017416952313017492927569814888272 * numpy.ones(6),
                0.010571417858227575428688503333985 * numpy.ones(3),
                0.0056691534816654586699581575228516 * numpy.ones(6),
                0.0065819977638522062148492158727196 * numpy.ones(6),
                0.009178924164927596766914931789319 * numpy.ones(6),
                0.0061344011647189098349072547281056 * numpy.ones(3),
                0.012518366498834426412812610334829 * numpy.ones(6),
                0.016310238807743203463864623783307 * numpy.ones(3),
                0.0036413404112807363179109181584966 * numpy.ones(6),
                0.0081728784412271315068178679342896 * numpy.ones(3),
                0.010313116635258498387921242243655 * numpy.ones(3),
                0.0056088471328313043216251727623595 * numpy.ones(3),
                0.015913215284088902811957883111148 * numpy.ones(3),
                0.00068867262504176093086598589556804 * numpy.ones(6),
                ])
        elif index == 30:
            bary = numpy.concatenate([
                _s111(
                    0.25905388452106734346190162559017,
                    0.047835123140772525935707948872268
                    ),
                _s21(0.0033187249366445024546421686949577),
                _s111(
                    0.39157218829125636736671675373777,
                    0.079659526931600594756749475577726
                    ),
                _s111(
                    0.3561402832962335761867221815975,
                    0.057693401273874186082249665297055
                    ),
                _s21(0.072372407224677937363653532290462),
                _s111(
                    0.28301249734958874321444197996882,
                    0.07726143757688407084955123286177
                    ),
                _s21(0.047157910242171889312037722954373),
                _s111(
                    0.24161376245152440520999753571001,
                    0.022758384295000031477474878057185
                    ),
                _s111(
                    0.25278124718793294503451217101498,
                    0.12381197877067458674981066703959
                    ),
                _s111(
                    0.18490610638391715516325237855078,
                    0.1158819672361005416523880538428
                    ),
                _s21(0.46803017365112542953577929820856),
                _s21(0.012686604674467735006027071136378),
                _s111(
                    0.19370437715364017299668765882381,
                    0.066517444781881577464741673473775
                    ),
                _s111(
                    0.076401248439387560468254911909187,
                    0.0044387813770613007568055630602655
                    ),
                _s21(0.1215915082227280044392545172266),
                _s111(
                    0.209048087452689630604774294769,
                    0.0046635793926885863685701332403891
                    ),
                _s111(
                    0.29854299405924263801892554849117,
                    0.0047036817644770118850262481870566
                    ),
                _s111(
                    0.33437904003403105409278809146273,
                    0.025182066703868666108167264025512
                    ),
                _s111(
                    0.12323080238069512477715681793169,
                    0.065776573824742854772656838431388
                    ),
                _s111(
                    0.33851416242984571855093801389233,
                    0.12612409498498941941507065213661
                    ),
                _s111(
                    0.35440360218068748258720004415744,
                    0.19487095092351842924645557434142
                    ),
                _s111(
                    0.26306829775780833706284381518364,
                    0.19100142457228306593144803193292
                    ),
                _s21(0.18240956151745298953057022289716),
                _s111(
                    0.43456617396466968128174894076231,
                    0.02753340612454985204798104771481
                    ),
                _s111(
                    0.16407098706987833574252433086152,
                    0.028063921981372955201358243854412
                    ),
                _s21(0.36228727935294983052988138602057),
                _s111(
                    0.042681999706080433705516301689041,
                    0.015902416268934672464118566688798
                    ),
                _s111(
                    0.093959794652729919461913932619057,
                    0.027294230652095704742840503081656
                    ),
                _s111(
                    0.13540885351993447169036013117013,
                    0.0056912114454160885997692932502306
                    ),
                _s21(0.43674324854846020764151179888262),
                _s21(0.27242804078392824107538918178154),
                _s111(
                    0.39622151473965911642760884527305,
                    0.005162347016621281058036811942699
                    ),
                _s21(0.49731933900030860282816091251702),
                _s111(
                    0.029484042597673956354915562753704,
                    0.00053370866069444535749021356719296
                    ),
                ])
            self.weights = numpy.concatenate([
                0.004214758463912433883152453188975 * numpy.ones(6),
                0.00017172990137104942468678980154104 * numpy.ones(3),
                0.0059249227450920147883356704504482 * numpy.ones(6),
                0.0059442201844244871419505014864301 * numpy.ones(6),
                0.0038019326684155175153829396834688 * numpy.ones(3),
                0.0068771843875359256357175391280748 * numpy.ones(6),
                0.0028444336676874951028018584478761 * numpy.ones(3),
                0.0040924989683471987760409642912016 * numpy.ones(6),
                0.0088726545060173444449149483165618 * numpy.ones(6),
                0.0075343229295464820076446697134052 * numpy.ones(6),
                0.0077841287323308469014853338770732 * numpy.ones(3),
                0.00084655014532888889814677262744939 * numpy.ones(3),
                0.0068040067561018738128716234210671 * numpy.ones(6),
                0.0012066965685421181155010859501611 * numpy.ones(6),
                0.0073869118349142885461174099139051 * numpy.ones(3),
                0.0019589485778932435157537933732866 * numpy.ones(6),
                0.0022744459009986780635765612529568 * numpy.ones(6),
                0.0053646848461863093012265838596212 * numpy.ones(6),
                0.0058700698030655262212872321479898 * numpy.ones(6),
                0.011303023542937485123353969678306 * numpy.ones(6),
                0.014316485269667588811858161811095 * numpy.ones(6),
                0.012926023450118294440080637252012 * numpy.ones(6),
                0.010752850965432809786603660918464 * numpy.ones(3),
                0.0062815965808916628801687659682508 * numpy.ones(6),
                0.0047019369941059667975799444504668 * numpy.ones(6),
                0.015596091325825796361722886942984 * numpy.ones(3),
                0.0018305068761488277462659432624979 * numpy.ones(6),
                0.00382188054709502025846421651054 * numpy.ones(6),
                0.0019198805574689176158305004731994 * numpy.ones(6),
                0.012404959028146008346000308716203 * numpy.ones(3),
                0.014914550762685758374897773272937 * numpy.ones(3),
                0.0026425679231632510990413400357792 * numpy.ones(6),
                0.0027913181197407686329486657759208 * numpy.ones(3),
                0.00033562171146640224065236065543195 * numpy.ones(6),
                ])
        elif index == 31:
            bary = numpy.concatenate([
                _s3(),
                _s21(0.48615246791305980749502825684294),
                _s111(
                    0.13909062847209461455690649557882,
                    0.10666968086407751050505742244438
                    ),
                _s21(0.41958808086426048984931948341944),
                _s111(
                    0.19723313078965306817875179488104,
                    0.0032611556637282505123819284932983
                    ),
                _s111(
                    0.23852415098431002503444667826602,
                    0.013139817840403620919779773530426
                    ),
                _s21(0.079206796506698645657074612240253),
                _s111(
                    0.02975072938598998751985786658237,
                    0.0023442211286868733673772054963996
                    ),
                _s21(0.017731349249714062780956468973933),
                _s111(
                    0.39578393329688537518355800879966,
                    0.0019614003698949870869893553009046
                    ),
                _s111(
                    0.35168005877608098360050792170723,
                    0.1526365657752196865156797968973
                    ),
                _s111(
                    0.20248219426175837173512130799879,
                    0.11342527834750284896977187892534
                    ),
                _s111(
                    0.42255395522221506716621355889,
                    0.054521384444683095426382848439811
                    ),
                _s111(
                    0.12468050884744401666377169353448,
                    0.060604929160772151103121935366406
                    ),
                _s111(
                    0.25554623640285731183017130349463,
                    0.037588952157517732370451852215149
                    ),
                _s111(
                    0.074503129935386320724306355593389,
                    0.0043863301700615852888398608689818
                    ),
                _s111(
                    0.41680727552365056682725655812886,
                    0.017273414324937101756342976184492
                    ),
                _s111(
                    0.19287172300619847291687804013997,
                    0.062578022014310871920833350966997
                    ),
                _s111(
                    0.047347818939552459451385294349203,
                    0.018407947031282917275717509441113
                    ),
                _s21(0.21843649679092809110388908679191),
                _s21(0.17245078718388606926746338618082),
                _s21(0.051095180822381385186897550490234),
                _s111(
                    0.34190500534752768934662834879714,
                    0.042903208942771387190820243344231
                    ),
                _s111(
                    0.32581395305750971778450701661229,
                    0.013189083779015684698044587029263
                    ),
                _s111(
                    0.28942594875713730987036549562486,
                    0.0011968816307095387402021589123163
                    ),
                _s111(
                    0.2987897956760133454300019671621,
                    0.21913276143017678236160397647493
                    ),
                _s111(
                    0.27447574774279502716293320084683,
                    0.081129891569330921100470047347482
                    ),
                _s111(
                    0.35322987389141561914984619459984,
                    0.095278114482016779136810285591491
                    ),
                _s21(0.49774862062647029311275000271722),
                _s111(
                    0.26636788247507180299581537514585,
                    0.1468887080021118560295241568475
                    ),
                _s21(0.38516846924568008554517859744433),
                _s21(0.29749149036799287470764273818236),
                _s111(
                    0.09604815210105114997668592653903,
                    0.027530630078587428956497597674641
                    ),
                _s21(0.0050828094627278708316254269460485),
                _s111(
                    0.16740782477835609118527352318765,
                    0.025875108219938122666084139460055
                    ),
                _s21(0.45034094173901821314635513779581),
                _s111(
                    0.13067763377028836277467568649085,
                    0.0057973186847397622504743490255103
                    ),
                ])
            self.weights = numpy.concatenate([
                0.0083431208138863716451842500150846 * numpy.ones(1),
                0.0040439883965193655109454020800179 * numpy.ones(3),
                0.0058062932483636674569302399220189 * numpy.ones(6),
                0.0088430826430883116704786445077798 * numpy.ones(3),
                0.0012838516987696961172806936114204 * numpy.ones(6),
                0.0028265142759277386727839666817174 * numpy.ones(6),
                0.0042237298511032305918995026024181 * numpy.ones(3),
                0.00050550904687225829324104734468596 * numpy.ones(6),
                0.00098686614502708321176285388431917 * numpy.ones(3),
                0.0013026954489734516227727429910979 * numpy.ones(6),
                0.010192579444341274262470441326086 * numpy.ones(6),
                0.007866870495738138652984691691726 * numpy.ones(6),
                0.0071208092414187916689902913930816 * numpy.ones(6),
                0.0049641144765076699299059654354795 * numpy.ones(6),
                0.0054364374730135536327491726066986 * numpy.ones(6),
                0.0011525829400842612783173807408691 * numpy.ones(6),
                0.0041783326918215304786407671259906 * numpy.ones(6),
                0.0063624230528668775786528758764828 * numpy.ones(6),
                0.0018127692882839107284968660329159 * numpy.ones(6),
                0.011237003508013174398521468050714 * numpy.ones(3),
                0.0095591959667284181398063232977943 * numpy.ones(3),
                0.0032590866501857710047571063808041 * numpy.ones(3),
                0.0066899592450456260425207567082523 * numpy.ones(6),
                0.0036397134849474400744209001676392 * numpy.ones(6),
                0.0009293703376931178577027761151765 * numpy.ones(6),
                0.012895586787984193472212205060013 * numpy.ones(6),
                0.0084758977171366182863034974445765 * numpy.ones(6),
                0.0096222755316967225902273459809358 * numpy.ones(6),
                0.0023276606420564320665204344760308 * numpy.ones(3),
                0.011599669031385622625288510119146 * numpy.ones(6),
                0.013739395196568263425534557699832 * numpy.ones(3),
                0.012359040509346901920745207354595 * numpy.ones(3),
                0.0036665867049268454154813460041555 * numpy.ones(6),
                0.00034440137927258108939381879431748 * numpy.ones(3),
                0.004445337991252401270686286100512 * numpy.ones(6),
                0.010433977295780582619573480449208 * numpy.ones(3),
                0.0018212527841224722251057920612272 * numpy.ones(6),
                ])
        elif index == 32:
            bary = numpy.concatenate([
                _s21(0.0014574454623699960928662778842528),
                _s111(
                    0.20683168631590207728382782302181,
                    0.15886794122904668335167432996975
                    ),
                _s111(
                    0.44676142653316705754340224483135,
                    0.014766442430854609038809583622015
                    ),
                _s21(0.49865089357884733467317138263189),
                _s111(
                    0.23246620764637869388233078609472,
                    0.051431127240233004536278838047479
                    ),
                _s111(
                    0.11683492670695798084402968185587,
                    0.019559740364592332299633003710662
                    ),
                _s111(
                    0.32824054584493401035362710309253,
                    0.012618566750157142476855595349733
                    ),
                _s111(
                    0.25935093200312021210061404563916,
                    0.021062156675639160581433452598795
                    ),
                _s111(
                    0.1654639982373495864465131078538,
                    0.1175435092543479168392925250652
                    ),
                _s111(
                    0.38246490738093727219101417291267,
                    0.029298850540485161765447511368629
                    ),
                _s111(
                    0.39781707264672087785504525338517,
                    0.003442821652274290588127036401571
                    ),
                _s21(0.48116811599583680465385465256298),
                _s111(
                    0.32541746964669014280029535198219,
                    0.090370281144864205393606641772195
                    ),
                _s21(0.45228575651160052699914299686186),
                _s111(
                    0.11333919148760312305027005272606,
                    0.050697635926740433626664113965624
                    ),
                _s111(
                    0.3110416328276417502244950214866,
                    0.047944478915561739305682504853682
                    ),
                _s111(
                    0.29883713538891299464399260564884,
                    0.0012766573802798304189189503653161
                    ),
                _s21(0.37632875804791900182308660196979),
                _s111(
                    0.15211228094378378643002415263567,
                    0.0042625634151115725540806091900697
                    ),
                _s111(
                    0.088661705033911549652583726772053,
                    0.0038525119887189992001816489583852
                    ),
                _s21(0.28301827984957322041438208655897),
                _s111(
                    0.31915258059306055052682569548033,
                    0.20190014863247661797831704906531
                    ),
                _s21(0.41367978432541339612384164928171),
                _s111(
                    0.16675569358213337607109843959011,
                    0.067710393221269095096865528910129
                    ),
                _s21(0.030269191040444305158298672945288),
                _s21(0.012810196480111965778556034789047),
                _s111(
                    0.24213792401943282841575674137468,
                    0.095572417110654464451918535782973
                    ),
                _s111(
                    0.27986174799782974750817077163953,
                    0.15037086457345822193050982534828
                    ),
                _s111(
                    0.067398400207256697201596044333727,
                    0.024349850684587075510488811093475
                    ),
                _s111(
                    0.37168796669342035410688335776626,
                    0.12938786736692337429221458160532
                    ),
                _s111(
                    0.40562814906588287128473423883937,
                    0.063227420456501166534839267828476
                    ),
                _s3(),
                _s21(0.062010649826264363373062757581289),
                _s111(
                    0.22198016573748106763882135701731,
                    0.0052666330207657266311770024644716
                    ),
                _s111(
                    0.17832568203786295913720771188737,
                    0.026158881979650974875697343558692
                    ),
                _s111(
                    0.042622190591409538319455252635676,
                    0.0055547678704322100910016934601372
                    ),
                _s21(0.229277996928070005670179524396),
                _s21(0.10195882690026968516835109295075),
                _s111(
                    0.01434430037928038266539914038517,
                    0.00062309802034501055088033433870033
                    ),
                ])
            self.weights = numpy.concatenate([
                0.000053323902333416106400254541987731 * numpy.ones(3),
                0.0083492191295379379607496132394322 * numpy.ones(6),
                0.0030406522577729978791953972317799 * numpy.ones(6),
                0.001414627247662839331887829984057 * numpy.ones(3),
                0.0052860436025324795676418401094526 * numpy.ones(6),
                0.0024932699246690685111871599878089 * numpy.ones(6),
                0.003020581056713626437496351816275 * numpy.ones(6),
                0.003644175576876181494317376340312 * numpy.ones(6),
                0.0074810896991601560555385628545728 * numpy.ones(6),
                0.004568154135572013594483137848108 * numpy.ones(6),
                0.0016852036405342876762487466186547 * numpy.ones(6),
                0.0055800151424132688782584690974637 * numpy.ones(3),
                0.0082092066464793067843961413896148 * numpy.ones(6),
                0.0083785907041262422595237160067135 * numpy.ones(3),
                0.0044656768476199326771769516189161 * numpy.ones(6),
                0.0062028347215228412773041773070433 * numpy.ones(6),
                0.00088889254315192676504862454858492 * numpy.ones(6),
                0.01266129456048602771725361565265 * numpy.ones(3),
                0.0014852293291851673203994954890674 * numpy.ones(6),
                0.0011127063455906364256688083694641 * numpy.ones(6),
                0.012882733177175360884631111259188 * numpy.ones(3),
                0.011973359102482773113618885833645 * numpy.ones(6),
                0.010880466026934249125420923745024 * numpy.ones(3),
                0.0058823160628198555212815224984492 * numpy.ones(6),
                0.0018592782531938896857646658701338 * numpy.ones(3),
                0.00080174105167307369825768250562602 * numpy.ones(3),
                0.0084487434719511223349188988017484 * numpy.ones(6),
                0.010211455916366077598470389560429 * numpy.ones(6),
                0.0024899591532165259138039824114071 * numpy.ones(6),
                0.0098885056550222806063784672836171 * numpy.ones(6),
                0.0074091131429364931200787471935183 * numpy.ones(6),
                0.013030070995762573565990654144053 * numpy.ones(1),
                0.003976348594925360105011977111671 * numpy.ones(3),
                0.0019455932409830669760362284180765 * numpy.ones(6),
                0.0042608365965001633997223569106397 * numpy.ones(6),
                0.0010553743696176169018919035445515 * numpy.ones(6),
                0.012535597725648880339379126948595 * numpy.ones(3),
                0.0066293991150591530370229835984074 * numpy.ones(3),
                0.00017008824774248790820827892340005 * numpy.ones(6),
                ])
        elif index == 33:
            bary = numpy.concatenate([
                _s21(0.086684239934569733403885167143719),
                _s21(0.40942796080264335789831473750982),
                _s111(
                    0.25815292878950676607686672575016,
                    0.086928635963286975610775750266817
                    ),
                _s111(
                    0.28496240197208291608117879251431,
                    0.11409203532283129909402218494055
                    ),
                _s21(0.37833353829624275368490350813993),
                _s21(0.49897971786005120720692194814208),
                _s111(
                    0.13406763384193414260892149363411,
                    0.091665331440051169079133873264866
                    ),
                _s111(
                    0.32691375020967157776175511558191,
                    0.21593939124341559037588409534933
                    ),
                _s21(0.25966066800872120209981351971012),
                _s111(
                    0.10084571598277561827694726457751,
                    0.022390316286583462861847746555183
                    ),
                _s111(
                    0.35523282540476398732937672183118,
                    0.099241937888028235915035812500651
                    ),
                _s111(
                    0.15840273419572484405682982860673,
                    0.020443675144174764098590795531386
                    ),
                _s21(0.47851197325638217631900795822411),
                _s111(
                    0.4012397484580546941894954458687,
                    0.058612118456650833465282430146349
                    ),
                _s111(
                    0.22732201316496659764736410630696,
                    0.01942491016557721264039335602839
                    ),
                _s111(
                    0.10400400013142704642525252037482,
                    0.050870192193764943302235508354503
                    ),
                _s111(
                    0.25161901936303670587634136664988,
                    0.0038062521288190870719914971813455
                    ),
                _s111(
                    0.19467152737204624402976673304363,
                    0.09395371036958993507480321518488
                    ),
                _s111(
                    0.17563963145803316213221778975729,
                    0.0036669776240120998599947063245873
                    ),
                _s111(
                    0.16542250395062627567949654189123,
                    0.050425703506065786761849786581334
                    ),
                _s111(
                    0.060482924090960285482021334627394,
                    0.0041121482809411533620253142203631
                    ),
                _s111(
                    0.33210356535752328062830463297139,
                    0.0043366708649928965603782236806839
                    ),
                _s111(
                    0.37941321684563875501696845549601,
                    0.14376417624054298010939641108396
                    ),
                _s21(0.023517443057913169160168272626804),
                _s111(
                    0.056568897248945651690568306749463,
                    0.021865818028370799119005418153264
                    ),
                _s21(0.30602386176036062032313753985622),
                _s21(0.14648905071664279703668095863974),
                _s21(0.054379203279496399107536698264657),
                _s111(
                    0.29098491626933011621441428450629,
                    0.16830733134250228878229520285562
                    ),
                _s111(
                    0.23684539448315156050309426989228,
                    0.048513979905208350923303097084203
                    ),
                _s111(
                    0.41454915297644301732280152430399,
                    0.0047372398809321298392971179394802
                    ),
                _s111(
                    0.30753618824094102801468181463031,
                    0.022135213907368180883525139396134
                    ),
                _s111(
                    0.11092666443347038081700973712993,
                    0.004393395109466625111034660994398
                    ),
                _s21(0.21819201586325396303785403094487),
                _s111(
                    0.31886312205636779071710155108764,
                    0.055160205356247328190026963273698
                    ),
                _s21(0.49199134445655327342218575518451),
                _s111(
                    0.3953677591013827776176428211142,
                    0.024399888469638662591600320838045
                    ),
                _s21(0.45339996161780528522617012208605),
                _s111(
                    0.21375783029735308284357643858972,
                    0.15063494182918608390955994202054
                    ),
                _s21(0.0047266605105681039478873612905816),
                _s111(
                    0.024779370217540628313928951140415,
                    0.0045023852882973229147156567184492
                    ),
                ])
            self.weights = numpy.concatenate([
                0.0032023640124142048739750967657622 * numpy.ones(3),
                0.0068346799344680854155768488467014 * numpy.ones(3),
                0.00501082345657489411288971055546 * numpy.ones(6),
                0.006266671514597953997320729841746 * numpy.ones(6),
                0.0089341618081367398654579625187736 * numpy.ones(3),
                0.0010410912372864206795195918106253 * numpy.ones(3),
                0.0051060923018633094664616213804899 * numpy.ones(6),
                0.009639246858840778084399681323672 * numpy.ones(6),
                0.0097388572366853671978795387616727 * numpy.ones(3),
                0.0024212380388684937883998177559395 * numpy.ones(6),
                0.0076973467270549324479462656148706 * numpy.ones(6),
                0.0029882446919490308137500787050675 * numpy.ones(6),
                0.0054643256840460344145840173631233 * numpy.ones(3),
                0.0061988540259433594751495123266679 * numpy.ones(6),
                0.0033169972939130725047893557203419 * numpy.ones(6),
                0.0037089296125892631657349261227552 * numpy.ones(6),
                0.0015426465184705627653738188610475 * numpy.ones(6),
                0.0066121724110684191526364376318673 * numpy.ones(6),
                0.0013551781957506493853195466887438 * numpy.ones(6),
                0.0049036283435529124798156335095857 * numpy.ones(6),
                0.00091576615718100729388391935449569 * numpy.ones(6),
                0.0018218162644265219030878226402802 * numpy.ones(6),
                0.0092423799017096075609741779767361 * numpy.ones(6),
                0.0013813690630279039675038718064968 * numpy.ones(3),
                0.001960631233693213500026487839849 * numpy.ones(6),
                0.012088739980916558016401500432842 * numpy.ones(3),
                0.0073774464967341880963588518235881 * numpy.ones(3),
                0.0032093504439416188195275454512005 * numpy.ones(3),
                0.010209917786396493163804847710144 * numpy.ones(6),
                0.00546377043737215336622131435287 * numpy.ones(6),
                0.0020111020097847252089391984830089 * numpy.ones(6),
                0.0042683347026344800840863270200627 * numpy.ones(6),
                0.0012878764088953661984957990951069 * numpy.ones(6),
                0.0099263761795441988993506939476073 * numpy.ones(3),
                0.0065998833272567722390003802831257 * numpy.ones(6),
                0.0038623620299660401788514807724095 * numpy.ones(3),
                0.004947432808808186575044622741984 * numpy.ones(6),
                0.0093337977696465760168257332065444 * numpy.ones(3),
                0.009181481557927099542980980107203 * numpy.ones(6),
                0.00029348273357801739472278653752721 * numpy.ones(3),
                0.00064400177434743147186589300111156 * numpy.ones(6),
                ])
        elif index == 34:
            bary = numpy.concatenate([
                _s111(
                    0.38512625183828101969680419254726,
                    0.021096760270101908015060797424162
                    ),
                _s111(
                    0.052820067345498798714885078556375,
                    0.0033881121503899924277959039102701
                    ),
                _s21(0.045982595154127250025235120269461),
                _s111(
                    0.050898765283131619036903666047845,
                    0.018208328799472352929384893281939
                    ),
                _s21(0.02083874751496854459875767454686),
                _s111(
                    0.12912437207260851035666365172043,
                    0.060510770643312189063795998961649
                    ),
                _s111(
                    0.087817902777524478493507945986305,
                    0.03591989136357774435488394712972
                    ),
                _s21(0.075762174186790875640220107594336),
                _s21(0.37660198051754571220352610478789),
                _s111(
                    0.27355373802481923697263040817808,
                    0.1640610286702125939455608418683
                    ),
                _s111(
                    0.094477065632374163099096602984735,
                    0.010683263523588402615098865968091
                    ),
                _s21(0.28735583982216495803624567374924),
                _s111(
                    0.41115223330360054298175094895014,
                    0.031132680306373656328092642202672
                    ),
                _s111(
                    0.3594475694982553255274100761827,
                    0.13319262322170576669505468876378
                    ),
                _s111(
                    0.38841008533043709290479191842792,
                    0.074297319237294153700985041421363
                    ),
                _s111(
                    0.32535504334205774221409824249099,
                    0.20516649781476186504903753351473
                    ),
                _s111(
                    0.17478047948106845162209604332071,
                    0.090219606125175470127476810116033
                    ),
                _s111(
                    0.19416461814956333715118429558775,
                    0.04391587290451137373943848868728
                    ),
                _s21(0.47441249890284560954301205127458),
                _s111(
                    0.24832787340457900742014954089071,
                    0.0689638299646016222582127841553
                    ),
                _s111(
                    0.22299204126120093971219669196546,
                    0.12518179814288659476377317186519
                    ),
                _s111(
                    0.14960442354895113916484301004967,
                    0.0046702985612832350654938081970465
                    ),
                _s21(0.11117270128584573407567388067239),
                _s21(0.49186496195768969613724859159279),
                _s111(
                    0.32839141369212664837070754396277,
                    0.048925922503864367545065093845449
                    ),
                _s21(0.15088772276287299186658091037198),
                _s111(
                    0.14238550032454284355271635258644,
                    0.024658249778207267858406454208758
                    ),
                _s111(
                    0.42384029972036885348122877849639,
                    0.0059599432497442501327953694319131
                    ),
                _s21(0.0041579573923190450007642821454266),
                _s21(0.24043491041381679306543528343322),
                _s111(
                    0.26769041934430564520491346765538,
                    0.028518722842285332488257357235445
                    ),
                _s21(0.41551159674483718682799222704665),
                _s21(0.49954218266151180551790778797915),
                _s111(
                    0.34541137895290114981568644555639,
                    0.015562643301015446933777788218238
                    ),
                _s111(
                    0.35344951965766379944866351063409,
                    0.00092571823393377039279135803617303
                    ),
                _s21(0.44856426587870827471387580887119),
                _s111(
                    0.30333563947178920523239966152255,
                    0.099286211618585356118152646378514
                    ),
                _s111(
                    0.021755560320183967673414530001464,
                    0.0039745058725184397403455181669523
                    ),
                _s3(),
                _s111(
                    0.20754686827352048533168144366297,
                    0.013779639405368733490778925623048
                    ),
                _s21(0.19438654357896620992861680102265),
                _s111(
                    0.2135602395179478542516577470156,
                    0.00077213912032697152213689119841515
                    ),
                _s111(
                    0.27938643178222376456868782533148,
                    0.0054617908942933969709742906708799
                    ),
                _s111(
                    0.096866922575547169785110400676201,
                    0.0000047096910185114333517329210283201
                    ),
                ])
            self.weights = numpy.concatenate([
                0.00010168597236713145917918943929861 * numpy.ones(6),
                0.000654047966587507558417734141509 * numpy.ones(6),
                0.0022989665688091104080085169976582 * numpy.ones(3),
                0.0015198116582616399702586923478535 * numpy.ones(6),
                0.0010996805664631989756818136999616 * numpy.ones(3),
                0.0046696376839459277753469832100408 * numpy.ones(6),
                0.002952115630562181900649888768771 * numpy.ones(6),
                0.0039973734536493486862999558412715 * numpy.ones(3),
                0.011740006591257594597396114783444 * numpy.ones(3),
                0.0098609303129047317031462902279874 * numpy.ones(6),
                0.0017262732588417000723275290073269 * numpy.ones(6),
                0.011761207639413641976017588121835 * numpy.ones(3),
                0.0050361825013339413029696500484145 * numpy.ones(6),
                0.0096578802294468123305173772307746 * numpy.ones(6),
                0.0076119150208858810612064138709982 * numpy.ones(6),
                0.011065863887599135196313949543339 * numpy.ones(6),
                0.006503839712649749614965255795077 * numpy.ones(6),
                0.0048326872534850940926992614595542 * numpy.ones(6),
                0.006484225213838527719652299177714 * numpy.ones(3),
                0.0065935144677494465913159465303958 * numpy.ones(6),
                0.0082852665131006666967910687292887 * numpy.ones(6),
                0.0013994293854030885665331842268853 * numpy.ones(6),
                0.0058334336795542317687539090881376 * numpy.ones(3),
                0.0037293709434181781422060646508304 * numpy.ones(3),
                0.0060952640324883482213230881835713 * numpy.ones(6),
                0.0076873106859636912940833092617121 * numpy.ones(3),
                0.003143206148159292408451855513968 * numpy.ones(6),
                0.0022313317933533029898637444936647 * numpy.ones(6),
                0.00022696135677094365139812751927838 * numpy.ones(3),
                0.010813699316269638891106967493837 * numpy.ones(3),
                0.0044357309734419606773044590628338 * numpy.ones(6),
                0.010647444842488946726436119754243 * numpy.ones(3),
                0.00076092616621810198263132006627822 * numpy.ones(3),
                0.0034883273457438320769888937537655 * numpy.ones(6),
                0.00073184587367011224932394773406912 * numpy.ones(6),
                0.0088609393001530497885314118510771 * numpy.ones(3),
                0.0082507974307529598571956607409103 * numpy.ones(6),
                0.00049796757319711846731972276887466 * numpy.ones(6),
                0.012100964930171137684661782846097 * numpy.ones(1),
                0.0028058181373137295814677588248264 * numpy.ones(6),
                0.0094001410371983655921353297437397 * numpy.ones(3),
                0.00057765057423475910923783062591089 * numpy.ones(6),
                0.0019626288054833393547132355726995 * numpy.ones(6),
                0.00028734535460813439989166698086861 * numpy.ones(6),
                ])
        elif index == 35:
            bary = numpy.concatenate([
                _s21(0.3014351686173180200524148644439),
                _s21(0.49164140398289072331460736357694),
                _s21(0.35779313579859713973298277598186),
                _s111(
                    0.055045479691531689579611117351936,
                    0.0020682202154058670016415942896076
                    ),
                _s111(
                    0.44537646491529750324099579373737,
                    0.0089049941861277531821820905320485
                    ),
                _s111(
                    0.15238400252431525193080872153866,
                    0.0021159775630398523960923020410077
                    ),
                _s21(0.26241396010934548382091805229993),
                _s111(
                    0.39751661002823135472949554364254,
                    0.0032612734291311844275561498600889
                    ),
                _s21(0.40281645058578187085235887962273),
                _s111(
                    0.1560787694034761239405739044613,
                    0.045167894813128594397727096006091
                    ),
                _s21(0.47308896550066417355947625525003),
                _s111(
                    0.19729567627655875831792682639624,
                    0.0086060217057567800570629121279222
                    ),
                _s111(
                    0.050093836899926241221776997385247,
                    0.013974353251568308049170014839065
                    ),
                _s111(
                    0.298859063014319120524487525992,
                    0.18747133249989443557156491616345
                    ),
                _s21(0.062868049044065832680233542800642),
                _s21(0.4415974337549303671477268240337),
                _s21(0.019627159983470849299636883494688),
                _s111(
                    0.33476334326533137911809235790699,
                    0.2327223440759484356922673482711
                    ),
                _s111(
                    0.209372824976375347510857623907,
                    0.029001884614684463660842035144076
                    ),
                _s111(
                    0.10656474486138343628019703947397,
                    0.053567482867753488601776312611099
                    ),
                _s111(
                    0.32136262766954331722108603295684,
                    0.0034961943491434885111939940911523
                    ),
                _s111(
                    0.24158702634798757088881075218504,
                    0.0011760427968417102046172834239476
                    ),
                _s111(
                    0.33315309972076488420869617256321,
                    0.11395392874414831587596111082982
                    ),
                _s111(
                    0.15511655481330346273505563403979,
                    0.084938599511554983583797145360595
                    ),
                _s111(
                    0.37101415167746192435395094121426,
                    0.15256518634721628527103649769338
                    ),
                _s21(0.040286868830193527612143047463),
                _s111(
                    0.27088366180777207232496448651032,
                    0.013618678109076597337976131540528
                    ),
                _s111(
                    0.35250674841546647950445568797316,
                    0.056226387443552348068810686446535
                    ),
                _s111(
                    0.40572013277819049302556980270183,
                    0.08175348943968816680635396246678
                    ),
                _s111(
                    0.21894968101499838765381347198724,
                    0.063642430189747335695159450217298
                    ),
                _s111(
                    0.20670924808556452648233984227418,
                    0.11539422083607405126296843810497
                    ),
                _s21(0.098526154668792450682903867781692),
                _s111(
                    0.25660108787909232496822702916498,
                    0.14818711628331062675506165946008
                    ),
                _s111(
                    0.13820629521124998435133828022587,
                    0.018597871061381355114356112917078
                    ),
                _s21(0.13953446298930737044476316588675),
                _s21(0.2215360105096920697277226373179),
                _s111(
                    0.42581020858548103029560200665348,
                    0.032971019926946322659123812719355
                    ),
                _s21(0.18075032213529961668311023452931),
                _s111(
                    0.27972592695224017304684186153074,
                    0.085379551572953806424156331276408
                    ),
                _s111(
                    0.022466683580211183953647165071965,
                    0.003633055842933938782072338151991
                    ),
                _s111(
                    0.083862258524455430386735391721244,
                    0.026020077662363919273085464819533
                    ),
                _s111(
                    0.28480034641446833344435940099593,
                    0.039456926156438076015372834758942
                    ),
                _s111(
                    0.35659410353657383112011774593851,
                    0.019770717704408124003588417784804
                    ),
                _s111(
                    0.097424231115138055861952300854933,
                    0.0050286974156273139122261259952037
                    ),
                _s21(0.0042772072684830632763542286437298),
                _s21(0.49991322248622513537238172173754),
                ])
            self.weights = numpy.concatenate([
                0.0057479340296701791281573065508254 * numpy.ones(3),
                0.0022963680920994285422540890110283 * numpy.ones(3),
                0.0079953508418643076266061200118748 * numpy.ones(3),
                0.00045419938237539701931339815246452 * numpy.ones(6),
                0.0019717493309841691952286686134243 * numpy.ones(6),
                0.00078266630543935302281339612849867 * numpy.ones(6),
                0.0087901858933277923698002331728089 * numpy.ones(3),
                0.0012646931023971697095664658608585 * numpy.ones(6),
                0.0083203666123256195924429803908734 * numpy.ones(3),
                0.003557256260909509655614058283622 * numpy.ones(6),
                0.005275188285161489914231671364736 * numpy.ones(3),
                0.0018428326631755046694244397742062 * numpy.ones(6),
                0.0013289799356541578355636587940591 * numpy.ones(6),
                0.0090234128551062183850134826236014 * numpy.ones(6),
                0.0029682066390516002096851310871601 * numpy.ones(3),
                0.0076051825305363510536683057928173 * numpy.ones(3),
                0.0010131210327380781996863296984869 * numpy.ones(3),
                0.0086055044940838447993912569060791 * numpy.ones(6),
                0.0036104679969938313898900699494671 * numpy.ones(6),
                0.0035235607672107945858109341824623 * numpy.ones(6),
                0.0013780765220063653047091090999738 * numpy.ones(6),
                0.00069502958595185644151735265368314 * numpy.ones(6),
                0.0081297205980779889536535339696769 * numpy.ones(6),
                0.0054896563598667055299310669387661 * numpy.ones(6),
                0.0082359873228582471797789893959744 * numpy.ones(6),
                0.0020648356265983247943124962730898 * numpy.ones(3),
                0.0027797735193625662231629443049078 * numpy.ones(6),
                0.0063317818814833934035575008701378 * numpy.ones(6),
                0.0069457107199415825055017156563808 * numpy.ones(6),
                0.005565941115686088797319006942461 * numpy.ones(6),
                0.0073376101567626745817678730576262 * numpy.ones(6),
                0.0050036396155573284191284593680021 * numpy.ones(3),
                0.0087114266235608622611548220359297 * numpy.ones(6),
                0.0027367494327220168972285176491668 * numpy.ones(6),
                0.0067190798124173888459733255991639 * numpy.ones(3),
                0.0091619089409372974037658067350699 * numpy.ones(3),
                0.00470165758093218034173765977786 * numpy.ones(6),
                0.00837190663701295860185128433997 * numpy.ones(3),
                0.0071750089599327115901370521392923 * numpy.ones(6),
                0.00047324697380048877766017221194308 * numpy.ones(6),
                0.0025346621113677792987428456705843 * numpy.ones(6),
                0.0050701573000113947427067796555698 * numpy.ones(6),
                0.0040540720868840226397804475193782 * numpy.ones(6),
                0.0012544783872410158196130090166232 * numpy.ones(6),
                0.00024033373740484100604262740998811 * numpy.ones(3),
                0.00062758434107056451114671085808747 * numpy.ones(3),
                ])
        elif index == 36:
            bary = numpy.concatenate([
                _s111(
                    0.25500233607541949512398192492591,
                    0.03229226397376931487287858188622
                    ),
                _s111(
                    0.27112151360334001811107274181685,
                    0.094140209804408110823992306066949
                    ),
                _s111(
                    0.09268860448487556523855437307129,
                    0.06499367210424015969256280843294
                    ),
                _s111(
                    0.32945437658455049271590293810876,
                    0.073713664453418720763804061931842
                    ),
                _s111(
                    0.28837884113044209474991357720998,
                    0.050064988146891895002181963064969
                    ),
                _s111(
                    0.39217241640517385644856715365417,
                    0.065896338474894637343624276669001
                    ),
                _s111(
                    0.0076222500155876565014384583058651,
                    0.0012341185277631950702921614604898
                    ),
                _s21(0.41235734031363032320408416512899),
                _s111(
                    0.1782079366279756831465047335389,
                    0.10569582732760981950379029364796
                    ),
                _s111(
                    0.31024893145991363321152659319331,
                    0.12387167039555084420968853882356
                    ),
                _s21(0.4390483832458435331647405626611),
                _s111(
                    0.082915702133125166901044285685691,
                    0.0028073976185661223652120566753582
                    ),
                _s111(
                    0.37963867934694315360652864417291,
                    0.1106897783260211745887032399063
                    ),
                _s21(0.23505232545821376610592419608828),
                _s111(
                    0.19806829211454337021860101491218,
                    0.002756422613705629574247549278643
                    ),
                _s111(
                    0.14794794532207304233683730518179,
                    0.072051699096015864153101440488593
                    ),
                _s111(
                    0.22920276655346734799882429630095,
                    0.014428611570684238868189058010906
                    ),
                _s111(
                    0.35747385034438452143846881569617,
                    0.0027425666974118116972640999308381
                    ),
                _s111(
                    0.12524967222943776240902776924218,
                    0.039024329842215862904065202155452
                    ),
                _s111(
                    0.020980564600406259183373944408834,
                    0.0074746634659332896688891807535419
                    ),
                _s111(
                    0.26913161174400975632479692262961,
                    0.1755717786015191928604308593436
                    ),
                _s21(0.15082173477989161121181791008303),
                _s111(
                    0.34621987242189356873119205039997,
                    0.17122675885075284887776010462274
                    ),
                _s111(
                    0.44238774265806850271276008451039,
                    0.034306143454093658790713219254332
                    ),
                _s111(
                    0.27463388372265885201103376219503,
                    0.0030296355118470474300145636130801
                    ),
                _s111(
                    0.22370128465514409439692358922009,
                    0.072486958455644065671930296487127
                    ),
                _s111(
                    0.2264729475116298008371642825192,
                    0.13381702334066388620664456221865
                    ),
                _s111(
                    0.16336196682937046306520101881712,
                    0.017674220238250244727740776877251
                    ),
                _s111(
                    0.19079127285437891237885745352034,
                    0.042208217039582026008283813808654
                    ),
                _s111(
                    0.35816587063143068182878537693962,
                    0.034416320396989541191627655463323
                    ),
                _s111(
                    0.099384989205052216371792358239105,
                    0.01745299026102043350792153181998
                    ),
                _s111(
                    0.45023289789165244461487407689061,
                    0.0012278348447230405457289117058109
                    ),
                _s21(0.46427377763467002371082454568976),
                _s111(
                    0.30486029298667931229814434581011,
                    0.23626319006630857722958637938601
                    ),
                _s111(
                    0.31358360627069516476830742039205,
                    0.015522928528731724982861457640047
                    ),
                _s111(
                    0.13480628826668292454020641349379,
                    0.0038480457708450127921411046128769
                    ),
                _s21(0.19116136187730902246629035539923),
                _s111(
                    0.065828192432819855015812422071014,
                    0.038072469323013046746649004011892
                    ),
                _s21(0.0293500490934270145562944875861),
                _s21(0.30814521899109328546075933536758),
                _s111(
                    0.05319284391945779635039490706026,
                    0.012170254153940130281843007793491
                    ),
                _s21(0.11117826478838099324275910603754),
                _s111(
                    0.040647962635958589709322113351126,
                    0.00085855904468391444106475254919908
                    ),
                _s111(
                    0.40644756721228486943689468907724,
                    0.012335837810327154256797019919775
                    ),
                _s21(0.38100234960344190051056077420231),
                _s21(0.49448099674770043490118158909823),
                ])
            self.weights = numpy.concatenate([
                0.0030616831670681794051376732899791 * numpy.ones(6),
                0.0049232719752900455260606407527321 * numpy.ones(6),
                0.0033237208433952315770482284839712 * numpy.ones(6),
                0.0047427364930337531703124254921969 * numpy.ones(6),
                0.0041138012926497477751177015397185 * numpy.ones(6),
                0.0052353425228937827953391587851771 * numpy.ones(6),
                0.00013033992845265019696432858063866 * numpy.ones(6),
                0.0076704945179086273098565856825717 * numpy.ones(3),
                0.0052277042839383569906168381079187 * numpy.ones(6),
                0.0066311375415485319144702122281638 * numpy.ones(6),
                0.0068762957093961478172151852823822 * numpy.ones(3),
                0.00069523463971707374574510331189471 * numpy.ones(6),
                0.0069183724873925471126402252308174 * numpy.ones(6),
                0.008067997716638237796745733816295 * numpy.ones(3),
                0.0010177068493549157286667694008545 * numpy.ones(6),
                0.0045309912970645343018541026675868 * numpy.ones(6),
                0.0023792804956385682868685389992833 * numpy.ones(6),
                0.0012168791097054403626811650636444 * numpy.ones(6),
                0.0031130272403783083819790673047488 * numpy.ones(6),
                0.00064818533281312525582010282492515 * numpy.ones(6),
                0.0083451647460545622752384456750151 * numpy.ones(6),
                0.0061710606814515137457975776991519 * numpy.ones(3),
                0.0089038028492518081625957780578287 * numpy.ones(6),
                0.0047638139693882780872377741060619 * numpy.ones(6),
                0.0012331434914374479983841638391279 * numpy.ones(6),
                0.0052810362353660750810300860244754 * numpy.ones(6),
                0.0072039581630509769674214733375883 * numpy.ones(6),
                0.0025500464663513038798351129230314 * numpy.ones(6),
                0.0040931211574026578149172156263739 * numpy.ones(6),
                0.0044693917436126364645808129074091 * numpy.ones(6),
                0.0021012224230352472056395187787744 * numpy.ones(6),
                0.00082208910731422006807610275045946 * numpy.ones(6),
                0.0066612336212475384540749356307419 * numpy.ones(3),
                0.010213345615529134321334461838031 * numpy.ones(6),
                0.0030402213548300950278881436299552 * numpy.ones(6),
                0.0011096101536170602471796121750917 * numpy.ones(6),
                0.0080932250756890799815324013100259 * numpy.ones(3),
                0.0026669011660034517595990037022139 * numpy.ones(6),
                0.0015674237111508730546185300476561 * numpy.ones(3),
                0.011008276161713617255886772347248 * numpy.ones(3),
                0.0013356087640650201317228672428898 * numpy.ones(6),
                0.0054235489241325678702158688457614 * numpy.ones(3),
                0.00026617782389344208391368574367552 * numpy.ones(6),
                0.0028876869913215207259343347916528 * numpy.ones(6),
                0.010506346006817103243921225387974 * numpy.ones(3),
                0.0028959157634685651437667668557156 * numpy.ones(3),
                ])
        elif index == 37:
            bary = numpy.concatenate([
                _s21(0.10230393215776038942336921827105),
                _s111(
                    0.31602409649663992132444343586255,
                    0.033829090682721268449822933145531
                    ),
                _s21(0.035020235875468501947471934532476),
                _s111(
                    0.10535941196028194987542464069424,
                    0.033145989556188382211881896964902
                    ),
                _s111(
                    0.074981021784113832257157119648981,
                    0.0020334760761173136690507308887381
                    ),
                _s111(
                    0.1194865933803552373373337100103,
                    0.055510403162458813953125341777516
                    ),
                _s111(
                    0.31923516504598916029116571847295,
                    0.23806084711702476533693995368031
                    ),
                _s111(
                    0.39178825772367013414193052935199,
                    0.075270548937509795797800404121681
                    ),
                _s21(0.25924349045654960210178277595049),
                _s111(
                    0.31754340721419831129637249452325,
                    0.053620831629113071062112349858354
                    ),
                _s111(
                    0.3267645245927317678154507714348,
                    0.091294767095042794653415499187448
                    ),
                _s111(
                    0.3187496668271801607966001687247,
                    0.13511574907028460762939269907494
                    ),
                _s111(
                    0.067457879081386864636374821246051,
                    0.040770123756678612942619718537768
                    ),
                _s111(
                    0.018379440921081184627033715785289,
                    0.002990277876563861959906641513652
                    ),
                _s111(
                    0.047656935862215985420451419435748,
                    0.017250960453807498868227052570945
                    ),
                _s111(
                    0.33572278703801007598718593613099,
                    0.18111640919553191219909471280818
                    ),
                _s21(0.13944123587229279152808070265133),
                _s21(0.4611684588908100347577237831536),
                _s111(
                    0.13567160791611764067782360237152,
                    0.09013027184829002960809478165719
                    ),
                _s111(
                    0.041769361435969812732607420826123,
                    0.004017097087231190245334506243355
                    ),
                _s21(0.37915371795503783203113389284347),
                _s111(
                    0.19126373743807623012081472000057,
                    0.080965320229254034403700906330752
                    ),
                _s111(
                    0.44892032995664287663429877444761,
                    0.016881458510130714006881120716625
                    ),
                _s111(
                    0.39708337913498853298483985216734,
                    0.041476678533636860634581046259786
                    ),
                _s111(
                    0.086226344062887231612686870821755,
                    0.013734586162212104739448389836039
                    ),
                _s21(0.017341175363902502822903580975),
                _s111(
                    0.28148795994228197454801015115455,
                    0.017918604330291788230907983487799
                    ),
                _s111(
                    0.19237683474425930131278756272775,
                    0.12847139005679287725204058119443
                    ),
                _s111(
                    0.25179424751743315932778737739543,
                    0.13732115973159557598642590154573
                    ),
                _s21(0.077988930180642711450465488650418),
                _s111(
                    0.36503479493632803478121195316509,
                    0.017441826494383764697977357831306
                    ),
                _s21(0.47924493555692034568432501757283),
                _s21(0.41206140547968950630502293859776),
                _s111(
                    0.25715663555025350877736770645366,
                    0.086198170278100264486990048719183
                    ),
                _s21(0.49841057441395697045007382252957),
                _s21(0.30912594868535053289933569022512),
                _s111(
                    0.1241816119341110627857762107449,
                    0.003411267800834652119298197794776
                    ),
                _s111(
                    0.17165072780936918349351448130016,
                    0.044555516711789470972751711322934
                    ),
                _s111(
                    0.39772310919828560748010180806439,
                    0.12209754888456262625301182744906
                    ),
                _s21(0.0036715668219832217918029162044199),
                _s111(
                    0.25789664566985651780429310061083,
                    0.19571037123642808600833472278447
                    ),
                _s111(
                    0.1441500744811936133110523090234,
                    0.018311391979833948640216477205167
                    ),
                _s111(
                    0.41225956100179108325089036713772,
                    0.0032622580352037452876055403983807
                    ),
                _s111(
                    0.25242846479143637152849272674502,
                    0.0034754415073156802927648280020431
                    ),
                _s111(
                    0.32936692748577454325973889137098,
                    0.0035387900825817599886210801720439
                    ),
                _s111(
                    0.20958148170051542247166905901425,
                    0.018844939055752260994674901623483
                    ),
                _s111(
                    0.18364016944623227334099095411648,
                    0.0036381021086962816022672781382663
                    ),
                _s111(
                    0.24169623246933752268131937040061,
                    0.04575328696131737766480630282354
                    ),
                _s21(0.18936426120019623397689321027703),
                ])
            self.weights = numpy.concatenate([
                0.0015139897718580163925513410451786 * numpy.ones(3),
                0.0020991920973494748879219895909269 * numpy.ones(6),
                0.0011046836690446153992533648924592 * numpy.ones(3),
                0.0019226251721820787724784008633998 * numpy.ones(6),
                0.00048044301000187699568215645014258 * numpy.ones(6),
                0.0030325474600894931053549135816003 * numpy.ones(6),
                0.0073189063292417020513202317221093 * numpy.ones(6),
                0.005388318842980985416566638322501 * numpy.ones(6),
                0.0075850923563705291487037660291836 * numpy.ones(3),
                0.0045692159189814848278286233510779 * numpy.ones(6),
                0.0057203599064584805098049290958191 * numpy.ones(6),
                0.0067120048490856740426997104755126 * numpy.ones(6),
                0.0023170775093131488999796971475839 * numpy.ones(6),
                0.00030866926197053260282205576045454 * numpy.ones(6),
                0.0012047027694410606635695822973198 * numpy.ones(6),
                0.0081637462521435939175320027515577 * numpy.ones(6),
                0.00525006015432878039824699342056 * numpy.ones(3),
                0.0058808346792778311784556126424908 * numpy.ones(3),
                0.0043645240281899114999569259313846 * numpy.ones(6),
                0.00054979655064731875557506026415517 * numpy.ones(6),
                0.0086445565740853366755609257457868 * numpy.ones(3),
                0.0052464381586921763475735533169637 * numpy.ones(6),
                0.0032584894563912941599868111744776 * numpy.ones(6),
                0.0047812349883803963730942991713056 * numpy.ones(6),
                0.0016087719028868608343062715856129 * numpy.ones(6),
                0.00087758725215081531225910908594101 * numpy.ones(3),
                0.0029825665052993322539302365125833 * numpy.ones(6),
                0.0060939488435528657106462390220189 * numpy.ones(6),
                0.0070873407460566485361930884356582 * numpy.ones(6),
                0.0035862336027496868497337254070107 * numpy.ones(3),
                0.0031497526280811912079238692407605 * numpy.ones(6),
                0.0049992428284758221422861737996442 * numpy.ones(3),
                0.009078842505539653648648609946981 * numpy.ones(3),
                0.0062518699531209620440824547446847 * numpy.ones(6),
                0.0014156781091945123893342702595185 * numpy.ones(3),
                0.010206636649173599581030321807386 * numpy.ones(3),
                0.00094930532254722946427351138770794 * numpy.ones(6),
                0.0040312520842928054905749973079593 * numpy.ones(6),
                0.0081722175482304027047711293578332 * numpy.ones(6),
                0.00017329533578085930875551508240657 * numpy.ones(3),
                0.0089586590287218429099732883666605 * numpy.ones(6),
                0.0024620987414923229849447070892849 * numpy.ones(6),
                0.0014212929213219290304715363224864 * numpy.ones(6),
                0.0013051564959186817691005804547067 * numpy.ones(6),
                0.001445909527788727993878020718992 * numpy.ones(6),
                0.0029345950580795215034177173929056 * numpy.ones(6),
                0.0011964590052698317399322571238026 * numpy.ones(6),
                0.004873504419922933871006644364971 * numpy.ones(6),
                0.0082906132570537271501653427749507 * numpy.ones(3),
                ])
        elif index == 38:
            bary = numpy.concatenate([
                _s21(0.45907213603176458246435226305336),
                _s111(
                    0.28376551470174071081775285162444,
                    0.052297113902011395792653636262807
                    ),
                _s21(0.39326800625536008887161685490857),
                _s111(
                    0.12851656496969396603708894900162,
                    0.099079762701841230091659280756317
                    ),
                _s21(0.46818258920693570436742521357158),
                _s111(
                    0.041307504142861442947527895044092,
                    0.002064405149305127390232704359356
                    ),
                _s111(
                    0.089271980053414487687220423695186,
                    0.036730175649348289528728823611487
                    ),
                _s111(
                    0.33720529231798726867880776700086,
                    0.22594582448778563812468932838853
                    ),
                _s111(
                    0.27385137412417876884163127435613,
                    0.032701181990053704912142566151332
                    ),
                _s111(
                    0.039916262629419480236709324769684,
                    0.013069509044329950644226919633182
                    ),
                _s21(0.017504504046267431378756512961672),
                _s111(
                    0.30810683275938431868662866211106,
                    0.07371761457101587514331943943538
                    ),
                _s21(0.14548011465286077639870973551104),
                _s111(
                    0.21962471313843564785243834640827,
                    0.002703882214543481675509477044714
                    ),
                _s111(
                    0.27346460711428144393647832660225,
                    0.219158263446028715003193404479
                    ),
                _s111(
                    0.06868300285555269396743358712642,
                    0.02072281188151401137980906844503
                    ),
                _s111(
                    0.18301359984017947616192388543367,
                    0.10062172089207602255388653984331
                    ),
                _s21(0.0030471742663605824001934196515523),
                _s111(
                    0.22879948322992527113739916388472,
                    0.069489817833509195680446064902791
                    ),
                _s21(0.058736160112022102583728238121706),
                _s111(
                    0.10188271148693912620437139617314,
                    0.068377355422513814088301691932677
                    ),
                _s111(
                    0.24356367339581708170719827968794,
                    0.014320441363316447634213302215118
                    ),
                _s111(
                    0.20588771651790122386589144050365,
                    0.038967856742294220042596183968835
                    ),
                _s111(
                    0.20186950888760312532158978505371,
                    0.14939646191794411417287377701137
                    ),
                _s111(
                    0.15911602101054522284614758654491,
                    0.062228211818953253430972617156012
                    ),
                _s111(
                    0.17794763585559181479553300407781,
                    0.017619241470485949705901195353563
                    ),
                _s111(
                    0.16203674794969208525148697426643,
                    0.0039477276652124153357567547372951
                    ),
                _s111(
                    0.29185498385852457923248010865898,
                    0.0029311798861728812161383527799786
                    ),
                _s111(
                    0.37187869174958914418080707858889,
                    0.0027779809824385259329705893488835
                    ),
                _s111(
                    0.13496937265246532531371790350104,
                    0.034322185835850192084000712112331
                    ),
                _s21(0.037383557501746567662958732953433),
                _s21(0.28550477829126408600697436719347),
                _s111(
                    0.25245342990642194472137022807325,
                    0.10872403934192483832716404154053
                    ),
                _s21(0.42033350177001423005917740430889),
                _s111(
                    0.32423516410830452563673290947146,
                    0.014949107908045373682707924156078
                    ),
                _s21(0.20768942165386299240854257245014),
                _s21(0.49219627295531079682666670364117),
                _s111(
                    0.43741619064662571006510949966458,
                    0.037324077159263422625562800664488
                    ),
                _s111(
                    0.45599880366120932438829635013292,
                    0.0029632719316369204765033320555992
                    ),
                _s111(
                    0.075487285700568336082153242582021,
                    0.0042029781607734603186810846424577
                    ),
                _s111(
                    0.32605412665294664076256852932457,
                    0.11427418738045437433721280358451
                    ),
                _s111(
                    0.4011995885054882178111566528747,
                    0.11288808531292145207613016639107
                    ),
                _s111(
                    0.26876066034126400519596994400147,
                    0.16019667789975580569238278831731
                    ),
                _s111(
                    0.40601470257387719695930544026618,
                    0.015155128640289381967377890771781
                    ),
                _s21(0.35656559398302617582716178586648),
                _s111(
                    0.3835103468343401980349466858887,
                    0.07004218150043782955028548618264
                    ),
                _s111(
                    0.016447064636943181064621482020815,
                    0.0034429885015807508054199107497701
                    ),
                _s111(
                    0.35429778616269475007748091821817,
                    0.037667614478677257016924184444153
                    ),
                _s111(
                    0.3419951757983160694881959886541,
                    0.16599842667671561153765721579527
                    ),
                _s111(
                    0.11528889107927841894306450821189,
                    0.012032531111165113037011713487608
                    ),
                _s111(
                    0.11652810181240256723345702106015,
                    0.00029822168152010235911031555613536
                    ),
                ])
            self.weights = numpy.concatenate([
                0.0031639644402839227298908428609207 * numpy.ones(3),
                0.0024850052311426517321248205426254 * numpy.ones(6),
                0.0058616582191086431077964854381819 * numpy.ones(3),
                0.0036964891607068277574099060788353 * numpy.ones(6),
                0.0039302169037011433084560806865839 * numpy.ones(3),
                0.00034770439899369939368174188412526 * numpy.ones(6),
                0.002073826873940778004686440445841 * numpy.ones(6),
                0.0073710410900112498522942627078269 * numpy.ones(6),
                0.00312539811443985276973423217093 * numpy.ones(6),
                0.00090617621445205684343413889600718 * numpy.ones(6),
                0.00069946976732960115417369895935593 * numpy.ones(3),
                0.0049956277697466920536097219938639 * numpy.ones(6),
                0.0052024225110438605290733322349531 * numpy.ones(3),
                0.00092566819884991224618131527941813 * numpy.ones(6),
                0.0082821664626929521715105354314857 * numpy.ones(6),
                0.0015368424658073105368511965963428 * numpy.ones(6),
                0.0052227407388143726957324024062719 * numpy.ones(6),
                0.00012404666265310360371433615627593 * numpy.ones(3),
                0.0047183927884199089052249645949856 * numpy.ones(6),
                0.0025063822165144870368436682333483 * numpy.ones(3),
                0.0036021170049528564213546920370779 * numpy.ones(6),
                0.0023702904673933514753009525302747 * numpy.ones(6),
                0.0036104898333757971884137225935053 * numpy.ones(6),
                0.006540454530745502548669401340144 * numpy.ones(6),
                0.0040792214803126023451285549012379 * numpy.ones(6),
                0.0023256389181058352827389755576104 * numpy.ones(6),
                0.0010057440220832629508843072070472 * numpy.ones(6),
                0.0011431617714117464945911561206274 * numpy.ones(6),
                0.0011884613225082550193495069454959 * numpy.ones(6),
                0.0029052189452109211836474167214887 * numpy.ones(6),
                0.0017075502039157281182790714243923 * numpy.ones(3),
                0.0089856095188049886748139170378525 * numpy.ones(3),
                0.0063267151975865725200942344701603 * numpy.ones(6),
                0.0079377042751204553085299904906044 * numpy.ones(3),
                0.0027256862207201200479880562173351 * numpy.ones(6),
                0.0077485531926674767056114973450406 * numpy.ones(3),
                0.0030358975545344329827397684291569 * numpy.ones(3),
                0.0045215915998051626774377449433979 * numpy.ones(6),
                0.0012891742774294219595797580242099 * numpy.ones(6),
                0.00081668391817278049605951896153405 * numpy.ones(6),
                0.006947779323552220307689844783423 * numpy.ones(6),
                0.0070986179394903163815922515555042 * numpy.ones(6),
                0.0076919203485779147605404796056596 * numpy.ones(6),
                0.0029252190187241046260390813329804 * numpy.ones(6),
                0.0099759936103825657999716886649749 * numpy.ones(3),
                0.0060566969964857717345687776300953 * numpy.ones(6),
                0.0003337056331188373011195781387518 * numpy.ones(6),
                0.0045294689209488320941813840929225 * numpy.ones(6),
                0.0084025374143197728849146012735994 * numpy.ones(6),
                0.0018109718524034119062834535719287 * numpy.ones(6),
                0.00029228566318282656607634810127754 * numpy.ones(6),
                ])
        elif index == 39:
            bary = numpy.concatenate([
                _s111(
                    0.23827520858455705038397719741145,
                    0.0043744824382612274196383751241506
                    ),
                _s111(
                    0.23143707588900051816654194469429,
                    0.26627272279608107339709739503404
                    ),
                _s111(
                    0.35068273415715828739582164961232,
                    0.23609944743651987040709561258668
                    ),
                _s111(
                    0.27780891771928405414889640035633,
                    0.0026623397478916491254344114112373
                    ),
                _s21(0.47240948215055669425312313089635),
                _s111(
                    0.32120348268261683383358629685352,
                    0.20743961193494545733536192077095
                    ),
                _s21(0.45750796867382429011401231492906),
                _s111(
                    0.45742421685298210211387501105706,
                    0.011909152476788598355337157938001
                    ),
                _s111(
                    0.19103834285218465520903205664018,
                    0.0021569623388998319001402730817768
                    ),
                _s21(0.3538472919168833466734091650188),
                _s21(0.49910614986029900954650242398451),
                _s111(
                    0.20635194073283629173806489407905,
                    0.10725749375164579347421286341623
                    ),
                _s111(
                    0.26123838278747194444895918719461,
                    0.11243839329148250035410179351293
                    ),
                _s21(0.11333165819636543872795693856663),
                _s111(
                    0.15938889962312034615969296443533,
                    0.11409032499671963997088040254304
                    ),
                _s111(
                    0.39886951322406686683127958458232,
                    0.061944841916062421762753059306575
                    ),
                _s111(
                    0.31227257927336089877941097982959,
                    0.12196288322528455881431048002157
                    ),
                _s111(
                    0.36443626235847638509268293460167,
                    0.038644796347540205306496418869573
                    ),
                _s111(
                    0.23056483424044315626815217349059,
                    0.017233971768566536009215235128396
                    ),
                _s21(0.015778704722834468442350599537996),
                _s111(
                    0.039268283381893976356782728481549,
                    0.014562370303059139096219186396631
                    ),
                _s111(
                    0.35825773354772274635994754803344,
                    0.15416110999734943716817815846501
                    ),
                _s111(
                    0.38050788254831254077247283949223,
                    0.10114013067894874119347423900886
                    ),
                _s111(
                    0.13188657558660895023702379087444,
                    0.073285816648122131939096838926252
                    ),
                _s111(
                    0.17263025506888606496244805085098,
                    0.01351746339721455784076503927526
                    ),
                _s111(
                    0.070775635890008583472203794998829,
                    0.043148646381222143866180680238771
                    ),
                _s111(
                    0.047382473120469954625202431514324,
                    0.0028483388383531039287482248666423
                    ),
                _s21(0.16400918476660571619343466674949),
                _s21(0.40847156072913401804565458442948),
                _s111(
                    0.29719324829574491099550682901776,
                    0.04039745419698131941607514780668
                    ),
                _s111(
                    0.44376523995359830397523719988839,
                    0.031002548343527617449598378264777
                    ),
                _s111(
                    0.22985489024982292933439723666506,
                    0.039285463900080663298515858079419
                    ),
                _s21(0.43654165510795662722306043069148),
                _s111(
                    0.25152542713890719223532378028988,
                    0.071564132049667623029381978385024
                    ),
                _s111(
                    0.2206287919327525128617009688086,
                    0.15439459178077457784060841823078
                    ),
                _s21(0.037032651339393698546820159045172),
                _s111(
                    0.41640318204681779345135626873083,
                    0.0027734733727749882802797390149083
                    ),
                _s111(
                    0.27803724358402760076882184328971,
                    0.1718458690352429704841103191598
                    ),
                _s111(
                    0.34207583221450852168147569345717,
                    0.0032524102580492533777234823262568
                    ),
                _s21(0.2117317495651855983232974105996),
                _s111(
                    0.37665152274202426758614792438954,
                    0.016081199477227620865872056790601
                    ),
                _s111(
                    0.11229791674123558929202317843003,
                    0.043659198704887485527565667310017
                    ),
                _s111(
                    0.18614115584177375069434082600742,
                    0.068295692361052732262259100393885
                    ),
                _s111(
                    0.32170339827518716900807278953162,
                    0.074722503110879761542174481920342
                    ),
                _s111(
                    0.16601268517886673259977226571593,
                    0.036155374364333908870856179401066
                    ),
                _s111(
                    0.019357695235428503585920745517842,
                    0.002833090699241162738003860396787
                    ),
                _s111(
                    0.071990409225125435881428998295629,
                    0.017949883530718709680914785753922
                    ),
                _s21(0.28695728845906866109550663173504),
                _s111(
                    0.13333271687989791500966184544381,
                    0.003192539870571775120092419097311
                    ),
                _s111(
                    0.29944552917299389954078191907998,
                    0.016292457419014677193560665540302
                    ),
                _s111(
                    0.11710788518915331947538029318979,
                    0.017915399181348695781387412134915
                    ),
                _s21(0.080676582882942615408231574639043),
                _s111(
                    0.085367601656265483759305391394938,
                    0.0035672493041341532013310710538541
                    ),
                _s21(0.0036432199207300034628450006252001),
                ])
            self.weights = numpy.concatenate([
                0.00078240303283956423874118232700277 * numpy.ones(6),
                0.005425925878952458884376059501337 * numpy.ones(6),
                0.0064335649832209141699779298417713 * numpy.ones(6),
                0.00076451063763166986110760180480544 * numpy.ones(6),
                0.0038004729776006585162670326660338 * numpy.ones(3),
                0.0062594886810263873207283050180193 * numpy.ones(6),
                0.0048090357401686065334371416842036 * numpy.ones(3),
                0.0021620874068945465643377797414776 * numpy.ones(6),
                0.00069085397832228005825083108176762 * numpy.ones(6),
                0.0080789945453473424513579123383499 * numpy.ones(3),
                0.00083847668182917052743023153085638 * numpy.ones(3),
                0.0043810960857082160611010917089701 * numpy.ones(6),
                0.0048977390707915226385893475660378 * numpy.ones(6),
                0.0037089356213761816567919420002911 * numpy.ones(3),
                0.0044309896508940548179863560755177 * numpy.ones(6),
                0.0046998784286279422357104212476523 * numpy.ones(6),
                0.0056949471589668717571259742494586 * numpy.ones(6),
                0.0037232299887488257850470124114754 * numpy.ones(6),
                0.0022711066220262451013372422819537 * numpy.ones(6),
                0.0006686069116893476495819287410296 * numpy.ones(3),
                0.0009501856551478382635461904407108 * numpy.ones(6),
                0.0074331006305705508167397120552734 * numpy.ones(6),
                0.006186430533741882270955505764133 * numpy.ones(6),
                0.0038374222817438528830850729477778 * numpy.ones(6),
                0.0018958729657530738865265265972127 * numpy.ones(6),
                0.0023158771686323537354188658780581 * numpy.ones(6),
                0.00048827076220379704873481858236667 * numpy.ones(6),
                0.0058035379196035618906751309999051 * numpy.ones(3),
                0.0076992066322762727193138977691469 * numpy.ones(3),
                0.0039332470075348944242742963039722 * numpy.ones(6),
                0.0037369831874342804157935161465231 * numpy.ones(6),
                0.0035306268703355609509228585229893 * numpy.ones(6),
                0.0067765032180855106279358142166106 * numpy.ones(3),
                0.004873301490613901092950841591629 * numpy.ones(6),
                0.006069692421926069587535345061572 * numpy.ones(6),
                0.0015920496320692565395466825613666 * numpy.ones(3),
                0.001101763580750840065630979151604 * numpy.ones(6),
                0.007309182815522299783383894356363 * numpy.ones(6),
                0.0011736195277390052042204533177438 * numpy.ones(6),
                0.0071255638341492700135129158662022 * numpy.ones(3),
                0.0028277464204102066432515679806619 * numpy.ones(6),
                0.0028118195580629150857478894253126 * numpy.ones(6),
                0.0043776217728844006198329784236964 * numpy.ones(6),
                0.0057797851749689305413366546442134 * numpy.ones(6),
                0.0033216883150330844141927693187156 * numpy.ones(6),
                0.00032342143857673339243114765827105 * numpy.ones(6),
                0.0015599585117589139172884565863101 * numpy.ones(6),
                0.0087287285098107616248334630217017 * numpy.ones(3),
                0.00086845965730168345490639181779516 * numpy.ones(6),
                0.0027516782565965268599502366061166 * numpy.ones(6),
                0.0020505862443132491117944585982003 * numpy.ones(6),
                0.0037076853680653037593553422033865 * numpy.ones(3),
                0.00078415012564603347596723080385263 * numpy.ones(6),
                0.0001749077815533339416223108576126 * numpy.ones(3),
                ])
        elif index == 40:
            bary = numpy.concatenate([
                _s111(
                    0.16717783835571051963376246575724,
                    0.0015279490490131733337319693917183
                    ),
                _s111(
                    0.37147809944523749626074442017884,
                    0.014971822472141024409886684219823
                    ),
                _s111(
                    0.15902388471237365821685629667575,
                    0.0068042587991277116749434101881472
                    ),
                _s111(
                    0.046217932928266555977785132837634,
                    0.035313974701554554297259679633073
                    ),
                _s111(
                    0.11486516738008069319561733674697,
                    0.0015839864476045920468400260750598
                    ),
                _s111(
                    0.33424600818471447530804955851607,
                    0.021516672084481753144593628514516
                    ),
                _s111(
                    0.072371164948390053040093274814338,
                    0.0020252616925630734951931132829451
                    ),
                _s111(
                    0.22368424666478873995749776522344,
                    0.015361203749086823994698514037345
                    ),
                _s111(
                    0.11112901699058496684636908556723,
                    0.010165943947914912867226370702714
                    ),
                _s111(
                    0.27644454180779090022633980742958,
                    0.030013853783829215137711342322996
                    ),
                _s111(
                    0.1692633554656757998035573552258,
                    0.018970042640021485971861750534705
                    ),
                _s111(
                    0.11965095524875834453526602185146,
                    0.026740784222126943483817175344592
                    ),
                _s111(
                    0.28301909712926866877146079203038,
                    0.00085381869353749576421839373610502
                    ),
                _s111(
                    0.070404550151146779275221119071937,
                    0.012258470043981681509584507896465
                    ),
                _s111(
                    0.21647209875141211747289770860616,
                    0.037170047113496947680197246145002
                    ),
                _s111(
                    0.0769268815632590626288472079731,
                    0.032194200663470392820086017038539
                    ),
                _s111(
                    0.42922273203855127214686001792993,
                    0.0089321527113153766254269780599656
                    ),
                _s21(0.015703745363793126601786941725552),
                _s21(0.29388996289744307070657793284023),
                _s21(0.37010286842374164975184695857811),
                _s111(
                    0.32514792488520486809325579835286,
                    0.22404521522749966650967234326754
                    ),
                _s111(
                    0.27827271028693636604218152826737,
                    0.18762621699067558324948572151168
                    ),
                _s111(
                    0.2852102535187100370863227626014,
                    0.009758521323868089630552017514736
                    ),
                _s111(
                    0.27272563182983105891574364428152,
                    0.058706290621673538518824589174707
                    ),
                _s111(
                    0.41546959477123946491752008305495,
                    0.030741999686160222974317286728017
                    ),
                _s111(
                    0.22057744653370911803847154732877,
                    0.0032567536371101491689822052053182
                    ),
                _s111(
                    0.34288562039435878660503784233979,
                    0.045187622960489076346909268802616
                    ),
                _s111(
                    0.23150386061168030327626800528008,
                    0.15081851022405672349522263356782
                    ),
                _s111(
                    0.038101712964918257785908248570516,
                    0.016673674037479828961790322283399
                    ),
                _s111(
                    0.35289110851016772428256675705798,
                    0.16210524214554394733483591857429
                    ),
                _s111(
                    0.40007606991958107069151584053398,
                    0.064525860099930314122766337266815
                    ),
                _s111(
                    0.18778434154653763297209499807307,
                    0.11437969546856907063157395975885
                    ),
                _s111(
                    0.20737693769034118226800699485912,
                    0.069571995852782021126280423446365
                    ),
                _s111(
                    0.10780383022637546988937726806675,
                    0.056034359700542548202663007444833
                    ),
                _s111(
                    0.14790146607312609862638050590311,
                    0.08158966578416703206891505452374
                    ),
                _s111(
                    0.16004266841135878852622500746535,
                    0.044629273771224669700718555567589
                    ),
                _s111(
                    0.37787666968877206563012487773812,
                    0.1087918765992173235808924647183
                    ),
                _s21(0.25269631485220673221675794136592),
                _s111(
                    0.039223369664808447656715193564727,
                    0.0032061485845241110097623221725454
                    ),
                _s111(
                    0.32430444507713815227341908124842,
                    0.082976548058919074176214822993003
                    ),
                _s21(0.21057712951480426709120641695334),
                _s21(0.0030575106285419621355963939637556),
                _s21(0.4032922222742628321359745151733),
                _s21(0.47675164680029289700450025370353),
                _s111(
                    0.35098484762250473274316967793719,
                    0.0035207350099882036762047091833696
                    ),
                _s21(0.064856245159702626212816239741417),
                _s21(0.16889637978296184143606030913156),
                _s3(),
                _s111(
                    0.42087082762569367816699589484148,
                    0.00059829472863675968600029811927026
                    ),
                _s21(0.096298034159157740435429809219389),
                _s111(
                    0.25265603856026470126032332920197,
                    0.09968371946019945101580182229927
                    ),
                _s21(0.49070156913020397839438373594675),
                _s21(0.49846457923659763589631351067494),
                _s111(
                    0.30167033356052759662788374519028,
                    0.13082391437504861842095186887726
                    ),
                _s21(0.43227197760680636353651654709095),
                _s21(0.13024160699020923749739572965558),
                _s21(0.45687173715432281105343195312627),
                _s111(
                    0.016061387966790736725359414721741,
                    0.0029955877400127188213258890486186
                    ),
                ])
            self.weights = numpy.concatenate([
                0.00037755236785356990322757708844313 * numpy.ones(6),
                0.0015516849698550359084385491295349 * numpy.ones(6),
                0.00086441456074659423647479386781886 * numpy.ones(6),
                0.00099718493187316731591746312990576 * numpy.ones(6),
                0.00040907703129804700764131593266948 * numpy.ones(6),
                0.0020436671600844307304135391154096 * numpy.ones(6),
                0.00041264433262762257309362927233909 * numpy.ones(6),
                0.0019146967543968340203976985872252 * numpy.ones(6),
                0.0011640031359899424702851939497132 * numpy.ones(6),
                0.00291867371065237459833053921265 * numpy.ones(6),
                0.0019423381723223563664420601603254 * numpy.ones(6),
                0.0020437160642329849773517378304998 * numpy.ones(6),
                0.0004561788584348895030045035043707 * numpy.ones(6),
                0.0011181526506985234067300543560843 * numpy.ones(6),
                0.0030364724116123745547307536555749 * numpy.ones(6),
                0.0019401856791114052131068187412846 * numpy.ones(6),
                0.0019149234843929558597017936842575 * numpy.ones(6),
                0.0006256066196748514371615149625715 * numpy.ones(3),
                0.0085945289589868128172207393838077 * numpy.ones(3),
                0.0085993323415620050084563487229761 * numpy.ones(3),
                0.0082321340562275838450012576646079 * numpy.ones(6),
                0.0075395204858701980611600587750656 * numpy.ones(6),
                0.0018551746248726809175949283354567 * numpy.ones(6),
                0.0044048213741680656634651166660365 * numpy.ones(6),
                0.0038072050061858979893524121622984 * numpy.ones(6),
                0.00096944457905454078100578257382062 * numpy.ones(6),
                0.0042240803334807804220246275356012 * numpy.ones(6),
                0.0066193255529404269738560673769695 * numpy.ones(6),
                0.0010014072888268371624462871187038 * numpy.ones(6),
                0.0075563110745214984193506996702024 * numpy.ones(6),
                0.0052987050439936187089717269724419 * numpy.ones(6),
                0.0055546897961284238957347067826876 * numpy.ones(6),
                0.0045059845609494478568113850477248 * numpy.ones(6),
                0.0031131351389202917001968041481465 * numpy.ones(6),
                0.0043898382162026902936883890450377 * numpy.ones(6),
                0.003238552372082920411687341695089 * numpy.ones(6),
                0.006574298111609529413017870416334 * numpy.ones(6),
                0.0080391355993214451993928844468303 * numpy.ones(3),
                0.00046193758037465469064280320768214 * numpy.ones(6),
                0.0056203298622696309354222215462745 * numpy.ones(6),
                0.007211823594908872263574588785191 * numpy.ones(3),
                0.00012296388513332560873334018113975 * numpy.ones(3),
                0.0080668049674519090876011140230359 * numpy.ones(3),
                0.0046789762284529159779004230718363 * numpy.ones(3),
                0.0012327695569297402309073597808901 * numpy.ones(6),
                0.0027172985719247752773272192160794 * numpy.ones(3),
                0.0062399420305303467044423373185357 * numpy.ones(3),
                0.0087854844528631139044953305211305 * numpy.ones(1),
                0.00044901023819614685494859515799402 * numpy.ones(6),
                0.0038909093796049571865225459857781 * numpy.ones(3),
                0.0056881832653567421434693475100953 * numpy.ones(6),
                0.0030289691062772661651305403011695 * numpy.ones(3),
                0.001205028062255519063349073804926 * numpy.ones(3),
                0.0067479088961222548109320676991965 * numpy.ones(6),
                0.0072248857349588951928103737993307 * numpy.ones(3),
                0.0051363917261404386280940486316984 * numpy.ones(3),
                0.0060844538754952708403763127659768 * numpy.ones(3),
                0.00027856062504863362656153077290914 * numpy.ones(6),
                ])
        elif index == 41:
            bary = numpy.concatenate([
                _s21(0.0018100550930816792181455727663449),
                _s21(0.31460775138387783999446595468217),
                _s21(0.37183245783943581904371856116814),
                _s111(
                    0.026772802208715575644116155958975,
                    0.03937930192202371429337033778205
                    ),
                _s111(
                    0.070874972726685938044887387930946,
                    0.021467386112313013444400067694118
                    ),
                _s111(
                    0.02276765890995240641335728720748,
                    0.01198511124275476302409355502637
                    ),
                _s111(
                    0.2127026585190979425575346305522,
                    0.096411202624268927933433027589243
                    ),
                _s111(
                    0.19707999494152021244839301290173,
                    0.12797822284643272491463268562628
                    ),
                _s21(0.48184896017583440656956857069247),
                _s111(
                    0.063913041436000726939862675044208,
                    0.042399106611156120901062962246614
                    ),
                _s111(
                    0.091625463275757714329820826597389,
                    0.011686997401174617539836947184562
                    ),
                _s111(
                    0.088917540355738020596058660917351,
                    0.063241532938583118630459956042788
                    ),
                _s111(
                    0.45923055890369903774129380073514,
                    0.014840944840779472668832879170775
                    ),
                _s111(
                    0.033884382234843877038407514292919,
                    0.0018400396156130491958817052266516
                    ),
                _s111(
                    0.32007521040661540801343551215441,
                    0.25805711373987160673744485661114
                    ),
                _s21(0.14228199830842082629649997053298),
                _s111(
                    0.42206849253188894655327590762935,
                    0.036219225434018616900113408281895
                    ),
                _s111(
                    0.047876441072928429446194896815604,
                    0.010543966044366349794044743549602
                    ),
                _s21(0.17463361486634168894048142237932),
                _s21(0.42519807841513996789079538486135),
                _s111(
                    0.36179902270607627469112983106453,
                    0.15102688749956549603687010799629
                    ),
                _s21(0.20997737233791780543938513470415),
                _s111(
                    0.06845694309222099282361257764115,
                    0.00248206090173231454121030108019
                    ),
                _s111(
                    0.23628847259645601691274740968919,
                    0.15495776702196281377766170014902
                    ),
                _s111(
                    0.29789541898504963895408305727251,
                    0.15196155366509877739190831498782
                    ),
                _s111(
                    0.11255725152392183531058671624445,
                    0.0023775457303720990464875952441714
                    ),
                _s111(
                    0.012065862971307616054991378548217,
                    0.0030492380858381186743187959303164
                    ),
                _s111(
                    0.39065309708631693165235500796633,
                    0.015075465750524257932222759281911
                    ),
                _s111(
                    0.42741072737397024098464451655922,
                    0.0028845541075684188323179736226094
                    ),
                _s21(0.26273997362534277134041787715106),
                _s111(
                    0.41184552209771978019242468902373,
                    0.1047651377258456064598216185206
                    ),
                _s21(0.49860363923111433028523304252744),
                _s111(
                    0.32222119451684366994821159286129,
                    0.014831767575572154518678439424394
                    ),
                _s111(
                    0.35452490163088131237412981985117,
                    0.036745193564627139867169905580208
                    ),
                _s111(
                    0.33202761729877971392087823315438,
                    0.20283297705625422215985233596356
                    ),
                _s111(
                    0.25302973908477228600902272760318,
                    0.065589317689696657506984221136736
                    ),
                _s111(
                    0.35703990489450974223197071789198,
                    0.0028334078818125141009749716253644
                    ),
                _s111(
                    0.28706068921165539638370722345368,
                    0.036052515766674376496147127543581
                    ),
                _s111(
                    0.13613135163113610660499195680144,
                    0.014156197127127483955427902774733
                    ),
                _s111(
                    0.1322241256617750545749789140314,
                    0.063772809672132732885343626183977
                    ),
                _s111(
                    0.10792421983955239460803835402987,
                    0.033553439545873064447476138995779
                    ),
                _s111(
                    0.15288690214138662149423598610157,
                    0.099579044603542747135230600108531
                    ),
                _s111(
                    0.18810177192960669125909479431962,
                    0.064649152025039321218759961823047
                    ),
                _s111(
                    0.25498242915088791756667013822044,
                    0.01489872687080796084270260762596
                    ),
                _s111(
                    0.22307808048127770528286468345289,
                    0.0028544472780138035160792464495396
                    ),
                _s21(0.39956652552697029660921180945476),
                _s111(
                    0.16405073602103255292130775565695,
                    0.0028307718211817070212425007931734
                    ),
                _s111(
                    0.28815111741898860859509676194938,
                    0.0028359467965889003515083365392174
                    ),
                _s111(
                    0.34095778887406159305735848265559,
                    0.10563268627817783184799001072941
                    ),
                _s111(
                    0.26829344806928391596394025647055,
                    0.20540889578274667850556814492906
                    ),
                _s111(
                    0.31940900352832468564698311379966,
                    0.067478631067926277748283767786637
                    ),
                _s111(
                    0.19229546518595229189086788810084,
                    0.014933107362951521295352411944572
                    ),
                _s111(
                    0.16041405340292364698238341662474,
                    0.035385293405442621831179434142974
                    ),
                _s111(
                    0.22125349898042174805121008150062,
                    0.036251995833857701246656630121245
                    ),
                _s111(
                    0.27119422406881354565727817057382,
                    0.10659969697303733210308559789584
                    ),
                _s111(
                    0.39232341129977041228543151656179,
                    0.066569707653674558153759730191275
                    ),
                _s21(0.46676105778551496920619985907799),
                _s21(0.10309548557603061178603678857082),
                ])
            self.weights = numpy.concatenate([
                0.000052209323125515952483256062203907 * numpy.ones(3),
                0.0064432696708827039027001878342356 * numpy.ones(3),
                0.0056462037725046143083090298366197 * numpy.ones(3),
                0.00094663399645289961740844892598064 * numpy.ones(6),
                0.0010778383293870836532404670103395 * numpy.ones(6),
                0.00055355140630463948775668212061983 * numpy.ones(6),
                0.0037447209152799688166870553502542 * numpy.ones(6),
                0.0040737392096223579909862219509848 * numpy.ones(6),
                0.0029751068247423738464067057240693 * numpy.ones(3),
                0.0017887545706646271701946597290686 * numpy.ones(6),
                0.0010314247381434697657830895685431 * numpy.ones(6),
                0.0025322072683270647353267872656321 * numpy.ones(6),
                0.0022524785579246364796122022411897 * numpy.ones(6),
                0.00028183750846594049371281096262836 * numpy.ones(6),
                0.0061093730696482188239615205218687 * numpy.ones(6),
                0.0042449557771111947536413747443463 * numpy.ones(3),
                0.0033089712153023544714697695851751 * numpy.ones(6),
                0.00075034037410059840067312214669489 * numpy.ones(6),
                0.005036251708976386900352708095373 * numpy.ones(3),
                0.0060101004771622411515890794037109 * numpy.ones(3),
                0.0062236592332095752419269682317926 * numpy.ones(6),
                0.0060772289969773640516985685867224 * numpy.ones(3),
                0.00049042087101432951192021573950062 * numpy.ones(6),
                0.0058906442657042353156260464677643 * numpy.ones(6),
                0.0062381201191140057762323842046705 * numpy.ones(6),
                0.00058891818812813349795866540217221 * numpy.ones(6),
                0.00023905156180476898675590801093811 * numpy.ones(6),
                0.0023381337598691202160321261398357 * numpy.ones(6),
                0.0010454329707756071968878431767348 * numpy.ones(6),
                0.0067671402084224723435638101357475 * numpy.ones(3),
                0.0059901891203095585772322520418084 * numpy.ones(6),
                0.0010236399703004335949955100569089 * numpy.ones(3),
                0.002288249587261252373016857515915 * numpy.ones(6),
                0.0036104870565558439466143940368033 * numpy.ones(6),
                0.007058657808550164555829640788964 * numpy.ones(6),
                0.0043785525737027569775286098477022 * numpy.ones(6),
                0.0010171383589503082953007433274045 * numpy.ones(6),
                0.0034554845247204603883551585328093 * numpy.ones(6),
                0.0016626192415771936180953861569058 * numpy.ones(6),
                0.0033511706532686987985376899965406 * numpy.ones(6),
                0.0023865646505328902212665456813466 * numpy.ones(6),
                0.0041753789931688428890888799877875 * numpy.ones(6),
                0.0038711214885426366049428479514011 * numpy.ones(6),
                0.002192339196043895014093104371291 * numpy.ones(6),
                0.00091273233652152079140230882265003 * numpy.ones(6),
                0.0072301097329516989158929252752935 * numpy.ones(3),
                0.00080283465261535068488709460001562 * numpy.ones(6),
                0.00097935960157862393091604916689658 * numpy.ones(6),
                0.0059550867608880760747265675231062 * numpy.ones(6),
                0.0067525153387175384632397914219197 * numpy.ones(6),
                0.0049326746082458918989542182770734 * numpy.ones(6),
                0.0019995171236440128210194459506808 * numpy.ones(6),
                0.0028581698928275127063633520021953 * numpy.ones(6),
                0.0032573829494494207517022705435398 * numpy.ones(6),
                0.0058722816583725265815057888666697 * numpy.ones(6),
                0.0050950784475005244247562375702472 * numpy.ones(6),
                0.0051338703252384472502820906971867 * numpy.ones(3),
                0.0039695690393596122823616274127965 * numpy.ones(3),
                ])
        elif index == 42:
            bary = numpy.concatenate([
                _s21(0.40237192510979328815564123689329),
                _s21(0.0013653925759541381426100407060897),
                _s111(
                    0.3305538536843532134263235690722,
                    0.001570425601462847873574198366295
                    ),
                _s21(0.49753172962520537269975877776702),
                _s111(
                    0.0096146691207101573205359143539356,
                    0.020431365732606569781308396513359
                    ),
                _s111(
                    0.057198306193975713708869684231522,
                    0.021903558615629309688509478820286
                    ),
                _s21(0.47990029575441913924092321658599),
                _s21(0.066466329055783719025703792704347),
                _s111(
                    0.45327097163794561466065258031135,
                    0.0014935512902909570214300021257703
                    ),
                _s111(
                    0.14073247390672036599019346074133,
                    0.055441499674364493180723755361478
                    ),
                _s111(
                    0.18378863030194105218400633851071,
                    0.1186935694213949307491677527776
                    ),
                _s21(0.043236389387746108758580549676078),
                _s111(
                    0.032112152486877870105143811172573,
                    0.0013930572171860837794157066353131
                    ),
                _s21(0.49030501420164670223736222249349),
                _s111(
                    0.042964072200676188480198294930469,
                    0.010299475181841361840030452273162
                    ),
                _s111(
                    0.010524842730074891544603908899763,
                    0.0025847099922881149693025519008533
                    ),
                _s111(
                    0.075379028562970311311420420528027,
                    0.039474312716157624670606215465357
                    ),
                _s111(
                    0.43633336208995605640374595289387,
                    0.012293838524552470152345280442423
                    ),
                _s111(
                    0.27049074651497698397818288760795,
                    0.0025898092478456274950087851499475
                    ),
                _s111(
                    0.3508136433810973700272284665632,
                    0.037556388872296805218018323300275
                    ),
                _s111(
                    0.37276551197429214988873521636513,
                    0.017094546720707277349519073626461
                    ),
                _s111(
                    0.21710735226334155451349426697109,
                    0.085758316294315819410561211201739
                    ),
                _s111(
                    0.19275636500900589633933511988896,
                    0.057985644964114331941519066433975
                    ),
                _s111(
                    0.42068503068641629172573660198072,
                    0.033848172882990274420784111796429
                    ),
                _s21(0.026998822628373124822688189712742),
                _s111(
                    0.21248518543709112026967537408026,
                    0.0025266913852085761097286996309678
                    ),
                _s111(
                    0.31925720872918739221459975885537,
                    0.010656542115527591520196347660107
                    ),
                _s21(0.13850480691032286285507760238703),
                _s21(0.46701320735576399397406257532003),
                _s111(
                    0.25163823551257262023947883458939,
                    0.013429673242258488430198923808594
                    ),
                _s111(
                    0.15221198145680938529521660950749,
                    0.088341778293265027519406148264097
                    ),
                _s111(
                    0.11114209479403333159964354017892,
                    0.032922118211160733047066687298779
                    ),
                _s111(
                    0.35815005885125032581791165651246,
                    0.19480076077732762948298405029594
                    ),
                _s111(
                    0.16379985695615011086690341502245,
                    0.031461015961891729830824166578061
                    ),
                _s111(
                    0.10211277333667435979437358591607,
                    0.067368310164525829869405494069245
                    ),
                _s111(
                    0.23355570301935160611932101765801,
                    0.18656524419997268713146365699031
                    ),
                _s111(
                    0.32141646283171138605797714057966,
                    0.063883483412395297676296902982187
                    ),
                _s111(
                    0.22333618489571259724178771286663,
                    0.032571221794949826162421456834491
                    ),
                _s111(
                    0.1905881189463203289326559374577,
                    0.013273482756157748726400857965572
                    ),
                _s111(
                    0.22735584024547954583131247812384,
                    0.13623377550912803400124005560337
                    ),
                _s21(0.1065804258588589073952758376067),
                _s111(
                    0.15676792278716006785214418181128,
                    0.0026153775338891857398008408451319
                    ),
                _s111(
                    0.13403333498677852517002838828898,
                    0.013481211205156754370064445274221
                    ),
                _s111(
                    0.29572117574545365623543149677162,
                    0.19188410435220618881731445716081
                    ),
                _s111(
                    0.10671927181300039247424339126155,
                    0.0025472816321698190165377194741646
                    ),
                _s111(
                    0.41256641948391657654544017789056,
                    0.10103729435575962439479627155977
                    ),
                _s111(
                    0.30785198299288609277563355754385,
                    0.24975074116528870465509617756079
                    ),
                _s111(
                    0.087020548699696697565120854127083,
                    0.013712625288557960069518089273335
                    ),
                _s111(
                    0.3894283316812197198976770727437,
                    0.003344595433402609966871760925933
                    ),
                _s21(0.31173539107258736916939551300132),
                _s111(
                    0.26129397031417700748814310809243,
                    0.055674630430605156666917334533113
                    ),
                _s111(
                    0.064975533109811901240670247256292,
                    0.0028037585935402048179489930985948
                    ),
                _s111(
                    0.27195000609813104592228980821749,
                    0.096350500077239797486713192868162
                    ),
                _s111(
                    0.29309036427254626210251966651541,
                    0.029199106600435966523786517748782
                    ),
                _s111(
                    0.39466482658781922132544916427517,
                    0.063574550662169118573832683387381
                    ),
                _s111(
                    0.35685221819298281597534415674046,
                    0.14392858006554380500837719829602
                    ),
                _s21(0.24486057586888178570860223935019),
                _s111(
                    0.34136135163887632082342626955661,
                    0.099208699753299723465361088508405
                    ),
                _s111(
                    0.2895129410795146086089159514496,
                    0.14038892142146281599242417081344
                    ),
                _s21(0.17483514581119072083932072414441),
                _s21(0.42756933415019040744201499721584),
                _s21(0.37427015667844769299026390146509),
                ])
            self.weights = numpy.concatenate([
                0.0035464400159771666962340527612761 * numpy.ones(3),
                0.000034719535436562089078710043551835 * numpy.ones(3),
                0.00053273169434053232025913630532291 * numpy.ones(6),
                0.0010141696428327357720574767790149 * numpy.ones(3),
                0.00043907702716982386267737757587166 * numpy.ones(6),
                0.0010370250465354264110580982882377 * numpy.ones(6),
                0.002844185957012238834469743944398 * numpy.ones(3),
                0.0018709815380228792571796815672242 * numpy.ones(3),
                0.00055955195022804377463706923837973 * numpy.ones(6),
                0.0025379205064647718257098939105595 * numpy.ones(6),
                0.0040581577981139574186977310197752 * numpy.ones(6),
                0.0013326250004079803556043883111094 * numpy.ones(3),
                0.00022664257520749128542584061021427 * numpy.ones(6),
                0.0021926816815858840991121779728007 * numpy.ones(3),
                0.00068185334314149642356658250867471 * numpy.ones(6),
                0.00018709306861560231959460410781152 * numpy.ones(6),
                0.0016866296698237316344446797650937 * numpy.ones(6),
                0.0019025497964578995081075055427877 * numpy.ones(6),
                0.00077336702529978049998484955485256 * numpy.ones(6),
                0.0031723878909002393211316290180665 * numpy.ones(6),
                0.002168524709861741656309169235359 * numpy.ones(6),
                0.0040741174736407135345212505271893 * numpy.ones(6),
                0.0033720890143961500658627262626054 * numpy.ones(6),
                0.0033374654153493006867868974202338 * numpy.ones(6),
                0.00095240595565468112173347014152002 * numpy.ones(3),
                0.00074430046071222707533639850781145 * numpy.ones(6),
                0.0017490056575092856048520491056517 * numpy.ones(6),
                0.004359850041892466468959842192335 * numpy.ones(3),
                0.0045325322550327508176216624811282 * numpy.ones(3),
                0.0019031683954380265016199786214044 * numpy.ones(6),
                0.0039581444907565927266012672532042 * numpy.ones(6),
                0.0021633300119175084266112907364381 * numpy.ones(6),
                0.0062667413965418204518070839014896 * numpy.ones(6),
                0.0024314709447528086539189694938141 * numpy.ones(6),
                0.0028265780976024542965910120800526 * numpy.ones(6),
                0.0064753562254832802873748235610856 * numpy.ones(6),
                0.0042484824512541860324387865699311 * numpy.ones(6),
                0.0029100859090622533835372303505772 * numpy.ones(6),
                0.0017543943309314071019005437626674 * numpy.ones(6),
                0.0053278143635325773115116557154912 * numpy.ones(6),
                0.0036960706615538763809116442504715 * numpy.ones(3),
                0.0007117355150533799892535900669087 * numpy.ones(6),
                0.0015685196181434368933770178779162 * numpy.ones(6),
                0.006909683703385332649004074427308 * numpy.ones(6),
                0.00060776688933841591212082329611747 * numpy.ones(6),
                0.0058599060505518001238729707720451 * numpy.ones(6),
                0.0077651860406296855657295157435621 * numpy.ones(6),
                0.001310040832986112085429096612836 * numpy.ones(6),
                0.001057097249172270163098249007835 * numpy.ones(6),
                0.0080973376090283895894316300654385 * numpy.ones(3),
                0.0040411976507506002926653039064974 * numpy.ones(6),
                0.00052424789134593807477469885768415 * numpy.ones(6),
                0.0052048680063201082803127583222162 * numpy.ones(6),
                0.0030333113262698715872816634393886 * numpy.ones(6),
                0.0048685523603978672544401311971399 * numpy.ones(6),
                0.0066028726180889637791397220147007 * numpy.ones(6),
                0.0072906010527867706385803535012068 * numpy.ones(3),
                0.0056531366883124027492484104256896 * numpy.ones(6),
                0.0061919290123439353712716412317378 * numpy.ones(6),
                0.0060555267820253706088134821421024 * numpy.ones(3),
                0.0067589446392297668701956874367544 * numpy.ones(3),
                0.007922044576591311385557734242534 * numpy.ones(3),
                ])
        elif index == 43:
            bary = numpy.concatenate([
                _s21(0.41178017085993385194101544969374),
                _s21(0.47557586709913148045619964597896),
                _s111(
                    0.22973943652724922221090629196654,
                    0.14545379538783254545396788785531
                    ),
                _s21(0.49917983673447802794571526086407),
                _s21(0.2962479881688698839478518356316),
                _s111(
                    0.36055807608878384924419210726749,
                    0.03241547484177850137887961130067
                    ),
                _s111(
                    0.095917227666943707613218067790621,
                    0.012412908224438337684287476053658
                    ),
                _s111(
                    0.18980151011108849289710692322325,
                    0.13785811789149299648537068528566
                    ),
                _s111(
                    0.40105845129260520468872264478266,
                    0.13533613115359614053833924495954
                    ),
                _s111(
                    0.27758253227610216648171214839283,
                    0.14352679189247660182688695039205
                    ),
                _s21(0.35318880488700073407544177436172),
                _s111(
                    0.306118376690334997880101898225,
                    0.03545159159651415070444373622507
                    ),
                _s21(0.49510598080276992333134259212415),
                _s111(
                    0.41828917017342478368947218205061,
                    0.035418533131322089677876843292273
                    ),
                _s21(0.03290981596053466761522804037944),
                _s111(
                    0.09650992433575238017682440983043,
                    0.055793868125798097489459219489449
                    ),
                _s111(
                    0.13283974111052923119398570119816,
                    0.063976234789808573618274556938945
                    ),
                _s111(
                    0.096964329798066351729024215526795,
                    0.030262579741771051575841825146069
                    ),
                _s21(0.48656375113068838174730293992424),
                _s111(
                    0.42875375369101297392393518667051,
                    0.014637616358968465956529385442583
                    ),
                _s21(0.062715048389761751504510024978092),
                _s111(
                    0.13685344532608310212524510475413,
                    0.013605250929398435154680514521313
                    ),
                _s111(
                    0.36356190909925285001928716317685,
                    0.0024760406221922225597492674986809
                    ),
                _s21(0.14474355977022148471601196730564),
                _s111(
                    0.33614357232761974758391625777852,
                    0.14292597337792934638376261887627
                    ),
                _s21(0.013364620608261697555512657840211),
                _s111(
                    0.09409885344990038057854875529763,
                    0.0024207777457972175421325294672996
                    ),
                _s111(
                    0.24271058666990053730942994625812,
                    0.1951503399895133076621945816127
                    ),
                _s111(
                    0.42523844283105769523604957109931,
                    0.064055101821656949489484902160386
                    ),
                _s111(
                    0.36427812697418732263553687983246,
                    0.013187611232483199816750862866052
                    ),
                _s111(
                    0.1776637264896365315742277317256,
                    0.061873082107012595735479193469917
                    ),
                _s111(
                    0.058845647949205092546865131135913,
                    0.0026124723427687684156060007403035
                    ),
                _s111(
                    0.25009960190982928387181471553683,
                    0.034190403560094950840876833089042
                    ),
                _s111(
                    0.031972586536534605331984573137376,
                    0.0026508313305474216926497167339774
                    ),
                _s111(
                    0.24172396135767907543636150834606,
                    0.0025915690897239313464661653860411
                    ),
                _s111(
                    0.18440470743920183100146783911369,
                    0.0142733464909807001170073735571
                    ),
                _s111(
                    0.061256389015271970187860622652456,
                    0.013856407685434492282208085099867
                    ),
                _s111(
                    0.013257141692271871088963444855027,
                    0.0024938262659978295231806520804869
                    ),
                _s21(0.0025649413700202923342105014440312),
                _s21(0.18777224277857851138689496978675),
                _s111(
                    0.13797314008145525893325690760939,
                    0.10044093105837189605760717403183
                    ),
                _s111(
                    0.033205679927311351594434949459772,
                    0.013697442352723626950898102329667
                    ),
                _s21(0.38256430878693553365233886877071),
                _s111(
                    0.19446619678110650131794975775816,
                    0.034143195889543975771021109245621
                    ),
                _s111(
                    0.23948263508520458147502438043399,
                    0.013804625294641210668899181123423
                    ),
                _s111(
                    0.36455422777634081033278312907258,
                    0.18649721429955501945511238429315
                    ),
                _s111(
                    0.13735640236619982748077346739604,
                    0.0025417593646837136406903623189061
                    ),
                _s111(
                    0.23223203731889045506527533067146,
                    0.063034933704680926479464720763072
                    ),
                _s111(
                    0.29266139507224689160354569239305,
                    0.06387691257936595102444230157325
                    ),
                _s111(
                    0.3176481163981567552579505887103,
                    0.24195284517547509062590029014886
                    ),
                _s111(
                    0.3003015966157783605350142528519,
                    0.014599271258976278094695144611404
                    ),
                _s111(
                    0.30056135371681098474665922234466,
                    0.00281811459692606506855326547554
                    ),
                _s111(
                    0.060523979398731234903411313200646,
                    0.034134499559379701912954550822417
                    ),
                _s111(
                    0.30078591923388939455179548891259,
                    0.19335547727641827434387718576435
                    ),
                _s111(
                    0.18955070397768522067973334643511,
                    0.097362969750150811796044244576555
                    ),
                _s111(
                    0.43013990336214889610538880435087,
                    0.0027883024051113908100410686249799
                    ),
                _s111(
                    0.18701857221766473506056737817309,
                    0.0027559624464836448093851788488073
                    ),
                _s111(
                    0.37918952250900750016756903598791,
                    0.097827107158867704972767657699656
                    ),
                _s111(
                    0.31296697698215932056546390740617,
                    0.099573685636870894197442144055445
                    ),
                _s111(
                    0.14180448701728903007288183143714,
                    0.0341143516662648087521958455595
                    ),
                _s21(0.45289226084595103617095678858064),
                _s111(
                    0.24818146001558851349936327969455,
                    0.10049744327292836435046582692688
                    ),
                _s111(
                    0.35879783540799872223411343110695,
                    0.061673803696551649468581772363764
                    ),
                _s21(0.092867046214047944518085778729407),
                _s21(0.25248493264909078369043974455353),
                ])
            self.weights = numpy.concatenate([
                0.0030799655691816325836512202136952 * numpy.ones(3),
                0.0024873482202464832365503711196722 * numpy.ones(3),
                0.0039126383408175682585199879403685 * numpy.ones(6),
                0.00060861572218077991162264798240214 * numpy.ones(3),
                0.0058780823337591549568280953907767 * numpy.ones(3),
                0.0026575597720771203532046312234442 * numpy.ones(6),
                0.0010403858576163626226810643389964 * numpy.ones(6),
                0.0038688222455241618504947323194077 * numpy.ones(6),
                0.0054360079215850587366658320008417 * numpy.ones(6),
                0.0051330210770134501527915183860799 * numpy.ones(6),
                0.0068684889665149535177586185661464 * numpy.ones(3),
                0.0027248852604064804038373862266246 * numpy.ones(6),
                0.001614307243915753657132498173042 * numpy.ones(3),
                0.0030819278277026026453793694837109 * numpy.ones(6),
                0.0010836224292710483115699663619806 * numpy.ones(3),
                0.0021960146608965101905778121175494 * numpy.ones(6),
                0.0027071716291894682915515012085397 * numpy.ones(6),
                0.0017732505663029985834679947087491 * numpy.ones(6),
                0.0027166050682951765829687069096522 * numpy.ones(3),
                0.0021194095320855991250844253229391 * numpy.ones(6),
                0.0020243293045496533470321835592698 * numpy.ones(3),
                0.0013822658273996519320236146275469 * numpy.ones(6),
                0.00082877960227598081370135516171246 * numpy.ones(6),
                0.0042454319210643227829964039065619 * numpy.ones(3),
                0.005760929970368311613302276564644 * numpy.ones(6),
                0.00047210210557920691904737753088959 * numpy.ones(3),
                0.00048664065357893428148800814810624 * numpy.ones(6),
                0.006061170604496654301054177666805 * numpy.ones(6),
                0.0039591603530013185673838333474277 * numpy.ones(6),
                0.0019515305614847015748449281479838 * numpy.ones(6),
                0.0031758640195784995248382584114296 * numpy.ones(6),
                0.00041776897840248034528108419326072 * numpy.ones(6),
                0.0028037985130037499620381517085299 * numpy.ones(6),
                0.00030743230012411995451602896493 * numpy.ones(6),
                0.00075804787451567809657090193569572 * numpy.ones(6),
                0.0016285601871914179347760113093859 * numpy.ones(6),
                0.0010015815151316952892485875189395 * numpy.ones(6),
                0.00018970350922109701593115139504011 * numpy.ones(6),
                0.000085732576679882198277005042682718 * numpy.ones(3),
                0.0053352753423961882197694059355971 * numpy.ones(3),
                0.0039225937188264304439200578005384 * numpy.ones(6),
                0.00073526444907111327075429780517715 * numpy.ones(6),
                0.0069496808690565782301338954617313 * numpy.ones(3),
                0.0026277500982665875036738449530496 * numpy.ones(6),
                0.0018416339190279123874687232091982 * numpy.ones(6),
                0.0059250327726737599568805600947878 * numpy.ones(6),
                0.00061233884081429031242706895284032 * numpy.ones(6),
                0.0038249634484788579580887840314782 * numpy.ones(6),
                0.0041283032222062560392852347731552 * numpy.ones(6),
                0.006692032041741229488106304357568 * numpy.ones(6),
                0.0020472829473929802328233366449517 * numpy.ones(6),
                0.00087801567916509065616995276674699 * numpy.ones(6),
                0.0015777851185458777046731303402116 * numpy.ones(6),
                0.0063074359109884056568133487334034 * numpy.ones(6),
                0.0042726398210531675072497002183799 * numpy.ones(6),
                0.00096983734623925983388711564713531 * numpy.ones(6),
                0.00073801827045786280192817844778589 * numpy.ones(6),
                0.0052618813601275176819778877273431 * numpy.ones(6),
                0.0053002719749982506825567236892445 * numpy.ones(6),
                0.0024558010102791692520256530213621 * numpy.ones(6),
                0.0053328054204010124725411162517315 * numpy.ones(3),
                0.0050556103951576838639219579001529 * numpy.ones(6),
                0.0044341909691302104883505729707572 * numpy.ones(6),
                0.0032907061458539955765032762707218 * numpy.ones(3),
                0.0073142091431203385404764277288737 * numpy.ones(3),
                ])
        elif index == 44:
            bary = numpy.concatenate([
                _s21(0.0017123075799638136895121702696589),
                _s111(
                    0.01975516761571974351422857470556,
                    0.010756719556684273996841403575671
                    ),
                _s111(
                    0.099516581381844660588237381611657,
                    0.00078058589986187171439914503388629
                    ),
                _s111(
                    0.15349454475521389848657446586064,
                    0.0019383329127905529842828651823285
                    ),
                _s111(
                    0.46060745200399893797476184757364,
                    0.0017861424805253371391154951570748
                    ),
                _s111(
                    0.029642060426444401449580575257606,
                    0.0017542251710527286582691333358851
                    ),
                _s111(
                    0.33384062482481181576210433921913,
                    0.00051729425966005439029876820081639
                    ),
                _s21(0.029167351465300611907976357415034),
                _s21(0.049549987912302506132814369017368),
                _s111(
                    0.11942668722741197261052886165152,
                    0.0075458709257060963506816857266597
                    ),
                _s111(
                    0.17331180669942895680964814481991,
                    0.011492546452805203327582669782099
                    ),
                _s111(
                    0.059397863385756653388396168102832,
                    0.024262967721104459836850617776141
                    ),
                _s111(
                    0.25944616462784317093231061266264,
                    0.0097976933024595940058660502442421
                    ),
                _s111(
                    0.26966487953111786089778648405519,
                    0.0018991508282082745299648290445182
                    ),
                _s111(
                    0.32342841041095645111497221201332,
                    0.0076746798561443045582585746398954
                    ),
                _s111(
                    0.20939760943519251101697756688524,
                    0.002707204307768565381246592814512
                    ),
                _s111(
                    0.010793737386739596633959626228363,
                    0.0027063965571698294841264945878515
                    ),
                _s111(
                    0.0416577600266789433427828409968,
                    0.010462734029639450895873172584785
                    ),
                _s111(
                    0.13930077441848229029512128936049,
                    0.021512860284919813351363336815925
                    ),
                _s21(0.4877239805831533674208505234338),
                _s111(
                    0.097140757455529047963043568865105,
                    0.022824535634951918878928858350252
                    ),
                _s111(
                    0.15112480973897059667318232965728,
                    0.10008853030478997846330041884439
                    ),
                _s21(0.49606180446547498931665457237687),
                _s111(
                    0.39393239553515449151936850417343,
                    0.0031595695366354523881048353649248
                    ),
                _s111(
                    0.12763745941732481313566103373729,
                    0.043567575280835105567536532145024
                    ),
                _s111(
                    0.21534278591705231735017346284812,
                    0.062815376245839494031917519812089
                    ),
                _s111(
                    0.18052388725602852975372800064306,
                    0.034691186616931116096078230640744
                    ),
                _s111(
                    0.37069747874972836853794096662522,
                    0.01595213576803113731312706513902
                    ),
                _s111(
                    0.23890434663459388544299879303359,
                    0.036565690980963945529753922700004
                    ),
                _s21(0.076596115605173873480894452985672),
                _s111(
                    0.40354450476491695048804616794253,
                    0.12286724840697114637357355671868
                    ),
                _s111(
                    0.25108861840106174498181275042029,
                    0.12206119203351075502150635066201
                    ),
                _s111(
                    0.31009246035300195724524354098638,
                    0.087827887032806748161559833919124
                    ),
                _s111(
                    0.41718893706811110468755886275192,
                    0.065681103056196455268309099965051
                    ),
                _s21(0.47595184605277128325997658069699),
                _s111(
                    0.084277197515868819848966112651716,
                    0.045720147184830505864113188027838
                    ),
                _s111(
                    0.28751728189049522759248548922039,
                    0.052803410895066976988859462268692
                    ),
                _s21(0.10909554195723774092741869605563),
                _s111(
                    0.34214218490699248634063749981477,
                    0.033262707550278518392991900900333
                    ),
                _s21(0.23124700369384099426412870400944),
                _s21(0.19590819822365059223003202553593),
                _s111(
                    0.29078133562748429923987231757242,
                    0.023090951165029763442807189827358
                    ),
                _s111(
                    0.077738529716100318558812353439734,
                    0.0084121228886457889554283546116549
                    ),
                _s111(
                    0.058353808798156726797388294443294,
                    0.0018558446904390173206187750733971
                    ),
                _s111(
                    0.35373621739118744807194779737164,
                    0.15186122149915331840223013650024
                    ),
                _s111(
                    0.19410368342192220550142277421835,
                    0.14669384649550368818431007870285
                    ),
                _s111(
                    0.11524468026207428787853510967263,
                    0.073058272423484747241113693777017
                    ),
                _s111(
                    0.25264891885682738685779437981035,
                    0.0833205145800466186582346257643
                    ),
                _s111(
                    0.37063812791103511671414892535759,
                    0.097670099519353697704984791383687
                    ),
                _s111(
                    0.19751076623532845040715550590117,
                    0.10212747323582017620463487878641
                    ),
                _s111(
                    0.21910217877663961150173681681716,
                    0.016940122646793473751474498249445
                    ),
                _s111(
                    0.1647494714555341753621342963377,
                    0.064607414654394189883030785404828
                    ),
                _s111(
                    0.34875688743061502363855680277301,
                    0.19429349186192532111681887418138
                    ),
                _s111(
                    0.355853330058662296129275069862,
                    0.058789587775321184870543589871329
                    ),
                _s111(
                    0.34871652898546943500991263699605,
                    0.24358497588517304325810093756688
                    ),
                _s21(0.4223218082071513691128748374919),
                _s21(0.4548817499595059604862078968651),
                _s111(
                    0.29070433554969610797160517312839,
                    0.18612596595616319976382359845807
                    ),
                _s21(0.1439001195404209205099721562311),
                _s111(
                    0.28986293481711371191986298581831,
                    0.23835701955255615001217202692543
                    ),
                _s21(0.40234044785473540134486965459046),
                _s111(
                    0.30938357544545799719768656032141,
                    0.13122961679197293893392872399713
                    ),
                _s111(
                    0.24697953909924447234886857340319,
                    0.1684589785524791567919746199932
                    ),
                _s111(
                    0.4353451951214652447864463294682,
                    0.013125469511223459830734779315977
                    ),
                _s21(0.29429576569083396443183644443979),
                _s21(0.35161578879576900186639188382097),
                _s111(
                    0.41370128931488312661055209723906,
                    0.034413440518580169077333504320989
                    ),
                ])
            self.weights = numpy.concatenate([
                0.000044304090914147461083608637148873 * numpy.ones(3),
                0.00041197564643695758682848689584431 * numpy.ones(6),
                0.00027324084600368188573579137691349 * numpy.ones(6),
                0.00053932844392989668164092111652383 * numpy.ones(6),
                0.00067097897058203822780354221256064 * numpy.ones(6),
                0.00022593153611393980724212417823821 * numpy.ones(6),
                0.00032206630941833743640016604738071 * numpy.ones(6),
                0.0011229773828237150111127423108677 * numpy.ones(3),
                0.0017451515655329449025228389316994 * numpy.ones(3),
                0.0009677673027951324028424826938538 * numpy.ones(6),
                0.0013466648921817605344666355515376 * numpy.ones(6),
                0.0013375082839812344984750556021468 * numpy.ones(6),
                0.0012094280366882308935155238788157 * numpy.ones(6),
                0.00059737614176188248114610324026525 * numpy.ones(6),
                0.0014682010637335697046411641960392 * numpy.ones(6),
                0.00080105448520124065354435012027682 * numpy.ones(6),
                0.0001855670423377325419005065519135 * numpy.ones(6),
                0.00072142202642612238530518988481701 * numpy.ones(6),
                0.0015321667240474278548336636590281 * numpy.ones(6),
                0.0027296527164899309709901277687046 * numpy.ones(3),
                0.0015090266204282936697273898386619 * numpy.ones(6),
                0.0034017224729490143862421418498012 * numpy.ones(6),
                0.00120450601341192323812630226734 * numpy.ones(3),
                0.0010359881638063848695363565725658 * numpy.ones(6),
                0.0023774876927798601099626944226953 * numpy.ones(6),
                0.0033953912341846749022349746696463 * numpy.ones(6),
                0.002704088382481818220007165186299 * numpy.ones(6),
                0.002017003390291287511535411262029 * numpy.ones(6),
                0.0027918811324257680474094096828913 * numpy.ones(6),
                0.0023657933107792473962125110896766 * numpy.ones(3),
                0.0043004588657304928988602122781004 * numpy.ones(6),
                0.0049838565291937565193527613575372 * numpy.ones(6),
                0.0046532962756669859796994360494392 * numpy.ones(6),
                0.0041278760422565010490154961789407 * numpy.ones(6),
                0.0035111942683412106167749241573608 * numpy.ones(3),
                0.0020578217623524905648156356616485 * numpy.ones(6),
                0.0038938204584288090324569734022558 * numpy.ones(6),
                0.002943468047195771546791754949162 * numpy.ones(3),
                0.0027707613965555407859412233583193 * numpy.ones(6),
                0.0067232428984730884936431030275869 * numpy.ones(3),
                0.004979858812518916310441457647736 * numpy.ones(3),
                0.0024801439985067616973064068500625 * numpy.ones(6),
                0.0008629021808335470432996124145408 * numpy.ones(6),
                0.00032093317757429485661872857295845 * numpy.ones(6),
                0.0049810086106174347098616192391821 * numpy.ones(6),
                0.0050195779180815229311867378599581 * numpy.ones(6),
                0.0027692328182935268195102268270983 * numpy.ones(6),
                0.0037015806596788268358161313559895 * numpy.ones(6),
                0.0050750374591619946093164799727134 * numpy.ones(6),
                0.0042676766915481803290805204452214 * numpy.ones(6),
                0.0017292273039881342495257366948476 * numpy.ones(6),
                0.0032247399122439484996999125300739 * numpy.ones(6),
                0.005287947712961897016355002830123 * numpy.ones(6),
                0.0040195842863650699159827576828669 * numpy.ones(6),
                0.0060608332073000481632600957576789 * numpy.ones(6),
                0.0049208372452467634789881377457718 * numpy.ones(3),
                0.0050214750285975961653994250803055 * numpy.ones(3),
                0.0057991370585632813923725957419318 * numpy.ones(6),
                0.0044465932407613083513924289332375 * numpy.ones(3),
                0.0063931475211100405316407424617316 * numpy.ones(6),
                0.0046466127645843993821018056947583 * numpy.ones(3),
                0.005654675008988430625533452205662 * numpy.ones(6),
                0.0053110494168999697079018267551565 * numpy.ones(6),
                0.0019732304214758985370024422809498 * numpy.ones(6),
                0.0064254720501743216199677359046688 * numpy.ones(3),
                0.0062655878836268250168776102555506 * numpy.ones(3),
                0.0035554794715669390910633920101491 * numpy.ones(6),
                ])
        elif index == 45:
            bary = numpy.concatenate([
                _s111(
                    0.01260480383073618169352002314306,
                    0.0072924295508433190860092162440195
                    ),
                _s111(
                    0.004605476631441322913693083038019,
                    0.0013499971821104595021124942421118
                    ),
                _s111(
                    0.067946496603911169897812327768347,
                    0.006859175938060519703313967658085
                    ),
                _s111(
                    0.069338004363440864040157454564223,
                    0.0012297679959348406009902890172741
                    ),
                _s111(
                    0.46827085987094802197116969398141,
                    0.008780100200694679959272001187448
                    ),
                _s111(
                    0.25049517033522151635779109706826,
                    0.0016830263611802216073204264322036
                    ),
                _s111(
                    0.017923101339837434656623428594374,
                    0.0017579725835351300112357187787217
                    ),
                _s111(
                    0.35945684199505284467011522134424,
                    0.0014756136531628703160871045952088
                    ),
                _s111(
                    0.22115398478494751870985718839346,
                    0.15366312746259246580382810346215
                    ),
                _s111(
                    0.23897319487869613799252490844711,
                    0.22762521613724917092142534043686
                    ),
                _s21(0.12980665330027676925579232756653),
                _s111(
                    0.22652288598117578105986245712184,
                    0.11420303775151179873357267012934
                    ),
                _s111(
                    0.23959052403779021901917735617855,
                    0.080375059468551332648632623387454
                    ),
                _s111(
                    0.17829140299259231148118204536714,
                    0.030618363682489993315584974575301
                    ),
                _s21(0.20634966563088314748925155384717),
                _s21(0.17508458880693992678753399039646),
                _s111(
                    0.28258668912041751457866758194858,
                    0.026942102781423060815210764374344
                    ),
                _s111(
                    0.1372571102587295761492554586383,
                    0.03279721354693261380798358674717
                    ),
                _s21(0.49950457686795946289871348708067),
                _s111(
                    0.10025441523656965740357654920637,
                    0.032727112733648902370805078926134
                    ),
                _s111(
                    0.051707476559965296585764972787764,
                    0.015160101313872511931613372179533
                    ),
                _s111(
                    0.26061485513765690322924350595772,
                    0.17088465295297748667405070295455
                    ),
                _s111(
                    0.22721432495995734560664426250011,
                    0.028663286204167158861011733426409
                    ),
                _s111(
                    0.31104048723505094070947646135304,
                    0.012695671565960890416585999988922
                    ),
                _s21(0.055220419060502929271886836493043),
                _s111(
                    0.18562387629905325188370704266923,
                    0.087901535248881014455651112376267
                    ),
                _s111(
                    0.065862643760817788659368669378353,
                    0.030720346176742388312942872303364
                    ),
                _s111(
                    0.30159047009536227383244825573472,
                    0.0028145506398763484375777952526873
                    ),
                _s111(
                    0.29254265765423633400640674103213,
                    0.079839416372291539930284386903505
                    ),
                _s111(
                    0.027302484752884941695684946407879,
                    0.013321197902767094745984491145921
                    ),
                _s21(0.29915276583659394425868888416397),
                _s111(
                    0.37157192057754425099229579764134,
                    0.009602186018115497123068276784766
                    ),
                _s111(
                    0.20100730399368711685422011083313,
                    0.0025024373852326889127392685266069
                    ),
                _s111(
                    0.29318312937226642134667664669254,
                    0.20173082806599649847619289682709
                    ),
                _s111(
                    0.23333189110362259535505409087726,
                    0.052677682300648582917945024396888
                    ),
                _s111(
                    0.17929914189700434866085770363203,
                    0.055769146110768553121590077042676
                    ),
                _s21(0.36698227725645825974458889598504),
                _s111(
                    0.29343533112909299134962123676814,
                    0.048236701309217153037741450468931
                    ),
                _s111(
                    0.3298405962756192713492005844647,
                    0.2336975835531220945406533611554
                    ),
                _s111(
                    0.1381756936115276356176859521807,
                    0.092662150608776075471210191038239
                    ),
                _s111(
                    0.35851286577610111804904261998402,
                    0.17345476810307231493454892726586
                    ),
                _s111(
                    0.42033778346654878110169049260874,
                    0.017850893646582775600505442849036
                    ),
                _s111(
                    0.27905113300348068705258340005615,
                    0.11994669249584712363358737641399
                    ),
                _s111(
                    0.093462672249371088063353386492546,
                    0.013984617251760835936215769103478
                    ),
                _s111(
                    0.38740898768901990746942401894233,
                    0.1193577643150774351467885799017
                    ),
                _s111(
                    0.40718469566152014064288423452878,
                    0.073584654985361707496446040442976
                    ),
                _s111(
                    0.18764058338970974686011423900594,
                    0.012716747848540712958933586823879
                    ),
                _s21(0.26492164239038054429020886347069),
                _s111(
                    0.15222552603827632527455209640559,
                    0.0024403804076878062427103035299982
                    ),
                _s111(
                    0.17414200178656339120180287456279,
                    0.12989062597729442096705713971281
                    ),
                _s21(0.39846906841394506406509491274352),
                _s111(
                    0.32140815517552619744993542730232,
                    0.14388063038730445357265777965525
                    ),
                _s111(
                    0.08880433617342649858219367517858,
                    0.05885980928835669300960199633859
                    ),
                _s111(
                    0.13042891555171220820608856797143,
                    0.060297341384809317560471346881513
                    ),
                _s21(0.45262019075149421835403988844769),
                _s111(
                    0.13683791272221002911004626163862,
                    0.013225452599358739575574311117718
                    ),
                _s111(
                    0.35285335917967059052260305770487,
                    0.057271259501988624359656255105434
                    ),
                _s111(
                    0.24638690021714218633840233231699,
                    0.010933647586285982188026616896802
                    ),
                _s21(0.42686104953330748500577400240201),
                _s3(),
                _s111(
                    0.34352004405947132394137323853713,
                    0.095925165089547494363764143626422
                    ),
                _s111(
                    0.41620896363646422927904442399878,
                    0.039670429715756720653343978505842
                    ),
                _s21(0.094572437133494494856778332633011),
                _s21(0.47298947086889735423796972940598),
                _s111(
                    0.039011810699985755266848630926926,
                    0.0027950477505293967180250456272175
                    ),
                _s111(
                    0.10718071413460607278898318186021,
                    0.0026666116526768370646136742346938
                    ),
                _s111(
                    0.42642381970632610063729240025728,
                    0.0025395300762577354307849974981285
                    ),
                _s21(0.03408514115679924541801096003401),
                _s111(
                    0.35196013136382962553154300401362,
                    0.028731375558288720110454742961242
                    ),
                _s21(0.48699899682578148833533003239266),
                ])
            self.weights = numpy.concatenate([
                0.00019430642477546114540474384505479 * numpy.ones(6),
                0.000059366480297833226859213763050794 * numpy.ones(6),
                0.00048068870750266316128504326245131 * numpy.ones(6),
                0.00023012493350022406032921192110457 * numpy.ones(6),
                0.00129446271706814368645722202097 * numpy.ones(6),
                0.00046052089829068583612298012538864 * numpy.ones(6),
                0.00015506900800286614875782498628692 * numpy.ones(6),
                0.00052643201689967487101796591685244 * numpy.ones(6),
                0.0040776028008771986479085962982326 * numpy.ones(6),
                0.0027273371555633171716216406685823 * numpy.ones(6),
                0.0032024171503525742635729254806445 * numpy.ones(3),
                0.0038178383600418168692377518167031 * numpy.ones(6),
                0.0032435899284886569717260817522179 * numpy.ones(6),
                0.0018849657183849070579318242446728 * numpy.ones(6),
                0.0046373488411219747471358959682485 * numpy.ones(3),
                0.0040870604654872975058210819937111 * numpy.ones(3),
                0.0021807223876069809496205654394384 * numpy.ones(6),
                0.0018345812391610436077069908718838 * numpy.ones(6),
                0.00048311352303202593125017302671545 * numpy.ones(3),
                0.0016273397954327090367712798836743 * numpy.ones(6),
                0.00082853282117272432197371636814537 * numpy.ones(6),
                0.0051074900867676870821701154429309 * numpy.ones(6),
                0.0022391133588915538439521323445791 * numpy.ones(6),
                0.0016401723644210384447685927036221 * numpy.ones(6),
                0.0016926101405552390034081430264421 * numpy.ones(3),
                0.0036716713931777372159446192564637 * numpy.ones(6),
                0.0014130049415994312303806168541055 * numpy.ones(6),
                0.00076444200086207835943035115093831 * numpy.ones(6),
                0.0039540316260532015773702017853888 * numpy.ones(6),
                0.00064483354495711560938135678054869 * numpy.ones(6),
                0.006833755375698125920392188568849 * numpy.ones(3),
                0.001554756902538077082294607890054 * numpy.ones(6),
                0.0006264797924812255561822366402375 * numpy.ones(6),
                0.0059918057501738611166075314882226 * numpy.ones(6),
                0.003058417814564012191271148650048 * numpy.ones(6),
                0.0029964079041531279696952032387586 * numpy.ones(6),
                0.0068925047985869562794340925435312 * numpy.ones(3),
                0.0032333735833172280654410756149283 * numpy.ones(6),
                0.006646364364979971474877723959709 * numpy.ones(6),
                0.0031949049329088624005422876045056 * numpy.ones(6),
                0.0061619916829777406134241700387564 * numpy.ones(6),
                0.0021838177312897504417461678565133 * numpy.ones(6),
                0.0048382287535972529666710549674985 * numpy.ones(6),
                0.0011709435195081118193000511064058 * numpy.ones(6),
                0.005407496856439837179070350200563 * numpy.ones(6),
                0.0045333891237490211741070953665009 * numpy.ones(6),
                0.0015268444445366377652559535597949 * numpy.ones(6),
                0.0063146569930329001611267812400687 * numpy.ones(3),
                0.0005985557587521112603516982542641 * numpy.ones(6),
                0.0042213374314072620160024080604321 * numpy.ones(6),
                0.0065583075001978874053755235658422 * numpy.ones(3),
                0.0055145622313394954276235079780303 * numpy.ones(6),
                0.0023371778099106668343769742364979 * numpy.ones(6),
                0.0027879889526796996342465528990632 * numpy.ones(6),
                0.0050168594859682129981700415431867 * numpy.ones(3),
                0.0014305749953705688554268289811251 * numpy.ones(6),
                0.0039720422418774463540113348827835 * numpy.ones(6),
                0.0016230214293910882471462784992365 * numpy.ones(6),
                0.0058825091261753016393580159592291 * numpy.ones(3),
                0.0069414696543850275965583717472967 * numpy.ones(1),
                0.0048308548778950909282989554918174 * numpy.ones(6),
                0.0035396788161465032816995968227414 * numpy.ones(6),
                0.0031371434621642701125645307177649 * numpy.ones(3),
                0.0040930651478196709689534178759353 * numpy.ones(3),
                0.00036706598941806133093445861124868 * numpy.ones(6),
                0.00057096990940418185655978548091091 * numpy.ones(6),
                0.00087611353916882224106990377996885 * numpy.ones(6),
                0.0012513146853204593394543879823602 * numpy.ones(3),
                0.0030871109749241656810064387182736 * numpy.ones(6),
                0.0030558057702988307263859691606889 * numpy.ones(3),
                ])
        elif index == 46:
            bary = numpy.concatenate([
                _s21(0.47916270948102112963674655694901),
                _s111(
                    0.39796976464070987374146326910061,
                    0.15412314916607718763386128751409
                    ),
                _s111(
                    0.41430760798817585247292119214798,
                    0.11313670592402336911221687663357
                    ),
                _s111(
                    0.32972824987662375141682422967114,
                    0.073018327779703224554945239571323
                    ),
                _s111(
                    0.22575924233005779173603560549228,
                    0.10980337680736350125723962255144
                    ),
                _s21(0.048683079583320023852666461154885),
                _s21(0.48684943445257851288308720301258),
                _s111(
                    0.32869934756052368461464058242351,
                    0.10283587537107290723838961214153
                    ),
                _s111(
                    0.05702451832928091702553909986946,
                    0.01045501534441005688570551365678
                    ),
                _s111(
                    0.21815349062760711114057703849767,
                    0.14752226154605636112875351735872
                    ),
                _s111(
                    0.20230634588346361682952356191044,
                    0.052126096807430030011130674006018
                    ),
                _s111(
                    0.21328320980880805000595729537009,
                    0.080084195766168846134720363241607
                    ),
                _s111(
                    0.0082133588712219659750885012117669,
                    0.0021607803546682809122092542336658
                    ),
                _s111(
                    0.46151536161230036434416581864494,
                    0.011363802349727234769544531942546
                    ),
                _s111(
                    0.35462246731270576707533967634695,
                    0.12977852361765098406747121478148
                    ),
                _s111(
                    0.31092064358603159769875464110206,
                    0.048937278607245348646142002993885
                    ),
                _s111(
                    0.30495578413323064938462348044207,
                    0.15841755900746588503712405238207
                    ),
                _s111(
                    0.25079443756390419607160674960593,
                    0.17550008543840904907256645799493
                    ),
                _s21(0.19629480922944954354331985945311),
                _s21(0.16579871656709604882560695663366),
                _s111(
                    0.43386922186448054471010795122971,
                    0.053126381310496243269366308787514
                    ),
                _s111(
                    0.077037698625319650942105904330187,
                    0.050908977358329517924066476964936
                    ),
                _s111(
                    0.27683618601952247195506771206059,
                    0.011121399597137709039082301349419
                    ),
                _s111(
                    0.17363716210342146647643822079576,
                    0.011779753680717899407831406402883
                    ),
                _s111(
                    0.33543368866137698848652488409585,
                    0.011808126303693292331432189564589
                    ),
                _s111(
                    0.37506421845239666006173672316455,
                    0.051243883990893601234481724864487
                    ),
                _s111(
                    0.34671128612990046081390225706415,
                    0.18254125647613275884858047888206
                    ),
                _s111(
                    0.2559685112606641443421715054101,
                    0.053901871450292819434717333645301
                    ),
                _s111(
                    0.29200692078853169466471209021405,
                    0.02734331346155362218688761534374
                    ),
                _s111(
                    0.23641723957957555432760661675021,
                    0.029742381879280235167878145552081
                    ),
                _s111(
                    0.46615089812943551630348676535391,
                    0.0021912946810572991343929558570205
                    ),
                _s111(
                    0.22348262914558564848090409390152,
                    0.012573665836604229269208442837794
                    ),
                _s111(
                    0.39736801387800062864950507114082,
                    0.011519669129512468209598208139418
                    ),
                _s111(
                    0.28567924324498999689170600059575,
                    0.21066803766676626497998081825783
                    ),
                _s21(0.22900711778704468509078510759346),
                _s21(0.29989810187972864792896516120506),
                _s111(
                    0.40099191229640677777210640365846,
                    0.0021708502824680004413958624077701
                    ),
                _s111(
                    0.15627351040407679102362607308105,
                    0.055766422379797911137112355090838
                    ),
                _s111(
                    0.047785893684681876566400181560147,
                    0.0020979082703016117353889276685849
                    ),
                _s111(
                    0.18224915993137408000943661488197,
                    0.029078197579592478733906760094412
                    ),
                _s21(0.36608126141142539383817281213756),
                _s111(
                    0.089347568191985939673501920182905,
                    0.01158296684402517921758240440782
                    ),
                _s111(
                    0.16521988123991276093953205920456,
                    0.087462361650272741648365976148341
                    ),
                _s111(
                    0.080528004739168421681531040536701,
                    0.002078503598219240465311275116347
                    ),
                _s111(
                    0.33054150156759908344676167623722,
                    0.23647778374478368931064679630579
                    ),
                _s111(
                    0.17507401325229365257054486064663,
                    0.12280901805371280570688487984232
                    ),
                _s111(
                    0.11339497302942240420030482225686,
                    0.055005028444483801174627200478575
                    ),
                _s111(
                    0.094544228406448840975163230740898,
                    0.028793533227112058462787171879365
                    ),
                _s111(
                    0.058958942360846319577833438172288,
                    0.026293331791867717396482084604345
                    ),
                _s111(
                    0.35544611710236966465463745323878,
                    0.028307327385223555124903285896573
                    ),
                _s111(
                    0.42062784444542893597095851299564,
                    0.028126938807138848819802698859599
                    ),
                _s21(0.030769343952081854005726466057613),
                _s111(
                    0.27634247716644574798797309338647,
                    0.0021009010930692362908788360671107
                    ),
                _s111(
                    0.13583440513813729124502621311192,
                    0.030522100402259810717066753382941
                    ),
                _s111(
                    0.12143244492038733789515612080398,
                    0.088143041981530291538804786171567
                    ),
                _s21(0.082498693763558409653969192049142),
                _s111(
                    0.12925448133294827568354596229906,
                    0.012580807272203798495784812298941
                    ),
                _s111(
                    0.27293931833032633582317093627509,
                    0.083876055657832100346937488127605
                    ),
                _s21(0.460823004347588137586530822286),
                _s111(
                    0.27622443167428703143713489781663,
                    0.12479160044470473025596191775494
                    ),
                _s111(
                    0.16653106199180106998808859444417,
                    0.0022302492894394474036461393299445
                    ),
                _s111(
                    0.33732386187895564837044622118328,
                    0.0022846671110009908452396269018243
                    ),
                _s21(0.26661570869769679859163499373998),
                _s21(0.12843943494751237412151995369469),
                _s111(
                    0.38998483768037880669859166073977,
                    0.083193302580065813710441938532391
                    ),
                _s3(),
                _s111(
                    0.21905740202960690188751108812894,
                    0.002407031091832555587141285503754
                    ),
                _s111(
                    0.12034134244309032583033862652419,
                    0.0024170636184783145191665268381002
                    ),
                _s21(0.39751157540204835249439751334302),
                _s111(
                    0.032143459785426795212588199786974,
                    0.012886871818453033026111720974283
                    ),
                _s21(0.012879698435885184778232430003051),
                _s111(
                    0.02397836258819743200201203732051,
                    0.0025873987120106583524102280903204
                    ),
                _s21(0.00049292224240589825770211789631932),
                ])
            self.weights = numpy.concatenate([
                0.0015063623878995840787742848595959 * numpy.ones(3),
                0.0044537081869330011462958733948107 * numpy.ones(6),
                0.0043200095246416356069403400666402 * numpy.ones(6),
                0.0030868796559389265310976590985829 * numpy.ones(6),
                0.0034302258154803293154984678806002 * numpy.ones(6),
                0.0012474691816453434023226787405125 * numpy.ones(3),
                0.0021895163833001088976822849286296 * numpy.ones(3),
                0.0037827619823912609405177880130901 * numpy.ones(6),
                0.00067043064006260809243287853344153 * numpy.ones(6),
                0.0041741147673405347157396403142038 * numpy.ones(6),
                0.0025583666569760223892447196932536 * numpy.ones(6),
                0.0031839869521232878496301561046634 * numpy.ones(6),
                0.0001317170927953910971186738806805 * numpy.ones(6),
                0.0016456723831060041792575859018257 * numpy.ones(6),
                0.0045954830932585576647150345281534 * numpy.ones(6),
                0.002868860156815623207459492166359 * numpy.ones(6),
                0.0047267340609202294049517674802243 * numpy.ones(6),
                0.0048087317688295177107851559139344 * numpy.ones(6),
                0.004532616518341832765856590283053 * numpy.ones(3),
                0.0041121335512160985008335327151929 * numpy.ones(3),
                0.0031927670836934522254826235207111 * numpy.ones(6),
                0.001748765050289432689250125736376 * numpy.ones(6),
                0.0014306703323830101244408070353807 * numpy.ones(6),
                0.0012681752789353323424683018055985 * numpy.ones(6),
                0.0015790337439116784450549182575903 * numpy.ones(6),
                0.0032154452084320582946231226327092 * numpy.ones(6),
                0.0053977058220169638545862040442811 * numpy.ones(6),
                0.0030627052030108834905417196409397 * numpy.ones(6),
                0.0023192184096374533548718690031638 * numpy.ones(6),
                0.0022961002901772868510970817721472 * numpy.ones(6),
                0.00073381033207120100787902455256655 * numpy.ones(6),
                0.0014437602767408076990074491016952 * numpy.ones(6),
                0.001658666323964262804518639773967 * numpy.ones(6),
                0.0055707293081132551380676481707032 * numpy.ones(6),
                0.0053362839978633247591789926771485 * numpy.ones(3),
                0.0064756219278226225841251386402228 * numpy.ones(3),
                0.00072143743446459932675274939117929 * numpy.ones(6),
                0.0026221278244147779937871777310052 * numpy.ones(6),
                0.00030357960337238272680280449256811 * numpy.ones(6),
                0.0021000870011484094570320184487611 * numpy.ones(6),
                0.0065283169454314222163143077525831 * numpy.ones(3),
                0.00096629861874133241763299878042671 * numpy.ones(6),
                0.003204672370228977644264384958221 * numpy.ones(6),
                0.00039229504538734212815337409106661 * numpy.ones(6),
                0.0062575812080364894296332922288324 * numpy.ones(6),
                0.0037812405140909231481317649352061 * numpy.ones(6),
                0.0023063619647721630748744073715299 * numpy.ones(6),
                0.0015964303372989401062356703566916 * numpy.ones(6),
                0.0012488360135939412716068355298553 * numpy.ones(6),
                0.0025556284082945359054032216093201 * numpy.ones(6),
                0.0026782849316905531295214781967434 * numpy.ones(6),
                0.00097170970013712346298731919466367 * numpy.ones(3),
                0.00064139001583378701299190034920906 * numpy.ones(6),
                0.0019107173540585363561823976373209 * numpy.ones(6),
                0.0030285648566489805418792360425348 * numpy.ones(6),
                0.0024386979014754190908252840206301 * numpy.ones(3),
                0.0011940194745137679980470233094651 * numpy.ones(6),
                0.0039897861442201054909716681331521 * numpy.ones(6),
                0.0041882040717763930426753738742527 * numpy.ones(3),
                0.0047741598764475305418668164513622 * numpy.ones(6),
                0.00056810129753066170737390082505091 * numpy.ones(6),
                0.00072928978989153362696702902228585 * numpy.ones(6),
                0.0059862742964751043786625674487063 * numpy.ones(3),
                0.0037734801565330675130108717738378 * numpy.ones(3),
                0.0043995561699860475327147397464839 * numpy.ones(6),
                0.0066350471232071649309321846175091 * numpy.ones(1),
                0.00067895899857733322628322313151918 * numpy.ones(6),
                0.00053144411340617999157724169354652 * numpy.ones(6),
                0.0064816472544865440176911821271735 * numpy.ones(3),
                0.00064752466341334005231245622217073 * numpy.ones(6),
                0.00044934671028004461890532674038755 * numpy.ones(3),
                0.00026132794600230058566502965550813 * numpy.ones(6),
                0.000014095220135950496697652632292852 * numpy.ones(3),
                ])
        elif index == 47:
            bary = numpy.concatenate([
                _s21(0.3491775041543909790546580008642),
                _s111(
                    0.013218234774261483880563808729232,
                    0.006926108118994644334700155718127
                    ),
                _s111(
                    0.18207414326155922162264149992621,
                    0.08163584711098549313367835050263
                    ),
                _s111(
                    0.15058092731770140208183825581054,
                    0.086196365988733004132401047583861
                    ),
                _s21(0.30289654411791439802927939071405),
                _s111(
                    0.076432982656707521762251965084946,
                    0.050293180530191378455443720453388
                    ),
                _s111(
                    0.013832616097629311784002845676612,
                    0.0010488989499174364759595447536553
                    ),
                _s111(
                    0.026745227243076112863827048104564,
                    0.010020497565352199950108325944513
                    ),
                _s21(0.050531902462301472961280761673615),
                _s111(
                    0.37119023675251458014283465483585,
                    0.20367425304517983233209726502225
                    ),
                _s111(
                    0.36889497279211839146594458654268,
                    0.15932668387393361126785014591288
                    ),
                _s21(0.42150939287870942170315458657603),
                _s111(
                    0.10707457989312000466039003490769,
                    0.053368702564695919885493789140984
                    ),
                _s111(
                    0.11521791216736023916459518580346,
                    0.085296064492316381580455524963577
                    ),
                _s21(0.4746270243429588183553206091011),
                _s111(
                    0.41796285396989951116360329967426,
                    0.051482297929066756975958549083355
                    ),
                _s21(0.079572842069243966745980365581543),
                _s111(
                    0.046299259915931931454742587191256,
                    0.027553631439473461886250244729638
                    ),
                _s111(
                    0.32502993731599943656491563891551,
                    0.1187020271689407972326767838452
                    ),
                _s111(
                    0.22302145685345382442214304094589,
                    0.085061577169692984295460202615737
                    ),
                _s111(
                    0.38224008025178311513671300816117,
                    0.1181275084411923325003327051435
                    ),
                _s111(
                    0.31868101680405751661410125655152,
                    0.25270835065274156757849867563484
                    ),
                _s111(
                    0.31610242919729198469524287196934,
                    0.20574594120817155864182448417888
                    ),
                _s21(0.023680229143040353443897153336515),
                _s111(
                    0.20171436026250172109727951558338,
                    0.028645737105414887282258468584719
                    ),
                _s111(
                    0.30409537445729691733416589142432,
                    0.028700825091621764132858755456921
                    ),
                _s111(
                    0.25094060009445292481777841566755,
                    0.029082235848793492425415417927053
                    ),
                _s111(
                    0.27011368548774695898543722418818,
                    0.12174223594066044670469195110747
                    ),
                _s111(
                    0.42300788677184292033602581998144,
                    0.027946918436998443915521784510104
                    ),
                _s111(
                    0.21883810840826200581887640795028,
                    0.12214006675837241710664033066324
                    ),
                _s21(0.26305643577037230963841107569773),
                _s111(
                    0.31555828243015628787067249778178,
                    0.1604129850910246252818044467694
                    ),
                _s111(
                    0.35999420507678677290272092836061,
                    0.052038405103005515782563411719102
                    ),
                _s111(
                    0.07614928056270968996813890753124,
                    0.02800809732643471896104111835132
                    ),
                _s111(
                    0.36168712504285654152408799094691,
                    0.028364199293252952420287693788423
                    ),
                _s21(0.44206939726383057448794820809999),
                _s21(0.48616010748368979260854540405771),
                _s111(
                    0.39695704297947914373845888090706,
                    0.081814955021604451980803985756895
                    ),
                _s111(
                    0.26358210465038579911281772367768,
                    0.21213389631656112528730298307919
                    ),
                _s111(
                    0.19253505421227116020692504100516,
                    0.052445204222916209419919532019483
                    ),
                _s111(
                    0.15505379800193294679355349091485,
                    0.029177676432732593138719051267614
                    ),
                _s111(
                    0.080813918367061020464083412636468,
                    0.011675295489491431741242417630664
                    ),
                _s21(0.1621461225224503712090069361361),
                _s111(
                    0.049490300569582763830851873154054,
                    0.011352955184895583919539707187113
                    ),
                _s21(0.3736911818051820485235093533903),
                _s111(
                    0.30136611917506083799335566121613,
                    0.052741208654800462247560743843546
                    ),
                _s21(0.21255126890481380911844937823889),
                _s111(
                    0.14555909412430643463148475579277,
                    0.053936384497133750706414368046369
                    ),
                _s111(
                    0.033923931145958407107827110094902,
                    0.0020268965838006197624064306738152
                    ),
                _s21(0.45983217414950114278038588622358),
                _s111(
                    0.46092039829835215208970734329879,
                    0.01137411397875084568952941056824
                    ),
                _s111(
                    0.20965063699274263075580658523667,
                    0.16471973579805898230045672591204
                    ),
                _s111(
                    0.33561883072763637451853163906892,
                    0.082485282777811566536055664871459
                    ),
                _s111(
                    0.11919835609480266621264490605608,
                    0.011866025809438346562063558216507
                    ),
                _s111(
                    0.1127517427175889428481769307165,
                    0.028988580123217919429028820086793
                    ),
                _s111(
                    0.24460448499371052160980060745837,
                    0.053609251970668332280372771625139
                    ),
                _s111(
                    0.26132789521301542225176186574986,
                    0.16496729258638806192431138547285
                    ),
                _s111(
                    0.27671022235884704094186716251008,
                    0.084133458635387945623559439025061
                    ),
                _s111(
                    0.33105769830551395851847559248892,
                    0.01171814078526501984577070408775
                    ),
                _s111(
                    0.17116596035352305218181401945689,
                    0.11982515927538520791444492887601
                    ),
                _s111(
                    0.21483809354041570634759942126295,
                    0.011808474050901097097460558176917
                    ),
                _s111(
                    0.16403259578823730507383316797078,
                    0.011898235715954000890737247152243
                    ),
                _s111(
                    0.27059952051927449539518883541036,
                    0.011847600207995976283509689032171
                    ),
                _s111(
                    0.39488152990542129587234901672367,
                    0.011529446726999968332180197223238
                    ),
                _s21(0.49891852914867343689305454109148),
                _s111(
                    0.097680144060963022562821552919871,
                    0.0022497058285754676732778334640376
                    ),
                _s111(
                    0.43166318847649614546565250880266,
                    0.0021780831191386290025657683691268
                    ),
                _s111(
                    0.14032386057258031434573715887238,
                    0.0022594310974725850456428439110756
                    ),
                _s21(0.12497450674891139718264380483543),
                _s111(
                    0.24379690010750172182086633436756,
                    0.0022489510449105549561676898083058
                    ),
                _s111(
                    0.061993993437915255570298247701841,
                    0.0022107829440679360539299194684872
                    ),
                _s111(
                    0.30295049284199041481259435863913,
                    0.002249931429543905231134559812413
                    ),
                _s111(
                    0.18926610050291029549352523762224,
                    0.0022640392620407878206101963750957
                    ),
                _s111(
                    0.36591902793680910877216118616738,
                    0.0022135119210075649837476010203879
                    ),
                _s21(0.0026050966947361450393817016621683),
                ])
            self.weights = numpy.concatenate([
                0.0043146569583335191928416558288118 * numpy.ones(3),
                0.00021278200652494503520676750690905 * numpy.ones(6),
                0.0022514827449286118831613255724 * numpy.ones(6),
                0.0022650701154880334365864459296083 * numpy.ones(6),
                0.0046991712312519045108139481482327 * numpy.ones(3),
                0.0014027295746802449698148589758296 * numpy.ones(6),
                0.000097809916812659808848384197703869 * numpy.ones(6),
                0.00040386298680349287282029663034889 * numpy.ones(6),
                0.0012364924097053459366789308343258 * numpy.ones(3),
                0.0051536553206130741934589788397222 * numpy.ones(6),
                0.0045094710563839332380248757971626 * numpy.ones(6),
                0.0045159973781931017025311421733998 * numpy.ones(3),
                0.0019145889672060224873888568247877 * numpy.ones(6),
                0.0026062651829993018541367641995678 * numpy.ones(6),
                0.0029444810662742586298524963681758 * numpy.ones(3),
                0.0030761572742866454978974091759861 * numpy.ones(6),
                0.0020830024863947298660339642673581 * numpy.ones(3),
                0.0010328123125118905165924135450689 * numpy.ones(6),
                0.0043543612304128450980284778859083 * numpy.ones(6),
                0.0033920879178230409638046640804537 * numpy.ones(6),
                0.0045331567046035262637793561666597 * numpy.ones(6),
                0.0053697193075646740400657855686553 * numpy.ones(6),
                0.0050112518902676699502000555320553 * numpy.ones(6),
                0.00066300786326677415212877612184054 * numpy.ones(3),
                0.0019533942747888751882501686996306 * numpy.ones(6),
                0.0022701532915994110723259782283175 * numpy.ones(6),
                0.0021386729830484870777580573012001 * numpy.ones(6),
                0.0042781341832724128287064793325955 * numpy.ones(6),
                0.0025070193444428064880216267092171 * numpy.ones(6),
                0.0039737311003984914735456718844516 * numpy.ones(6),
                0.0055982544326030759989760431223675 * numpy.ones(3),
                0.0047404136863378226558762346593905 * numpy.ones(6),
                0.0031609021577844786273243368963539 * numpy.ones(6),
                0.0013064930317837199100776856536544 * numpy.ones(6),
                0.0024224165705940803464167788260539 * numpy.ones(6),
                0.0045518145135698885309983429960709 * numpy.ones(3),
                0.0025072813505808650659272690871947 * numpy.ones(3),
                0.0041337585903884473624377879918706 * numpy.ones(6),
                0.0051051331812112277965493883493761 * numpy.ones(6),
                0.0026976233339834156557930933192635 * numpy.ones(6),
                0.0018809643979072983211112721640693 * numpy.ones(6),
                0.00091421901027426237659098695196042 * numpy.ones(6),
                0.0040644299258672032178497638828798 * numpy.ones(3),
                0.00070106560286908116563842952918107 * numpy.ones(6),
                0.0055051135413514141373586069309779 * numpy.ones(3),
                0.003190105648715051436395570544654 * numpy.ones(6),
                0.0049446615998697491865097666478942 * numpy.ones(3),
                0.0024338258411055340119124251598528 * numpy.ones(6),
                0.00025084745145746957978885278232108 * numpy.ones(6),
                0.0040657898075941602101518612653261 * numpy.ones(3),
                0.0017108131277605055948580984738201 * numpy.ones(6),
                0.0045280121334294017779981887937264 * numpy.ones(6),
                0.0040081392301490437888476750433182 * numpy.ones(6),
                0.0011195403999073294455905602741588 * numpy.ones(6),
                0.0016556490812276695490161147635792 * numpy.ones(6),
                0.0030745572112735754912785136882923 * numpy.ones(6),
                0.0048300704241971473390133650616931 * numpy.ones(6),
                0.0038986974038214249817258777400309 * numpy.ones(6),
                0.0016458519911547444785505661006025 * numpy.ones(6),
                0.0037119048709274736384122003829731 * numpy.ones(6),
                0.0014185756078616090601015039156472 * numpy.ones(6),
                0.0012890581955709145109213289099414 * numpy.ones(6),
                0.0015614900549057867686283948220674 * numpy.ones(6),
                0.0016952718909460044498760030036021 * numpy.ones(6),
                0.00074883086120397942849860954395469 * numpy.ones(3),
                0.00045280227256791726315439594618803 * numpy.ones(6),
                0.00074533983106636671568553223119497 * numpy.ones(6),
                0.00053234814725867168787023717252183 * numpy.ones(6),
                0.003632139855663754981931331678907 * numpy.ones(3),
                0.00065761968362639610357428781649619 * numpy.ones(6),
                0.00036206905113067701539554677153421 * numpy.ones(6),
                0.00070601194567049312732844208051029 * numpy.ones(6),
                0.00060187710221997117773774997688409 * numpy.ones(6),
                0.00073263934654966526183069298492737 * numpy.ones(6),
                0.000089253661418059160787189703718929 * numpy.ones(3),
                ])
        elif index == 48:
            bary = numpy.concatenate([
                _s21(0.34772445871080945921858087263637),
                _s111(
                    0.36956729159016621347427614942414,
                    0.21570064416961684136937961138523
                    ),
                _s111(
                    0.3084818979933622488395569991715,
                    0.10830149992200391895429803833053
                    ),
                _s111(
                    0.3196427292245513327190901671665,
                    0.034695952771136395220571115533313
                    ),
                _s21(0.37343127125532209807137430019027),
                _s111(
                    0.0013460639309496434449913949535816,
                    0.0023979705918984515228087151480167
                    ),
                _s111(
                    0.32808481780883128583518538926516,
                    0.019170784080321566960562147174653
                    ),
                _s111(
                    0.31981290827271948514679785328427,
                    0.27636156531080965313873694597749
                    ),
                _s111(
                    0.15892660136182003023150105360172,
                    0.0016411683438174247902165307454497
                    ),
                _s111(
                    0.19627597719723563119763203960289,
                    0.053323193421323463215766804575727
                    ),
                _s21(0.48933905356285323189755068911384),
                _s111(
                    0.29334218459018172786281846228406,
                    0.082928836401138320605645137556541
                    ),
                _s21(0.048410924825140003549531873101234),
                _s21(0.080159548747975685998971495847438),
                _s111(
                    0.12472732588139755964792292690437,
                    0.023319017149181456732602590123587
                    ),
                _s111(
                    0.31159324594327961154094500023497,
                    0.13871309395938550187253910064601
                    ),
                _s21(0.4790902659592135284977427301879),
                _s21(0.46664241534843466585864544118995),
                _s111(
                    0.1105414410318807724923808541231,
                    0.010540639408122495091558873512415
                    ),
                _s111(
                    0.29741719627591727905453883419723,
                    0.010347250567507892711563605899565
                    ),
                _s111(
                    0.29406805865227855377342870582411,
                    0.053781881658776657835797077086018
                    ),
                _s111(
                    0.11140121760631086717111461667156,
                    0.069760358020446372996210595748298
                    ),
                _s111(
                    0.44374345386362933341008860035434,
                    0.015502302779250982036104173300093
                    ),
                _s111(
                    0.3676357119294025251719024329336,
                    0.039622507258816878821754942534388
                    ),
                _s111(
                    0.04838341889362429270014331747278,
                    0.011810700923664012376069823553521
                    ),
                _s111(
                    0.42619854603484827823495662248263,
                    0.035287892346372193986783028257057
                    ),
                _s111(
                    0.34640370750834257306753607827504,
                    0.064715093001975105538156759324634
                    ),
                _s111(
                    0.076168644079201458271328490462095,
                    0.010422571567642056697267948968709
                    ),
                _s111(
                    0.31985580598920702216025018816017,
                    0.2271049175306688832046709924165
                    ),
                _s111(
                    0.11497674997632783103739697051906,
                    0.043154888080570560164983517001927
                    ),
                _s111(
                    0.15534118748633498969652388576385,
                    0.048974948218390376682219659456535
                    ),
                _s111(
                    0.15366413947806527019328559579452,
                    0.0094545488420266123295429943693395
                    ),
                _s111(
                    0.079853820771823415057216664214541,
                    0.001894305879112616250909584511808
                    ),
                _s21(0.011704574467113419693380165432123),
                _s111(
                    0.19786011212854839726667211452286,
                    0.01202832240369283999959763054187
                    ),
                _s111(
                    0.31691538038078183286172803130958,
                    0.0019044677779362753109165814597736
                    ),
                _s111(
                    0.16456397716341348003463484724313,
                    0.026205223792356400527934433412733
                    ),
                _s111(
                    0.36112983597993795082271628018747,
                    0.094470921524026361692054346092907
                    ),
                _s111(
                    0.41982597745422942253932528739348,
                    0.09632334903378368129660417430849
                    ),
                _s111(
                    0.24317522457466311883338098170251,
                    0.053490317263427349155838774774726
                    ),
                _s111(
                    0.14598322311364346252423416424198,
                    0.079524478600837120920709233492172
                    ),
                _s111(
                    0.02723928326798026400026025802696,
                    0.011030754781471895242965502314791
                    ),
                _s111(
                    0.24531370896477844967907467936246,
                    0.01218623296397879904889604395132
                    ),
                _s21(0.2678465137111688515074985190996),
                _s111(
                    0.11655421907345294491298292941702,
                    0.002113218262477547388525105532789
                    ),
                _s111(
                    0.053585682870044646865439262181293,
                    0.027218990658071144278383544971874
                    ),
                _s111(
                    0.085135518150016901363844301959048,
                    0.026256259012472017664350753627913
                    ),
                _s111(
                    0.4084926257856330763559307971186,
                    0.062475097450370377288537586231733
                    ),
                _s111(
                    0.38225104509481382758638848053974,
                    0.17377620987619342492790827622788
                    ),
                _s111(
                    0.025974032206589028985905447816123,
                    0.0020491260430695296428935965616937
                    ),
                _s111(
                    0.38568918876006436248719897524023,
                    0.018890041661542007346634550479093
                    ),
                _s111(
                    0.19074647943140314199705742149581,
                    0.082618952945698380026475861936674
                    ),
                _s111(
                    0.049314651604851817057279859185433,
                    0.0022799567867530885140044089538634
                    ),
                _s111(
                    0.42648254549804660960050559038719,
                    0.0047722521203025807622353368283012
                    ),
                _s111(
                    0.15136475267056509687327616076048,
                    0.11620235641158798962523053133241
                    ),
                _s111(
                    0.24108556039098643770541227087692,
                    0.085525264242367015354689057254146
                    ),
                _s111(
                    0.21203621120547973015996857160719,
                    0.029776788322806571532698105681461
                    ),
                _s111(
                    0.077342028084253878962970082698106,
                    0.050253537035871350214574119059837
                    ),
                _s111(
                    0.26592702365453135474205632070986,
                    0.21630020710235432060500142819858
                    ),
                _s21(0.21313827248019083386472426606489),
                _s111(
                    0.01027811240403521337792704128211,
                    0.0023182076844770317105552550801207
                    ),
                _s21(0.43362666136661098376959718637171),
                _s111(
                    0.36374868048981512450620730464798,
                    0.0064301877897616796009993759675447
                    ),
                _s111(
                    0.36982020510671371195440026435992,
                    0.13311111506547137968496744179646
                    ),
                _s21(0.15990091444053889623264586185501),
                _s21(0.028063663492678269556170797287994),
                _s111(
                    0.26223685562554154020644420340496,
                    0.16825248537156288995388320897387
                    ),
                _s111(
                    0.20586527419899884232721988359442,
                    0.0023776437008401596162292156160215
                    ),
                _s111(
                    0.32106129923917116232439622478747,
                    0.17927517691473483569458311796422
                    ),
                _s111(
                    0.26700892535281116636621424148074,
                    0.028463351612019049647948281285426
                    ),
                _s111(
                    0.25860454101035106910572440188761,
                    0.0022264670918958109330618409104635
                    ),
                _s111(
                    0.19962476267736045586731428102304,
                    0.12036436957224947280704771240853
                    ),
                _s111(
                    0.20889086438440644887001385363093,
                    0.16433446168671559607793204014215
                    ),
                _s21(0.10804712924230600644035619507315),
                _s111(
                    0.25333360060180893373128751394599,
                    0.12388745883767666667415588852488
                    ),
                _s21(0.49717563723846608315257060863406),
                _s111(
                    0.46379612274656633860354149770257,
                    0.00016001038938425443311552267327327
                    ),
                _s111(
                    0.38371779905054644691598940501265,
                    0.00010850920174364424913717116352104
                    ),
                ])
            self.weights = numpy.concatenate([
                0.004365138753782736588387763532721 * numpy.ones(3),
                0.0037199314122967066649237027104757 * numpy.ones(6),
                0.0030326109445702034449349190324064 * numpy.ones(6),
                0.0018470108648334892583979500487599 * numpy.ones(6),
                0.0039691750453025813664361239115348 * numpy.ones(3),
                0.000023772908084865779657340845432758 * numpy.ones(6),
                0.0015449674526409663225223655293904 * numpy.ones(6),
                0.0050318869826993904790280485786729 * numpy.ones(6),
                0.00038873576680098335179002800670132 * numpy.ones(6),
                0.0022929693353194649487983103987857 * numpy.ones(6),
                0.0018221145367845095413207306278064 * numpy.ones(3),
                0.0032590081957709469128477082913635 * numpy.ones(6),
                0.0012114472500649791090372338192252 * numpy.ones(3),
                0.0019309544896270690998024562387319 * numpy.ones(3),
                0.001299434393401518008061313240849 * numpy.ones(6),
                0.0040442317991415074851237003429078 * numpy.ones(6),
                0.002640024375641072353889254549975 * numpy.ones(3),
                0.0032841920864060679401283582412589 * numpy.ones(3),
                0.00086447910709451905824359191798643 * numpy.ones(6),
                0.0012659714364757416887984495397779 * numpy.ones(6),
                0.0027425885952415113590652000522662 * numpy.ones(6),
                0.0021154075536657444564587203452719 * numpy.ones(6),
                0.0016231545272856842416399344573229 * numpy.ones(6),
                0.0025412737334052936150968257690468 * numpy.ones(6),
                0.00062263547211217268605960399934217 * numpy.ones(6),
                0.0026418403799762712261685623608238 * numpy.ones(6),
                0.0032653171027573867557654463773852 * numpy.ones(6),
                0.00076551458592241779961600497824702 * numpy.ones(6),
                0.0051641360944715229935894634481422 * numpy.ones(6),
                0.0018554389812758353651935553820972 * numpy.ones(6),
                0.0022798100113623112887710333853045 * numpy.ones(6),
                0.0010267946171960897441063499706024 * numpy.ones(6),
                0.00033284312112662229422513601197115 * numpy.ones(6),
                0.00032294495843481555534957163991891 * numpy.ones(3),
                0.0012157493114517174528204742123372 * numpy.ones(6),
                0.00058931157289567738799989525642262 * numpy.ones(6),
                0.0017342658044948606086221400688213 * numpy.ones(6),
                0.0040498750757604590132385197811776 * numpy.ones(6),
                0.0042742322977128438717633479963017 * numpy.ones(6),
                0.0028335085475802439343267491864638 * numpy.ones(6),
                0.002743812545883741128792092508501 * numpy.ones(6),
                0.00047302145635505195586232787053933 * numpy.ones(6),
                0.0013733177441610052918846835798366 * numpy.ones(6),
                0.0054310463039951051770483299662245 * numpy.ones(3),
                0.0004233009914606136797065923668114 * numpy.ones(6),
                0.0010664611871638438230133205587207 * numpy.ones(6),
                0.001382401073311702664124279345323 * numpy.ones(6),
                0.0036585117930269911187890085944361 * numpy.ones(6),
                0.0052100048548665774797246777954799 * numpy.ones(6),
                0.00020731056059076486342891119378285 * numpy.ones(6),
                0.0019583162868631561789802565135373 * numpy.ones(6),
                0.0032089225299187570565373788686298 * numpy.ones(6),
                0.00031441515461318483006083292790225 * numpy.ones(6),
                0.0010785145551204093923463476712953 * numpy.ones(6),
                0.0036004269766152113324326297432672 * numpy.ones(6),
                0.003661191417842613589867188497473 * numpy.ones(6),
                0.0021617123138750800896979298565002 * numpy.ones(6),
                0.0018010562763861775633387583260045 * numpy.ones(6),
                0.0054475592942758763722635639483587 * numpy.ones(6),
                0.0051525605444838104514570362127744 * numpy.ones(3),
                0.00013978520152437230779237575086904 * numpy.ones(6),
                0.0049930051499418295575304474207131 * numpy.ones(3),
                0.0011962011851953004816245810067358 * numpy.ones(6),
                0.0049525658235835790024417166362793 * numpy.ones(6),
                0.004243498491854159858762059419321 * numpy.ones(3),
                0.00088760056102461919691531901760991 * numpy.ones(3),
                0.0050552058738964203121492416897065 * numpy.ones(6),
                0.00060559106189920555570614397513677 * numpy.ones(6),
                0.0053577902109237166326079371660466 * numpy.ones(6),
                0.0023655352245002668955573742612554 * numpy.ones(6),
                0.00063617605564497544016787761705707 * numpy.ones(6),
                0.0041167794317583617645309382462914 * numpy.ones(6),
                0.0047435661701724518494737922666682 * numpy.ones(6),
                0.0031739283489719255519214112325817 * numpy.ones(3),
                0.0046022674251481818283838939984268 * numpy.ones(6),
                0.0012950799113895639663827350585878 * numpy.ones(3),
                0.00023500892248057527905261792531808 * numpy.ones(6),
                0.00023987367893110875148858899312939 * numpy.ones(6),
                ])
        elif index == 49:
            bary = numpy.concatenate([
                _s111(
                    0.17470985430717311394573984855263,
                    0.040008155056936224332958516723801
                    ),
                _s111(
                    0.38696505716382622439800609577711,
                    0.0086624172314517809018104428858362
                    ),
                _s111(
                    0.12606093343540540718508978043096,
                    0.0013011091483538096627332303317969
                    ),
                _s21(0.17837301622963160245881203427969),
                _s111(
                    0.2839249520005005219395417313958,
                    0.2456285053783457921738284984191
                    ),
                _s111(
                    0.33024826933446234724512561532079,
                    0.25074163617360918924271554410812
                    ),
                _s21(0.072840653606011441774813462188164),
                _s111(
                    0.25604594106815780073777138683611,
                    0.0015097374938803929479736041384535
                    ),
                _s111(
                    0.26039531613447623317352912761836,
                    0.10898338021292514796399831696298
                    ),
                _s21(0.42898835563940702767222570315939),
                _s111(
                    0.17071269637934203790887148079842,
                    0.1019687169087760640187168122413
                    ),
                _s111(
                    0.37464405809512494610593895814093,
                    0.019772288574859819555146029474527
                    ),
                _s111(
                    0.4351684572939831291631867061549,
                    0.0013541852880936327570320832408324
                    ),
                _s111(
                    0.13909239122910394136601523091313,
                    0.087642852515997409655702444217936
                    ),
                _s111(
                    0.28015065430122799427420583288414,
                    0.078415578965199478061381269385937
                    ),
                _s111(
                    0.34633118742755955749897617750584,
                    0.062747772638091691181750340331695
                    ),
                _s21(0.38116883784557031451454556853725),
                _s111(
                    0.11424479335878181837937710081964,
                    0.0071057304538390595162377544408432
                    ),
                _s111(
                    0.26014778271440394776023321331883,
                    0.040890474228270071123057244761755
                    ),
                _s111(
                    0.33859545036972212866007943487871,
                    0.032494883027560273446852131597024
                    ),
                _s111(
                    0.44052494716069208892179906086833,
                    0.0087290442899361611375901598769119
                    ),
                _s111(
                    0.085723483740763625517886035193873,
                    0.0015747586640191503001195929682239
                    ),
                _s21(0.2007383675971185893151280851638),
                _s111(
                    0.28598408663977853085732770365976,
                    0.20153220553570834090461059975155
                    ),
                _s111(
                    0.30224306485270802342030074004581,
                    0.051418163514358964636240130086242
                    ),
                _s111(
                    0.1145574257288858714909092310646,
                    0.019028748823558390651663845244439
                    ),
                _s111(
                    0.21229643001953076269850900164246,
                    0.040788059207840106101969862343129
                    ),
                _s21(0.046278223175869237589549923778028),
                _s21(0.29714308886324797750314629443625),
                _s111(
                    0.14652842970990201245169282499374,
                    0.028514611896576519802821941965943
                    ),
                _s111(
                    0.37408164610882701375977555193828,
                    0.15370865861164023288218928312174
                    ),
                _s111(
                    0.081706830399755079405994477722645,
                    0.024605320384088840162325053186091
                    ),
                _s21(0.010842121769007150534725046286933),
                _s111(
                    0.2941474788808447090823091621213,
                    0.023535501929982482343332330885678
                    ),
                _s111(
                    0.17507036948265155472722787799839,
                    0.13771672126012072029272142434755
                    ),
                _s21(0.10285338166312913209054555371543),
                _s111(
                    0.047620537805773032444267523725112,
                    0.010975310089231762213865934431625
                    ),
                _s111(
                    0.31757248151332692463716126070082,
                    0.1604656572027247772116914906645
                    ),
                _s111(
                    0.28813805399507747519107650478994,
                    0.13294926773769973138990058674476
                    ),
                _s111(
                    0.26797088930329615007796561113548,
                    0.008754951231311353670281398617561
                    ),
                _s21(0.48417356341669489101694326843861),
                _s111(
                    0.23937905920302690614113799223895,
                    0.020743761983356521452348535612853
                    ),
                _s111(
                    0.19185327156970310253803935259433,
                    0.01976028968407493892033733395874
                    ),
                _s111(
                    0.38732424253418504153765474766158,
                    0.042939285776866986285250838623074
                    ),
                _s111(
                    0.35353241582841083999089927149662,
                    0.11789241815829039711701717133484
                    ),
                _s111(
                    0.052177489014856154906828555592733,
                    0.025374729462583658155644865632684
                    ),
                _s111(
                    0.32370211615131400594758217787274,
                    0.0935293820848459614364090960853
                    ),
                _s111(
                    0.37205936067484861776013353947085,
                    0.0019871034880730657811630942342162
                    ),
                _s111(
                    0.41369253460710497234793028784417,
                    0.10869019629350585534293618319367
                    ),
                _s111(
                    0.25086544364143718334668895740184,
                    0.17490233415096355473534193382083
                    ),
                _s111(
                    0.23684770353878356875409066779974,
                    0.06696591337164844954030779899705
                    ),
                _s111(
                    0.10995768304264821108381509267487,
                    0.040979133770732177318460744032073
                    ),
                _s111(
                    0.42973190289518155391801908130346,
                    0.024466415075845744605454960479233
                    ),
                _s111(
                    0.14376803743337255017665637597753,
                    0.056688088341008345204237954346032
                    ),
                _s111(
                    0.077410677998690816733404047343677,
                    0.0096524833615463131610528440437161
                    ),
                _s111(
                    0.15570600566534693897233724745043,
                    0.011099567627321117997169955187323
                    ),
                _s111(
                    0.3128533024496512905773479327589,
                    0.0019740717879364078996339861311026
                    ),
                _s111(
                    0.026090266138179864269261499439275,
                    0.01105886728057323507399589015529
                    ),
                _s111(
                    0.16567286081799158433167107427954,
                    0.0024625444223661337714248197094087
                    ),
                _s21(0.49326119556605273496159794202774),
                _s111(
                    0.074225754676311409846855961346731,
                    0.046395012982637350939376499713998
                    ),
                _s111(
                    0.028833917237014066202882056790053,
                    0.0020785999248741445544508967015953
                    ),
                _s111(
                    0.18682266067231509396444985771883,
                    0.068353296798321981662724994784735
                    ),
                _s111(
                    0.32708392818568737601342082534572,
                    0.010676304302980427864595576368219
                    ),
                _s111(
                    0.10389110343382679653355131076642,
                    0.069169283463613255295201153175359
                    ),
                _s111(
                    0.21553428824295463757262368797902,
                    0.10066880943464366221450450149799
                    ),
                _s21(0.40573247197552248683179918336621),
                _s111(
                    0.34315161302738998488799982276295,
                    0.20121793752844419256083933676968
                    ),
                _s111(
                    0.22240875872332544415389091604099,
                    0.14345324849263682159854043453268
                    ),
                _s111(
                    0.053480976914394243630846605195786,
                    0.0021251605376942635500295773332177
                    ),
                _s21(0.13042164492513155739074958008957),
                _s111(
                    0.43878896884502699425680059090761,
                    0.051422502449401989467464647761887
                    ),
                _s21(0.4614348056302962608030886921859),
                _s21(0.027573960243111313171581571411679),
                _s21(0.23184244537547400520409359685337),
                _s21(0.35172177662777146824478383066669),
                _s21(0.0022124213099818633010519716967457),
                _s111(
                    0.39193210137555588463970906804864,
                    0.078152632925578181728753668795839
                    ),
                _s111(
                    0.21108217989703460094137890828527,
                    0.0064969244459217452254521572818476
                    ),
                _s111(
                    0.011709857310951749891488305667151,
                    0.0020903064191528306423590100667678
                    ),
                _s21(0.49878956469151284312709763011092),
                _s111(
                    0.20638335375998988473074346596068,
                    0.00013986351814746797526502225778382
                    ),
                ])
            self.weights = numpy.concatenate([
                0.0015179493827508213114659732947223 * numpy.ones(6),
                0.0009640363088870214883579998184213 * numpy.ones(6),
                0.00028135036653219960440760330724168 * numpy.ones(6),
                0.0032298776744524847178047902458543 * numpy.ones(3),
                0.0046862213374651167825630656646937 * numpy.ones(6),
                0.0045178989983832158869339744370649 * numpy.ones(6),
                0.0015790917860183804506508431780044 * numpy.ones(3),
                0.00042641348129921866997375891535468 * numpy.ones(6),
                0.0032468366301052171593609009764572 * numpy.ones(6),
                0.0039309054788202352881917108092192 * numpy.ones(3),
                0.0027113498937414436834843241034594 * numpy.ones(6),
                0.0016327582269424586233004363566329 * numpy.ones(6),
                0.00047903017617903526332649812576226 * numpy.ones(6),
                0.0026002718294388297276002981325151 * numpy.ones(6),
                0.0030738076620057483277061046405966 * numpy.ones(6),
                0.0029878198502272939791152427799343 * numpy.ones(6),
                0.0049958592249661608448166769280989 * numpy.ones(3),
                0.00069259134967855788057321242474465 * numpy.ones(6),
                0.0022641894958965860686532471750589 * numpy.ones(6),
                0.0021683619372713905326688489832259 * numpy.ones(6),
                0.0012603344824756224563404737737502 * numpy.ones(6),
                0.00030192800544203797296835148081379 * numpy.ones(6),
                0.0040779189028797265914562248451031 * numpy.ones(3),
                0.0045616230274563116513731471352509 * numpy.ones(6),
                0.0026431427639030697708258800897681 * numpy.ones(6),
                0.0011314599292834434349353229611448 * numpy.ones(6),
                0.0021088731883818778511348329688981 * numpy.ones(6),
                0.0011717865482744910367339140286932 * numpy.ones(3),
                0.0050840537579920432445421110557325 * numpy.ones(3),
                0.0015608187011213255833546828062745 * numpy.ones(6),
                0.0046308511124103565544815766148062 * numpy.ones(6),
                0.0011648378245233953418983017287693 * numpy.ones(6),
                0.00029099343979125299793735646969871 * numpy.ones(3),
                0.0019078475871128063055909730738128 * numpy.ones(6),
                0.003620735933035441771986005237527 * numpy.ones(6),
                0.0024679826387392238910812958436766 * numpy.ones(3),
                0.00060470591661133424986723710372684 * numpy.ones(6),
                0.0044715780288225674533568785420637 * numpy.ones(6),
                0.0043217865559445603325388004005922 * numpy.ones(6),
                0.0011834459680028650438959826431482 * numpy.ones(6),
                0.0024475548918751539996868099405226 * numpy.ones(3),
                0.0017098461192591258043031422442127 * numpy.ones(6),
                0.0015260667775319674130816907713276 * numpy.ones(6),
                0.0027529020694635879528657889711598 * numpy.ones(6),
                0.0043632272367912435199986099483834 * numpy.ones(6),
                0.0009662853562325522268314130814568 * numpy.ones(6),
                0.0037815526680686423619164800635178 * numpy.ones(6),
                0.00060356931274693017588456011655511 * numpy.ones(6),
                0.004255163960838959651364142063852 * numpy.ones(6),
                0.004849053322409156023403933608078 * numpy.ones(6),
                0.0030254918153351583666722842815792 * numpy.ones(6),
                0.0018302659109491988291671863954219 * numpy.ones(6),
                0.0023112911452724827762723890957858 * numpy.ones(6),
                0.0022976181282931051426351159671405 * numpy.ones(6),
                0.00077374260145531970790743802292086 * numpy.ones(6),
                0.0011091107906632695132813104499629 * numpy.ones(6),
                0.00059172800707874940161563541003125 * numpy.ones(6),
                0.00046632907833982510514261674663255 * numpy.ones(6),
                0.00050605451421417245451892205321647 * numpy.ones(6),
                0.0016834711767858439565587456713885 * numpy.ones(3),
                0.0015462052435902458077809739396969 * numpy.ones(6),
                0.00022381677252706233778626171419472 * numpy.ones(6),
                0.0029561851838939697150170697887149 * numpy.ones(6),
                0.0014290418491702153822077169253153 * numpy.ones(6),
                0.0021458693484175845449866560026665 * numpy.ones(6),
                0.0035847674194569904228072183121255 * numpy.ones(6),
                0.005365422540602800704617431316119 * numpy.ones(3),
                0.0054527357782272263168610682616803 * numpy.ones(6),
                0.0043816657756571581941513425296911 * numpy.ones(6),
                0.00030717453596445567296367281998968 * numpy.ones(6),
                0.0034992680317616632376625209719692 * numpy.ones(3),
                0.0032381770164589147964246126866413 * numpy.ones(6),
                0.0039600197915087516011590363123491 * numpy.ones(3),
                0.00083195088890381481473975469357511 * numpy.ones(3),
                0.0053644090708949186573259370221044 * numpy.ones(3),
                0.006126773607554056538430828980122 * numpy.ones(3),
                0.000064611429908450439856021862114718 * numpy.ones(3),
                0.0041053828589275478155890999172918 * numpy.ones(6),
                0.0010683587323421039182719089759252 * numpy.ones(6),
                0.0001423198916149831062073835626988 * numpy.ones(6),
                0.00078855516545326865646794205230635 * numpy.ones(3),
                0.00016055847055223060984911012924419 * numpy.ones(6),
                ])
        else:
            assert index == 50
            bary = numpy.concatenate([
                _s21(0.25635586338479089235345337173654),
                _s111(
                    0.4433230868583420661077936032547,
                    0.00092245610607042072056740744921666
                    ),
                _s111(
                    0.025936893875790419143576184473258,
                    0.022437991262313005453542213083744
                    ),
                _s111(
                    0.24014410526326950653581197245268,
                    0.22279072829373281690334440934239
                    ),
                _s111(
                    0.13455261802664274627360994606217,
                    0.023120213970246226060924643867322
                    ),
                _s111(
                    0.045035572467331431764630960300307,
                    0.00085677979500894595451692790971656
                    ),
                _s111(
                    0.40871456656777191840753062400224,
                    0.0048086376608719361888044005854543
                    ),
                _s111(
                    0.17146345773998324858376269529448,
                    0.023493012359438530659159210605824
                    ),
                _s21(0.0016305781421461421964611794535412),
                _s111(
                    0.21010592470503317311765956803812,
                    0.093396647129363525216061396683746
                    ),
                _s111(
                    0.13929682011858504377154326138332,
                    0.041796636651148120260187447093626
                    ),
                _s111(
                    0.070382053912759211404055225007813,
                    0.049357311892861232034445516215324
                    ),
                _s21(0.49032686924845125981757877162256),
                _s111(
                    0.37260648815649406877202409079376,
                    0.15281410701216105481351354651187
                    ),
                _s111(
                    0.14791571727764409349538319733569,
                    0.065115429930621510027435127559073
                    ),
                _s111(
                    0.21402464358164372634133779332607,
                    0.023462143176020752343896396228687
                    ),
                _s111(
                    0.33026740305814401596755673571545,
                    0.052032278234411030684581762938457
                    ),
                _s111(
                    0.36998486354440394872232900383194,
                    0.0014945191055322731198374229304144
                    ),
                _s111(
                    0.067511569043336129124865214446112,
                    0.012043746934792498279581492829306
                    ),
                _s111(
                    0.43709884118009677128127146550562,
                    0.056684785543461209864012996533698
                    ),
                _s21(0.42258716850752065574380179381836),
                _s111(
                    0.41451222115778534700073832862263,
                    0.11486014204481563852540120841065
                    ),
                _s111(
                    0.36238051226892748752278385619202,
                    0.03126483618670886592836480199125
                    ),
                _s111(
                    0.31225128551120350696259703286246,
                    0.10224664014554168426180958006543
                    ),
                _s21(0.48057567643403007031974967915307),
                _s111(
                    0.39268932170368790006291715324816,
                    0.0483493018730728789062023861195
                    ),
                _s111(
                    0.026063854162288056346783051110947,
                    0.0099795769296426198913762769870566
                    ),
                _s111(
                    0.39839099584495546366329878616841,
                    0.083912940769928322505597427470674
                    ),
                _s111(
                    0.045880391320914956012098757134895,
                    0.0067281317012842103784308347222412
                    ),
                _s111(
                    0.25109830786516684924025834614286,
                    0.17730802483236840760881691737914
                    ),
                _s111(
                    0.183136933934611227877808741766,
                    0.042684087735085587651488525166963
                    ),
                _s111(
                    0.35404630327836232396485273668492,
                    0.074392768040038105189286099484577
                    ),
                _s111(
                    0.29921098711016716046469861225902,
                    0.24356894779703624729137644759716
                    ),
                _s111(
                    0.2908392329854116760501913204993,
                    0.076391569745987766072814759903787
                    ),
                _s111(
                    0.1861461318533069878290004207922,
                    0.12422458966888609555403404458714
                    ),
                _s111(
                    0.13362717497122467011448385130779,
                    0.0097790083705654904827843134073224
                    ),
                _s111(
                    0.06897633514821163925285036071902,
                    0.027943202048164079602502648106399
                    ),
                _s111(
                    0.26149101926980808607752870902388,
                    0.024808268510852052812559478835895
                    ),
                _s111(
                    0.043979456059116789067485804775854,
                    0.021126000427883826331613590940522
                    ),
                _s111(
                    0.24695723170731077076966773274437,
                    0.068149870254045245645421810008183
                    ),
                _s111(
                    0.19160730297984118595948069852465,
                    0.068795306900434617649547086511681
                    ),
                _s111(
                    0.31739598230853128042332616534552,
                    0.15635347712002139067957877191837
                    ),
                _s111(
                    0.39383930629984276059039306480199,
                    0.01578877113567470023880440670988
                    ),
                _s111(
                    0.31814466821278573006120272006377,
                    0.023446407967230234678882038916882
                    ),
                _s111(
                    0.28036448755799682750738931599361,
                    0.13592981347138613297816191027361
                    ),
                _s111(
                    0.22997573008415667745763318386899,
                    0.044695597013042020722419665438548
                    ),
                _s111(
                    0.34948812503713508382673551996045,
                    0.24605625068915952949434249373866
                    ),
                _s111(
                    0.096408357906222464877197346187873,
                    0.010763472737716835342197623887545
                    ),
                _s111(
                    0.10033800796659960653906435074785,
                    0.026797769725969450261905800938793
                    ),
                _s111(
                    0.28943784308397893210759283405063,
                    0.19990165663953924196040799420904
                    ),
                _s21(0.4595867718244118790760372825268),
                _s111(
                    0.17861645395509190959079769361239,
                    0.0096597928476049357604609869895832
                    ),
                _s21(0.11333323413507438418364497906229),
                _s111(
                    0.28615392263125203165445781094914,
                    0.043616289025622319606539499667131
                    ),
                _s111(
                    0.2293288959608979256854774739612,
                    0.0096337955497693002589136330274757
                    ),
                _s111(
                    0.35421289616647160920903580523351,
                    0.11645700215175254837476126565884
                    ),
                _s111(
                    0.34494449154498538698237514319401,
                    0.0094768771491044251917721765961197
                    ),
                _s111(
                    0.431673230823553892961330163352,
                    0.02796661629599724814584848814553
                    ),
                _s111(
                    0.22079722375183962932663549554841,
                    0.14726584419679686065993682661439
                    ),
                _s21(0.010939013064663018058015508903292),
                _s111(
                    0.15317283509768431426472684421849,
                    0.0018705669110342341426443222489978
                    ),
                _s21(0.078646361951399820046024272098613),
                _s111(
                    0.0091986167658109381134941756174822,
                    0.0020671401785266431228719998747056
                    ),
                _s111(
                    0.023545234257410343061622497195595,
                    0.0020063808420237499453300191212904
                    ),
                _s111(
                    0.34572847971846620948741398821291,
                    0.1965919922795322893721792930149
                    ),
                _s111(
                    0.25016930881215731352596628009597,
                    0.10747886610223062197660508749046
                    ),
                _s111(
                    0.25554986260163291921787392056505,
                    0.0018618093532656415646868207114604
                    ),
                _s111(
                    0.10329343526262501435355623158876,
                    0.05020729259155965149799405162955
                    ),
                _s21(0.1702588823641314358072673446768),
                _s111(
                    0.11348277919054020264668708516045,
                    0.079098766756874565359854851785091
                    ),
                _s111(
                    0.202001519139572872108006224181,
                    0.0018288696800150928796266146998579
                    ),
                _s21(0.1996334593368172280383502393609),
                _s111(
                    0.15540729334166614108614591495978,
                    0.097498460674939456244633031191421
                    ),
                _s111(
                    0.28373312584990880800207002163612,
                    0.010271042737936499112420844970971
                    ),
                _s21(0.40066893254883004766449626095416),
                _s21(0.14061387405638423652463649160293),
                _s111(
                    0.11013761071196110862096690113627,
                    0.0019450667745236041442457725518546
                    ),
                _s21(0.35119036444541556070870347215726),
                _s21(0.042250362515138067157622417020802),
                _s111(
                    0.31162698639629240760918445839735,
                    0.0020127149292620447855048635533134
                    ),
                _s21(0.49868389831205793232488510274182),
                _s111(
                    0.4560656275499565553771507565422,
                    0.0098798242146205913479393294910199
                    ),
                _s21(0.29647986751297095527472829692039),
                _s111(
                    0.07405907467985092181313317616693,
                    0.0022261525470100421399714648358003
                    ),
                ])
            self.weights = numpy.concatenate([
                0.0033286951747132279123492864983453 * numpy.ones(3),
                0.00034732815949247716775761110709498 * numpy.ones(6),
                0.00030005010190990768147692126332646 * numpy.ones(6),
                0.0027304020592297865721122567792925 * numpy.ones(6),
                0.0011329714517846248235723081876289 * numpy.ones(6),
                0.00013652184099670511056952042535743 * numpy.ones(6),
                0.00076937606678432250468879499764765 * numpy.ones(6),
                0.0012969731946121357422045151897163 * numpy.ones(6),
                0.000036593224246369764117276863414231 * numpy.ones(3),
                0.0028058551621104313785063896117009 * numpy.ones(6),
                0.0017088242229934075069894095813209 * numpy.ones(6),
                0.0014835269432186778403695422530588 * numpy.ones(6),
                0.0016513963572791381733831138207843 * numpy.ones(3),
                0.0042541944332406005438980600721035 * numpy.ones(6),
                0.0021977995867183790800425034996552 * numpy.ones(6),
                0.0015485819478916375812737756528323 * numpy.ones(6),
                0.0026032805301535696135328744321822 * numpy.ones(6),
                0.00046617026410325438787372176076813 * numpy.ones(6),
                0.00065757715637351260925303268337316 * numpy.ones(6),
                0.0027864654094993725771943691572702 * numpy.ones(6),
                0.004259096678191470414522102679138 * numpy.ones(3),
                0.0041296569531416851894730156362364 * numpy.ones(6),
                0.0021542829173274647003494295251102 * numpy.ones(6),
                0.003475952201671089115896281707904 * numpy.ones(6),
                0.0024242799574098931250056610943746 * numpy.ones(3),
                0.0027066142543970948994776510136556 * numpy.ones(6),
                0.0003857792869821268431829777993308 * numpy.ones(6),
                0.0034396664292256190303443713927315 * numpy.ones(6),
                0.00043802693709037091369961704263301 * numpy.ones(6),
                0.0044060145278102687123499153063631 * numpy.ones(6),
                0.0020175008084553776352123506022532 * numpy.ones(6),
                0.0032382786045496547235654495167762 * numpy.ones(6),
                0.0046662107991939713663350642014653 * numpy.ones(6),
                0.0031592258177548723781412807006987 * numpy.ones(6),
                0.0034735845184462476910072280690917 * numpy.ones(6),
                0.00089539288040979812704666584412168 * numpy.ones(6),
                0.001080478483625377330560219489673 * numpy.ones(6),
                0.0018129744491099434327540980517936 * numpy.ones(6),
                0.00076515713048633167868677672775288 * numpy.ones(6),
                0.0028576817042790386648379351112937 * numpy.ones(6),
                0.0026558120905516908801336776112557 * numpy.ones(6),
                0.0044774376211064567033958341483738 * numpy.ones(6),
                0.0016535172612734334805128299877321 * numpy.ones(6),
                0.0019449896357145405264029362687799 * numpy.ones(6),
                0.0040545175807175370268042844227495 * numpy.ones(6),
                0.0023471970098679356791737946264535 * numpy.ones(6),
                0.005419796681891528521952310932314 * numpy.ones(6),
                0.00082141911555617243433300051895105 * numpy.ones(6),
                0.0013331230747819641329418490748497 * numpy.ones(6),
                0.0046726023038406813431924015852921 * numpy.ones(6),
                0.0035827112110621060024809133149837 * numpy.ones(3),
                0.0010473492443510246753163540590247 * numpy.ones(6),
                0.0027795988004446621224959331897008 * numpy.ones(3),
                0.0025570128208941351724792655059087 * numpy.ones(6),
                0.0011537659281314454307887463422756 * numpy.ones(6),
                0.0042237575830890532674535071698789 * numpy.ones(6),
                0.0013374175145030963635916919741648 * numpy.ones(6),
                0.0023814233045742072462775016255997 * numpy.ones(6),
                0.0042828524520190392122940883462578 * numpy.ones(6),
                0.00029995218691480666052959252395431 * numpy.ones(3),
                0.00044069395853764165813923075585858 * numpy.ones(6),
                0.0020678716670797737633981001897741 * numpy.ones(3),
                0.00011564525775079982244294178161948 * numpy.ones(6),
                0.00018242510180398216544034175162022 * numpy.ones(6),
                0.005135115847268768004534510408375 * numpy.ones(6),
                0.0039050557101546513702065538883118 * numpy.ones(6),
                0.00052821755983390689324171425704485 * numpy.ones(6),
                0.0018990634758470579055735050333984 * numpy.ones(6),
                0.0040277109885183108159818424750644 * numpy.ones(3),
                0.0023340593121636970794674024480338 * numpy.ones(6),
                0.00048267067797968960872305841256777 * numpy.ones(6),
                0.0045706738891261794937934566702863 * numpy.ones(3),
                0.0032490194493213969887464643881846 * numpy.ones(6),
                0.0013197669167676519743129401063918 * numpy.ones(6),
                0.005165355927453872464227868822851 * numpy.ones(3),
                0.0034852633534526555897586509106253 * numpy.ones(3),
                0.0003967292610444128373299498762537 * numpy.ones(6),
                0.0058568144539184509563857892539525 * numpy.ones(3),
                0.0011989487838384482035290357824589 * numpy.ones(3),
                0.00058559344102113444885270994252071 * numpy.ones(6),
                0.00074932205410521045063207466063199 * numpy.ones(3),
                0.0014007748847511846907649392824895 * numpy.ones(6),
                0.0057736389035923078386875324785005 * numpy.ones(3),
                0.00037050554881417209197128012650276 * numpy.ones(6),
                ])

        self.degree = index
        self.points = bary[:, [1, 2]]
        return

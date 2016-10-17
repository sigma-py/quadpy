# -*- coding: utf-8 -*-
#
import math
import numpy
import sympy

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


def show(tet, scheme, ball_scale=1.0, alpha=0.3):
    '''Shows the quadrature points on a given tetrahedron. The size of the
    balls around the points coincides with their weights.
    '''
    from matplotlib import pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.set_aspect('equal')

    edges = numpy.array([
        [tet[0], tet[1]],
        [tet[0], tet[2]],
        [tet[0], tet[3]],
        [tet[1], tet[2]],
        [tet[1], tet[3]],
        [tet[2], tet[3]],
        ])
    for edge in edges:
        plt.plot(edge[:, 0], edge[:, 1], edge[:, 2], '-k')

    transformed_pts = \
        + numpy.outer(
            (1.0 - scheme.points[:, 0]
                 - scheme.points[:, 1]
                 - scheme.points[:, 2]),
            tet[0]
            ) \
        + numpy.outer(scheme.points[:, 0], tet[1]) \
        + numpy.outer(scheme.points[:, 1], tet[2]) \
        + numpy.outer(scheme.points[:, 2], tet[3])

    tet_vol = volume(tet)
    phi, theta = numpy.mgrid[0:numpy.pi:101j, 0:2*numpy.pi:101j]
    x = numpy.sin(phi)*numpy.cos(theta)
    y = numpy.sin(phi)*numpy.sin(theta)
    z = numpy.cos(phi)
    for tp, weight in zip(transformed_pts, scheme.weights):
        color = 'b' if weight >= 0 else 'r'
        # highlight ball center
        plt.plot([tp[0]], [tp[1]], [tp[2]], '.' + color)

        # plot ball
        # scale the circle volume according to the weight
        r = ball_scale \
            * (tet_vol * abs(weight) / (4.0/3.0 * numpy.pi))**(1.0/3.0)

        ax.plot_surface(
            r*x + tp[0], r*y + tp[1], r*z + tp[2],
            color=color,
            alpha=alpha,
            linewidth=0
            )

    # http://stackoverflow.com/a/21765085/353337
    alpha = 1.3
    max_range = alpha * 0.5 * numpy.array([
        tet[:, 0].max() - tet[:, 0].min(),
        tet[:, 1].max() - tet[:, 1].min(),
        tet[:, 2].max() - tet[:, 2].min(),
        ]).max()
    mid_x = 0.5 * (tet[:, 0].max() + tet[:, 0].min())
    mid_y = 0.5 * (tet[:, 1].max() + tet[:, 1].min())
    mid_z = 0.5 * (tet[:, 2].max() + tet[:, 2].min())
    ax.set_xlim(mid_x - max_range, mid_x + max_range)
    ax.set_ylim(mid_y - max_range, mid_y + max_range)
    ax.set_zlim(mid_z - max_range, mid_z + max_range)
    return


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


def integrate(f, tetrahedron, scheme):
    g = _transform_to_unit_tetrahedron(f, tetrahedron)
    out = math.fsum([
        weight * g(point)
        for point, weight in zip(scheme.points, scheme.weights)
        ])
    return volume(tetrahedron) * out


def _s4():
    return numpy.array([
        [0.25, 0.25, 0.25, 0.25]
        ])


def _s31(a):
    b = 1.0 - 3*a
    return numpy.array([
        [a, a, a, b],
        [a, a, b, a],
        [a, b, a, a],
        [b, a, a, a],
        ])


def _s22(a):
    b = 0.5 - a
    return numpy.array([
        [a, a, b, b],
        [a, b, a, b],
        [b, a, a, b],
        [a, b, b, a],
        [b, a, b, a],
        [b, b, a, a],
        ])


def _s211(a, b):
    c = 1.0 - 2*a - b
    return numpy.array([
        [a, a, b, c],
        [a, b, a, c],
        [b, a, a, c],
        [a, b, c, a],
        [b, a, c, a],
        [b, c, a, a],
        [a, a, c, b],
        [a, c, a, b],
        [c, a, a, b],
        [a, c, b, a],
        [c, a, b, a],
        [c, b, a, a],
        ])


def _s1111(a, b, c):
    d = 1.0 - a - b - c
    return numpy.array([
        [a, b, c, d],
        [a, b, d, c],
        [a, c, b, d],
        [a, c, d, b],
        [a, d, b, c],
        [a, d, c, b],
        [b, a, c, d],
        [b, a, d, c],
        [b, c, a, d],
        [b, c, d, a],
        [b, d, a, c],
        [b, d, c, a],
        [c, a, b, d],
        [c, a, d, b],
        [c, b, a, d],
        [c, b, d, a],
        [c, d, a, b],
        [c, d, b, a],
        [d, a, b, c],
        [d, a, c, b],
        [d, b, a, c],
        [d, b, c, a],
        [d, c, a, b],
        [d, c, b, a],
        ])


class Keast(object):
    '''
    P. Keast,
    Moderate degree tetrahedral quadrature formulas,
    CMAME 55: 339-348
    (1986)

    Abstract:
    Quadrature formulas of degrees 4 to 8 for numerical integration over the
    tetrahedron are constructed. The formulas are fully symmetric with respect
    to the tetrahedron, and in some cases are the minimum point rules with this
    symmetry.

    https://people.sc.fsu.edu/~jburkardt/datasets/quadrature_rules_tet/quadrature_rules_tet.html
    '''
    def __init__(self, index):
        if index == 0:
            self.weights = numpy.array([
                1.0
                ])
            bary = _s4()
            self.degree = 1
        elif index == 1:
            self.weights = 0.25 * numpy.ones(4)
            bary = _s31(0.1381966011250105)
            self.degree = 2
        elif index == 2:
            self.weights = numpy.concatenate([
                -0.8 * numpy.ones(1),
                0.45 * numpy.ones(4),
                ])
            bary = numpy.concatenate([
                _s4(),
                _s31(1.0/6.0),
                ])
            self.degree = 3
        elif index == 3:
            self.weights = numpy.concatenate([
                0.2177650698804054 * numpy.ones(4),
                0.0214899534130631 * numpy.ones(6),
                ])
            bary = numpy.concatenate([
                _s31(0.1438564719343852),
                _s22(0.5),
                ])
            self.degree = 3
        elif index == 4:
            self.weights = numpy.concatenate([
                -148.0 / 1875.0 * numpy.ones(1),
                343.0 / 7500.0 * numpy.ones(4),
                56.0 / 375.0 * numpy.ones(6),
                ])
            bary = numpy.concatenate([
                _s4(),
                _s31(1.0/14.0),
                _s22(0.3994035761667992),
                ])
            self.degree = 4
        elif index == 5:
            self.weights = numpy.concatenate([
                2.0/105.0 * numpy.ones(6),
                0.0885898247429807 * numpy.ones(4),
                0.1328387466855907 * numpy.ones(4),
                ])
            bary = numpy.concatenate([
                _s22(0.5),
                _s31(0.1005267652252045),
                _s31(0.3143728734931922),
                ])
            self.degree = 4
        elif index == 6:
            self.weights = numpy.concatenate([
                6544.0 / 36015.0 * numpy.ones(1),
                81.0 / 2240.0 * numpy.ones(4),
                161051.0 / 2304960.0 * numpy.ones(4),
                338.0 / 5145.0 * numpy.ones(6),
                ])
            bary = numpy.concatenate([
                _s4(),
                _s31(1.0/3.0),
                _s31(1.0/11.0),
                _s22(0.0665501535736643),
                ])
            self.degree = 5
        elif index == 7:
            self.weights = numpy.concatenate([
                0.0399227502581679 * numpy.ones(4),
                0.0100772110553207 * numpy.ones(4),
                0.0553571815436544 * numpy.ones(4),
                27.0/560.0 * numpy.ones(12),
                ])
            bary = numpy.concatenate([
                _s31(0.2146028712591517),
                _s31(0.0406739585346113),
                _s31(0.3223378901422757),
                _s211(0.0636610018750175, 0.2696723314583159)
                ])
            self.degree = 6
        elif index == 8:
            self.weights = numpy.concatenate([
                0.1095853407966528 * numpy.ones(1),
                0.0635996491464850 * numpy.ones(4),
                -0.3751064406859797 * numpy.ones(4),
                0.0293485515784412 * numpy.ones(4),
                0.0058201058201058 * numpy.ones(6),
                0.1653439153439105 * numpy.ones(12),
                ])
            bary = numpy.concatenate([
                _s4(),
                _s31(0.0782131923303186),
                _s31(0.1218432166639044),
                _s31(0.3325391644464206),
                _s22(0.5),
                _s211(0.1, 0.2),
                ])
            self.degree = 7
        elif index == 9:
            self.weights = numpy.concatenate([
                -0.2359620398477557 * numpy.ones(1),
                0.0244878963560562 * numpy.ones(4),
                0.0039485206398261 * numpy.ones(4),
                0.0263055529507371 * numpy.ones(6),
                0.0829803830550589 * numpy.ones(6),
                0.0254426245481023 * numpy.ones(12),
                0.0134324384376852 * numpy.ones(12),
                ])
            bary = numpy.concatenate([
                _s4(),
                _s31(0.1274709365666390),
                _s31(0.0320788303926323),
                _s22(0.0497770956432810),
                _s22(0.1837304473985499),
                _s211(0.2319010893971509, 0.5132800333608811),
                _s211(0.0379700484718286, 0.1937464752488044),
                ])
            self.degree = 7
        else:
            raise ValueError('Illegal Keast index')
        self.points = bary[:, 1:]
        return


def _newton_cotes(n, point_fun):
    '''
    Construction after

    P. Silvester,
    Symmetric quadrature formulae for simplexes
    Math. Comp., 24, 95-100 (1970),
    <http://www.ams.org/journals/mcom/1970-24-109/S0025-5718-1970-0258283-6/S0025-5718-1970-0258283-6.pdf>
    '''
    def get_poly(t, m, n):
        f = 1
        for k in range(m):
            f *= (t - point_fun(k, n)) / (point_fun(m, n) - point_fun(k, n))
        return f
    degree = n
    num_points = (n+1) * (n**2 + 5*n+6) // 6
    bary = numpy.empty((num_points, 4))
    weights = numpy.empty(num_points)
    idx = 0
    for i in range(n + 1):
        for j in range(n + 1 - i):
            for k in range(n + 1 - i - j):
                l = n - i - j - k
                bary[idx] = point_fun(
                    numpy.array([i, j, k, l], dtype=float), n
                    )
                # Compute weight.
                # Define the polynomial which to integrate over the
                # tetrahedron.
                t = sympy.DeferredVector('t')
                g = sympy.expand(
                    get_poly(t[0], i, n)
                    * get_poly(t[1], j, n)
                    * get_poly(t[2], k, n)
                    * get_poly(t[3], l, n)
                    )
                # tranform it into a polynomial class
                gpoly = sympy.poly_from_expr(
                    g, (t[0], t[1], t[2], t[3])
                    )[0]
                # The integral of monomials over a tetrahedron are well-known,
                # see Silvester.
                weights[idx] = numpy.sum([
                     c * numpy.prod([math.factorial(k) for k in m]) * 6.0
                     / math.factorial(numpy.sum(m) + 3)
                     for m, c in zip(gpoly.monoms(), gpoly.coeffs())
                     ])
                idx += 1
    points = bary[:, [1, 2, 3]]
    return points, weights, degree


class NewtonCotesClosed(object):
    def __init__(self, n):
        self.points, self.weights, self.degree = \
            _newton_cotes(n, lambda k, n: k / float(n))
        return


class NewtonCotesOpen(object):
    def __init__(self, n):
        self.points, self.weights, self.degree = \
            _newton_cotes(n, lambda k, n: (k+1) / float(n+4))
        return


class Zienkiewicz(object):
    '''
    Olgierd Zienkiewicz,
    The Finite Element Method,
    Sixth Edition,
    Butterworth-Heinemann, 2005,
    ISBN: 0750663200,
    https://people.sc.fsu.edu/~jburkardt/datasets/quadrature_rules_tet/quadrature_rules_tet.html
    '''
    def __init__(self, index):
        if index == 4:
            self.weights = numpy.array([
                0.25,
                0.25,
                0.25,
                0.25
                ])
            self.points = numpy.array([
                [0.5854101966249685, 0.1381966011250105, 0.1381966011250105],
                [0.1381966011250105, 0.5854101966249685, 0.1381966011250105],
                [0.1381966011250105, 0.1381966011250105, 0.5854101966249685],
                [0.1381966011250105, 0.1381966011250105, 0.1381966011250105],
                ])
            self.degree = 2
        elif index == 5:
            self.weights = numpy.array([
                -0.8,
                0.45,
                0.45,
                0.45,
                0.45,
                ])
            self.points = numpy.array([
                [0.25, 0.25, 0.25],
                [0.5, 1.0/6.0, 1.0/6.0],
                [1.0/6.0, 0.5, 1.0/6.0],
                [1.0/6.0, 1.0/6.0, 0.5],
                [1.0/6.0, 1.0/6.0, 1.0/6.0],
                ])
            self.degree = 3
        else:
            raise ValueError('Illegal closed Newton-Cotes index')


class ShunnHam(object):
    '''
    Lee Shunn, Frank Ham,
    Symmetric quadrature rules for tetrahedra based on a cubic
    close-packed lattice arrangement,
    Journal of Computational and Applied Mathematics,
    2012,
    <http://dx.doi.org/10.1016/j.cam.2012.03.032>.

    Abstract:
    A family of quadrature rules for integration over tetrahedral volumes is
    developed. The underlying structure of the rules is based on the cubic
    close-packed (CCP) lattice arrangement using 1, 4, 10, 20, 35, and 56
    quadrature points. The rules are characterized by rapid convergence,
    positive weights, and symmetry. Each rule is an optimal approximation in
    the sense that lower-order terms have zero contribution to the truncation
    error and the leading-order error term is minimized. Quadrature formulas up
    to order 9 are presented with relevant numerical examples.
    '''
    def __init__(self, index):
        if index == 1:
            self.weights = numpy.array([
                1.0
                ])
            bary = numpy.array([
                [0.25, 0.25, 0.25, 0.25]
                ])
            self.degree = 1
        elif index == 2:
            self.weights = 0.25 * numpy.ones(4)
            bary = _s31(0.1381966011250110)
            self.degree = 2
        elif index == 3:
            self.weights = numpy.concatenate([
                0.0476331348432089 * numpy.ones(4),
                0.1349112434378610 * numpy.ones(6),
                ])
            bary = numpy.concatenate([
                _s31(0.0738349017262234),
                _s22(0.0937556561159491),
                ])
            self.degree = 3
        elif index == 4:
            self.weights = numpy.concatenate([
                0.0070670747944695 * numpy.ones(4),
                0.0469986689718877 * numpy.ones(12),
                0.1019369182898680 * numpy.ones(4),
                ])
            bary = numpy.concatenate([
                _s31(0.0323525947272439),
                _s211(0.0603604415251421, 0.2626825838877790),
                _s31(0.3097693042728620),
                ])
            self.degree = 5
        elif index == 5:
            self.weights = numpy.concatenate([
                0.0021900463965388 * numpy.ones(4),
                0.0143395670177665 * numpy.ones(12),
                0.0250305395686746 * numpy.ones(6),
                0.0479839333057554 * numpy.ones(12),
                0.0931745731195340 * numpy.ones(1)
                ])
            bary = numpy.concatenate([
                _s31(0.0267367755543735),
                _s211(0.0391022406356488, 0.7477598884818090),
                _s22(0.0452454000155172),
                _s211(0.2232010379623150, 0.0504792790607720),
                numpy.array([[0.25, 0.25, 0.25, 0.25]]),
                ])
            self.degree = 6
        elif index == 6:
            self.weights = numpy.concatenate([
                0.0010373112336140 * numpy.ones(4),
                0.0096016645399480 * numpy.ones(12),
                0.0164493976798232 * numpy.ones(12),
                0.0153747766513310 * numpy.ones(12),
                0.0293520118375230 * numpy.ones(12),
                0.0366291366405108 * numpy.ones(4),
                ])
            bary = numpy.concatenate([
                _s31(0.0149520651530592),
                _s211(0.0340960211962615, 0.1518319491659370),
                _s211(0.0462051504150017, 0.5526556431060170),
                _s211(0.2281904610687610, 0.0055147549744775),
                _s211(0.3523052600879940, 0.0992057202494530),
                _s31(0.1344783347929940)
                ])
            self.degree = 7
        else:
            raise ValueError('Illegal Shunn-Ham index')

        self.points = bary[:, 1:]

        return


class ZhangCuiLiu(object):
    '''
    Linbo Zhang, Tao Cui and Hui Liu,
    A set of symmetric quadrature rules on triangles and tetrahedra,
    Journal of Computational Mathematics
    Vol. 27, No. 1 (January 2009), pp. 89-96.

    Abstract:
    We present a program for computing symmetric quadrature rules on triangles
    and tetrahedra. A set of rules are obtained by using this program.
    Quadrature rules up to order 21 on triangles and up to order 14 on
    tetrahedra have been obtained which are useful for use in finite element
    computations. All rules presented here have positive weights with points
    lying within the integration domain.
    '''
    def __init__(self, index):
        if index == 1:
            self.weights = numpy.concatenate([
                0.0063971477799023213214514203351730 * numpy.ones(4),
                0.0401904480209661724881611584798178 * numpy.ones(4),
                0.0243079755047703211748691087719226 * numpy.ones(4),
                0.0548588924136974404669241239903914 * numpy.ones(4),
                0.0357196122340991824649509689966176 * numpy.ones(6),
                0.0071831906978525394094511052198038 * numpy.ones(12),
                0.0163721819453191175409381397561191 * numpy.ones(12),
                ])
            bary = numpy.concatenate([
                _s31(.0396754230703899012650713295393895),
                _s31(.3144878006980963137841605626971483),
                _s31(.1019866930627033000000000000000000),
                _s31(.1842036969491915122759464173489092),
                _s22(.0634362877545398924051412387018983),
                _s211(
                    .0216901620677280048026624826249302,
                    .7199319220394659358894349533527348
                    ),
                _s211(
                    .2044800806367957142413355748727453,
                    .5805771901288092241753981713906204
                    ),
                ])
            self.degree = 8
        elif index == 2:
            self.weights = numpy.concatenate([
                .0040651136652707670436208836835636 * numpy.ones(4),
                .0022145385334455781437599569500071 * numpy.ones(4),
                .0058134382678884505495373338821455 * numpy.ones(4),
                .0196255433858357215975623333961715 * numpy.ones(4),
                .0003875737905908214364538721248394 * numpy.ones(4),
                .0116429719721770369855213401005552 * numpy.ones(12),
                .0052890429882817131317736883052856 * numpy.ones(12),
                .0018310854163600559376697823488069 * numpy.ones(12),
                .0082496473772146452067449669173660 * numpy.ones(12),
                .0030099245347082451376888748208987 * numpy.ones(24),
                .0008047165617367534636261808760312 * numpy.ones(24),
                .0029850412588493071187655692883922 * numpy.ones(24),
                .0056896002418760766963361477811973 * numpy.ones(24),
                .0041590865878545715670013980182613 * numpy.ones(24),
                .0007282389204572724356136429745654 * numpy.ones(24),
                .0054326500769958248216242340651926 * numpy.ones(24),
                ])
            bary = numpy.concatenate([
                _s31(.3272533625238485639093096692685289),
                _s31(.0447613044666850808837942096478842),
                _s31(.0861403311024363536537208740298857),
                _s31(.2087626425004322968265357083976176),
                _s31(.0141049738029209600635879152102928),
                _s211(
                    .1021653241807768123476692526982584,
                    .5739463675943338202814002893460107
                    ),
                _s211(
                    .4075700516600107157213295651301783,
                    .0922278701390201300000000000000000
                    ),
                _s211(
                    .0156640007402803585557586709578084,
                    .7012810959589440327139967673208426
                    ),
                _s211(
                    .2254963562525029053780724154201103,
                    .4769063974420887115860583354107011
                    ),
                _s1111(
                    .3905984281281458000000000000000000,
                    .2013590544123922168123077327235092,
                    .0161122880710300298578026931548371
                    ),
                _s1111(
                    .1061350679989021455556139029848079,
                    .0327358186817269284944004077912660,
                    .0035979076537271666907971523385925
                    ),
                _s1111(
                    .5636383731697743896896816630648502,
                    .2302920722300657454502526874135652,
                    .1907199341743551862712487790637898
                    ),
                _s1111(
                    .3676255095325860844092206775991167,
                    .2078851380230044950717102125250735,
                    .3312104885193449000000000000000000
                    ),
                _s1111(
                    .7192323689817295295023401840796991,
                    .1763279118019329762157993033636973,
                    .0207602362571310090754973440611644
                    ),
                _s1111(
                    .5278249952152987298409240075817276,
                    .4372890892203418165526238760841918,
                    .0092201651856641949463177554949220
                    ),
                _s1111(
                    .5483674544948190728994910505607746,
                    .3447815506171641228703671870920331,
                    .0867217283322215394629438740085828
                    ),
                ])
            self.degree = 14
        else:
            raise ValueError('Illegal Zhang index')

        self.points = bary[:, [1, 2, 3]]
        return



class Yu(object):
    '''
    Yu Jinyun,
    Symmetyric Gaussian quadrature formulae for tetrahedronal regions,
    Computer Methods in Applied Mechanics and Engineering, 43 (1984) 349-353.

    Abstract:
    Quadrature formulae of degrees 2 to 6 are presented for the numerical
    integration of a function over tetrahedronal regions. The formulae
    presented are of Gaussian type and fully symmetric with respect to the four
    vertices of the tetrahedron.
    '''
    def __init__(self, index):
        if index == 1:
            self.weights = 0.25 * numpy.ones(4)
            bary = _s31(0.138196601125015)
            self.degree = 2
        elif index == 2:
            self.weights = numpy.concatenate([
                -0.8 * numpy.ones(1),
                0.45 * numpy.ones(4)
                ])
            bary = numpy.concatenate([
                _s4(),
                _s31(1.0/6.0)
                ])
            self.degree = 3
        elif index == 3:
            self.weights = numpy.concatenate([
                0.5037379410012282E-01 * numpy.ones(4),
                0.6654206863329239E-01 * numpy.ones(12)
                ])
            bary = numpy.concatenate([
                _s31(0.7611903264425430E-01),
                _s211(0.4042339134672644, 0.1197005277978019)
                ])
            self.degree = 4
        elif index == 4:
            self.weights = numpy.concatenate([
                0.1884185567365411 * numpy.ones(1),
                0.6703858372604275E-01 * numpy.ones(4),
                0.4528559236327399E-01 * numpy.ones(12)
                ])
            bary = numpy.concatenate([
                _s4(),
                _s31(0.8945436401412733E-01),
                _s211(0.4214394310662522, 0.1325810999384657),
                ])
            self.degree = 5
        elif index == 5:
            self.weights = numpy.concatenate([
                0.9040129046014750E-01 * numpy.ones(1),
                0.1911983427899124E-01 * numpy.ones(4),
                0.4361493840666568E-01 * numpy.ones(12),
                0.2581167596199161E-01 * numpy.ones(12)
                ])
            bary = numpy.concatenate([
                _s4(),
                _s31(0.5742691731735682E-01),
                _s211(0.2312985436519147, 0.5135188412556341E-01),
                _s211(0.4756909881472290E-01, 0.2967538129690260),
                ])
            self.degree = 6
        else:
            raise ValueError('Illegal closed Yu index')

        self.points = bary[:, [1, 2, 3]]
        return

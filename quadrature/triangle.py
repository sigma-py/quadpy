# -*- coding: utf-8 -*-
#
import math
import numpy


def volume(triangle):
    # It doesn't matter much which cross product we take for computing the
    # triangle volumes; deliberately take
    #
    #   <e0 x e1, e0 x e1> = <e0, e0> <e1, e1> - <e0, e1>^2.
    #
    e0 = triangle[1] - triangle[0]
    e1 = triangle[2] - triangle[1]
    e0_dot_e0 = numpy.dot(e0, e0)
    e0_dot_e1 = numpy.dot(e0, e1)
    e1_dot_e1 = numpy.dot(e1, e1)
    return 0.5 * numpy.sqrt(e0_dot_e0 * e1_dot_e1 - e0_dot_e1**2)


def show(triangle, scheme, circle_scale=1.0):
    '''Shows the quadrature points on a given triangle. The size of the circles
    around the points coincides with their weights.
    '''
    from matplotlib import pyplot as plt

    plt.plot(triangle[:, 0], triangle[:, 1], '-k')
    plt.plot(
        [triangle[-1, 0], triangle[0, 0]],
        [triangle[-1, 1], triangle[0, 1]],
        '-k')

    transformed_pts = \
        + numpy.outer(
            (1.0 - scheme.points[:, 0] - scheme.points[:, 1]),
            triangle[0]
            ) \
        + numpy.outer(scheme.points[:, 0], triangle[1]) \
        + numpy.outer(scheme.points[:, 1], triangle[2])

    # plt.plot(transformed_pts[:, 0], transformed_pts[:, 1], 'or')
    triangle_vol = volume(triangle)
    for tp, weight in zip(transformed_pts, scheme.weights):
        color = 'b' if weight >= 0 else 'r'
        # highlight circle center
        plt.plot([tp[0]], [tp[1]], '.' + color)
        # plot circle
        # scale the circle volume according to the weight
        radius = circle_scale \
            * numpy.sqrt(triangle_vol * abs(weight) / numpy.pi)
        circ = plt.Circle((tp[0], tp[1]), radius, color=color, alpha=0.5)
        plt.gcf().gca().add_artist(circ)

    plt.axis('equal')
    return


def _transform_to_unit_triangle(f, triangle):
    '''Transformation

      x = x0 * N0(xi, eta) + x1 * N1(xi, eta) + x2 * N2(xi, eta)

    with

      N0(xi, eta) = 1 - xi - eta,
      N1(xi, eta) = xi,
      N2(xi, eta) = eta.
    '''
    return lambda xi: f(
        + triangle[0] * (1.0 - xi[0] - xi[1])
        + triangle[1] * xi[0]
        + triangle[2] * xi[1]
        )


def integrate(f, triangle, rule):
    # w * f(x(xi)) * |det(J)|
    g = _transform_to_unit_triangle(f, triangle)
    out = math.fsum([
        weight * g(point)
        for point, weight in zip(rule.points, rule.weights)
        ])
    return volume(triangle) * out


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
        self.weights = [1.0]
        self.points = numpy.array([
            [1.0/3.0, 1.0/3.0]
            ])
        self.degree = 1
        return


class Vertex(object):
    def __init__(self):
        self.weights = numpy.array([1.0/3.0, 1.0/3.0, 1.0/3.0])
        self.points = _s21(0.0)
        self.degree = 1
        return


class SevenPoint(object):
    def __init__(self):
        self.weights = numpy.concatenate([
            0.45 * numpy.ones(1),
            0.05 * numpy.ones(3),
            2.0 / 15.0 * numpy.ones(3),
            ])
        self.points = numpy.concatenate([
            _s3(),
            _s21(0.0),
            _s21(0.5),
            ])
        self.degree = 3
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
    LC: TA335.S77.
    '''
    def __init__(self, index):
        if index == 1:
            self.weights = numpy.array([1.0/3.0, 1.0/3.0, 1.0/3.0])
            self.points = _s21(1.0/6.0)
            self.degree = 2
        elif index == 2:
            self.weights = numpy.array([1.0/3.0, 1.0/3.0, 1.0/3.0])
            self.points = _s21(0.5)
            self.degree = 2
        elif index == 3:
            self.weights = numpy.concatenate([
                -0.5625 * numpy.ones(1),
                25.0 / 48.0 * numpy.ones(3),
                ])
            self.points = numpy.concatenate([
                _s3(),
                _s21(0.2),
                ])
            self.degree = 3
        elif index == 4:
            self.weights = 1.0/6.0 * numpy.ones(6)
            self.points = _s111(0.659027622374092, 0.231933368553031)
            self.degree = 3
        elif index == 5:
            self.weights = numpy.concatenate([
                0.109951743655322 * numpy.ones(3),
                0.223381589678011 * numpy.ones(3),
                ])
            self.points = numpy.concatenate([
                _s21(0.091576213509771),
                _s21(0.445948490915965),
                ])
            self.degree = 4
        elif index == 6:
            self.weights = numpy.concatenate([
                0.375 * numpy.ones(1),
                5.0 / 48.0 * numpy.ones(6),
                ])
            self.points = numpy.concatenate([
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
            self.points = numpy.concatenate([
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
            self.points = numpy.concatenate([
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
            self.points = numpy.concatenate([
                _s21(0.063089014491502),
                _s21(0.249286745170910),
                _s111(0.636502499121399, 0.310352451033785),
                ])
            self.degree = 6
        elif index == 10:
            self.weights = numpy.concatenate([
                -0.149570044467670 * numpy.ones(1),
                0.175615257433204 * numpy.ones(3),
                0.053347235608839 * numpy.ones(3),
                0.077113760890257 * numpy.ones(6),
                ])
            self.points = numpy.concatenate([
                _s3(),
                _s21(0.260345966079038),
                _s21(0.065130102902216),
                _s111(0.638444188569809, 0.312865496004875),
                ])
            self.degree = 7
        else:
            raise ValueError('Illegal Strang index')

        return


class Toms584_19(object):
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
        self.weights = numpy.concatenate([
            0.0378610912003147 * numpy.ones(1),
            0.0376204254131829 * numpy.ones(3),
            0.0783573522441174 * numpy.ones(3),
            0.1162714796569659 * numpy.ones(3),
            0.0134442673751655 * numpy.ones(3),
            0.0375097224552317 * numpy.ones(6),
            ])

        self.points = numpy.concatenate([
            _s3(),
            _s21(0.1012865073234563),
            _s21(0.4701420641051151),
            _s21(0.2321023267750504),
            _s21(0.0294808608844396),
            _s111(0.7384168123405100, 0.2321023267750504),
            ])
        self.degree = 8
        return


class Toms612_19(object):
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
    def __init__(self):
        self.weights = numpy.concatenate([
            9.71357962827961025E-002 * numpy.ones(1),
            3.13347002271398278E-002 * numpy.ones(3),
            7.78275410047754301E-002 * numpy.ones(3),
            7.96477389272090969E-002 * numpy.ones(3),
            2.55776756586981006E-002 * numpy.ones(3),
            4.32835393772893970E-002 * numpy.ones(6),
            ])
        self.points = numpy.concatenate([
            _s3(),
            _s21(0.48968251919873701),
            _s21(0.43708959149293553),
            _s21(0.18820353561903219),
            _s21(4.47295133944529688E-002),
            _s111(0.74119859878449801, 3.68384120547362581E-002),
            ])
        self.degree = 9
        return


class Toms612_28(object):
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
    def __init__(self):
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
        self.points = numpy.concatenate([
            _s3(),
            _s21(0.02598914092828833),
            _s21(0.09428750264792270),
            _s21(0.4946367750172147),
            _s21(0.2073433826145142),
            _s21(0.4389078057004907),
            _s111(0.6779376548825902, 0.04484167758913055),
            _s111(0.8588702812826364, 0.0),
            ])
        self.degree = 11
        return


class Toms706_37(object):
    '''
    See
    https://people.sc.fsu.edu/~jburkardt/datasets/quadrature_rules_tri/quadrature_rules_tri.html
    and

    Berntsen and Espelid,
    Algorithm 706: DCUTRI: An Algorithm for Adaptive Cubature over a Collection
    of Triangles,
    ACM Trans. Math. Softw.,
    Sept. 1992,
    10.1145/131766.131772,
    <http://dl.acm.org/citation.cfm?id=131772>.
    '''
    def __init__(self):
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
        self.points = numpy.concatenate([
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
                0.635867859433872768286976979827
                ),
            ])
        self.degree = 13
        return


class Dunavant(object):
    '''
    Triangle integration schemes from

    D. A. Dunavant,
    High Degree Efficient Symmetrical Gaussian Quadrature Rules for the
    Triangle,
    Article in International Journal for Numerical Methods in Engineering,
    21(6):1129-1148, June 1985,
    10.1002/nme.1620210612.
    '''
    def __init__(self, index):
        if index == 1:
            self.weights = [1.0]
            bary = numpy.array([
                [1.0/3.0, 1.0/3.0, 1.0/3.0]
                ])
            self.degree = 1
        elif index == 2:
            self.weights = 1.0/3.0 * numpy.ones(3)
            bary = self.double_mix(2.0/3.0, 1.0/6.0)
            self.degree = 2
        elif index == 3:
            self.weights = numpy.concatenate([
                numpy.array([-0.5625]),
                25.0 / 48.0 * numpy.ones(3),
                ])
            bary = numpy.concatenate([
                numpy.array([[1.0/3.0, 1.0/3.0, 1.0/3.0]]),
                self.double_mix(0.6, 0.2),
                ])
            self.degree = 3
        elif index == 4:
            self.weights = numpy.concatenate([
                0.223381589678011 * numpy.ones(3),
                0.109951743655322 * numpy.ones(3),
                ])
            bary = numpy.concatenate([
                self.double_mix(0.108103018168070, 0.445948490915965),
                self.double_mix(0.816847572980459, 0.091576213509771),
                ])
            self.degree = 4
        elif index == 5:
            self.weights = numpy.concatenate([
                0.225 * numpy.ones(1),
                0.132394152788506 * numpy.ones(3),
                0.125939180544827 * numpy.ones(3),
                ])
            bary = numpy.concatenate([
                numpy.array([[1.0/3.0, 1.0/3.0, 1.0/3.0]]),
                self.double_mix(0.059715871789770, 0.4701420641051),
                self.double_mix(0.797426985353087, 0.101286507323456),
                ])
            self.degree = 5
        elif index == 6:
            self.weights = numpy.concatenate([
                0.116786275726379 * numpy.ones(3),
                0.050844906370207 * numpy.ones(3),
                0.082851075618374 * numpy.ones(6),
                ])
            bary = numpy.concatenate([
                self.double_mix(0.501426509658179, 0.249286745170910),
                self.double_mix(0.873821971016996, 0.063089014491502),
                self.triple_mix(
                    0.053145049844817, 0.310352451033784, 0.636502499121399
                    ),
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
                numpy.array([[1.0/3.0, 1.0/3.0, 1.0/3.0]]),
                self.double_mix(0.479308067841920, 0.260345966079040),
                self.double_mix(0.869739794195568, 0.065130102902216),
                self.triple_mix(
                    0.048690315425316, 0.312865496004874, 0.638444188569810
                    ),
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
                numpy.array([[1.0/3.0, 1.0/3.0, 1.0/3.0]]),
                self.double_mix(0.081414823414554, 0.459292588292723),
                self.double_mix(0.658861384496480, 0.170569307751760),
                self.double_mix(0.898905543365938, 0.050547228317031),
                self.triple_mix(
                    0.008394777409958, 0.263112829634638, 0.728492392955404
                    ),
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
                numpy.array([[1.0/3.0, 1.0/3.0, 1.0/3.0]]),
                self.double_mix(0.020634961602525, 0.489682519198738),
                self.double_mix(0.125820817014127, 0.437089591492937),
                self.double_mix(0.623592928761935, 0.188203535619033),
                self.double_mix(0.910540973211095, 0.044729513394453),
                self.triple_mix(
                    0.036838412054736, 0.221962989160766, 0.741198598784498
                    ),
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
                numpy.array([[1.0/3.0, 1.0/3.0, 1.0/3.0]]),
                self.double_mix(0.028844733232685, 0.485577633383657),
                self.double_mix(0.781036849029926, 0.109481575485037),
                self.triple_mix(
                    0.141707219414880, 0.307939838764121, 0.550352941820999
                    ),
                self.triple_mix(
                    0.025003534762686, 0.246672560639903, 0.728323904597411
                    ),
                self.triple_mix(
                    0.009540815400299, 0.066803251012200, 0.923655933587500
                    ),
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
                self.double_mix(-0.069222096541517, 0.534611048270758),
                self.double_mix(0.202061394068290, 0.398969302965855),
                self.double_mix(0.593380199137435, 0.203309900431282),
                self.double_mix(0.761298175434837, 0.119350912282581),
                self.double_mix(0.935270103777448, 0.032364948111276),
                self.triple_mix(
                    0.050178138310495, 0.356620648261293, 0.593201213428213
                    ),
                self.triple_mix(
                    0.021022016536166, 0.171488980304042, 0.807489003159792
                    ),
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
                self.double_mix(0.023565220452390, 0.488217389773805),
                self.double_mix(0.120551215411079, 0.439724392294460),
                self.double_mix(0.457579229975768, 0.271210385012116),
                self.double_mix(0.744847708916828, 0.127576145541586),
                self.double_mix(0.957365299093579, 0.021317350453210),
                self.triple_mix(
                    0.115343494534698, 0.275713269685514, 0.608943235779788
                    ),
                self.triple_mix(
                    0.022838332222257, 0.281325580989940, 0.695836086787803
                    ),
                self.triple_mix(
                    0.025734050548330, 0.116251915907597, 0.858014033544073
                    ),
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
                numpy.array([[1.0/3.0, 1.0/3.0, 1.0/3.0]]),
                self.double_mix(0.009903630120591, 0.495048184939705),
                self.double_mix(0.062566729780852, 0.468716635109574),
                self.double_mix(0.170957326397447, 0.414521336801277),
                self.double_mix(0.541200855914337, 0.229399572042831),
                self.double_mix(0.771151009607340, 0.114424495196330),
                self.double_mix(0.950377217273082, 0.024811391363459),
                self.triple_mix(
                    0.094853828379579, 0.268794997058761, 0.636351174561660
                    ),
                self.triple_mix(
                    0.018100773278807, 0.291730066734288, 0.690169159986905
                    ),
                self.triple_mix(
                    0.022233076674090, 0.126357385491669, 0.851409537834241
                    ),
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
                self.double_mix(0.022072179275643, 0.488963910362179),
                self.double_mix(0.164710561319092, 0.417644719340454),
                self.double_mix(0.453044943382323, 0.273477528308839),
                self.double_mix(0.645588935174913, 0.177205532412543),
                self.double_mix(0.876400233818255, 0.061799883090873),
                self.double_mix(0.961218077502598, 0.019390961248701),
                self.triple_mix(
                    0.057124757403648, 0.172266687821356, 0.770608554774996
                    ),
                self.triple_mix(
                    0.092916249356972, 0.336861459796345, 0.570222290846683
                    ),
                self.triple_mix(
                    0.014646950055654, 0.298372882136258, 0.686980167808088
                    ),
                self.triple_mix(
                    0.001268330932872, 0.118974497696957, 0.879757171370171
                    ),
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
                self.double_mix(-0.013945833716486, 0.506972916858243),
                self.double_mix(0.137187291433955, 0.431406354283023),
                self.double_mix(0.444612710305711, 0.277693644847144),
                self.double_mix(0.747070217917492, 0.126464891041254),
                self.double_mix(0.858383228050628, 0.070808385974686),
                self.double_mix(0.962069659517853, 0.018965170241073),
                self.triple_mix(
                    0.133734161966621, 0.261311371140087, 0.604954466893291
                    ),
                self.triple_mix(
                    0.036366677396917, 0.575586555512814, 0.388046767090269
                    ),
                self.triple_mix(
                    -0.010174883126571, 0.285712220049916, 0.724462663076655
                    ),
                self.triple_mix(
                    0.036843869875878, 0.215599664072284, 0.747556466051838
                    ),
                self.triple_mix(
                    0.012459809331199, 0.103575616576386, 0.883964574092416
                    ),
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
                numpy.array([[1.0/3.0, 1.0/3.0, 1.0/3.0]]),
                self.double_mix(0.005238916103123, 0.497380541948438),
                self.double_mix(0.173061122901295, 0.413469438549352),
                self.double_mix(0.059082801866017, 0.470458599066991),
                self.double_mix(0.518892500060958, 0.240553749969521),
                self.double_mix(0.704068411554854, 0.147965794222573),
                self.double_mix(0.849069624685052, 0.075465187657474),
                self.double_mix(0.966807194753950, 0.016596402623025),
                self.triple_mix(
                    0.103575692245252, 0.296555596579887, 0.599868711174861
                    ),
                self.triple_mix(
                    0.020083411655416, 0.337723063403079, 0.642193524941505
                    ),
                self.triple_mix(
                    -0.004341002614139, 0.204748281642812, 0.799592720971327
                    ),
                self.triple_mix(
                    0.041941786468010, 0.189358492130623, 0.768699721401368
                    ),
                self.triple_mix(
                    0.014317320230681, 0.085283615682657, 0.900399064086661
                    ),
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
                numpy.array([[1.0/3.0, 1.0/3.0, 1.0/3.0]]),
                self.double_mix(0.005658918886452, 0.497170540556774),
                self.double_mix(0.035647354750751, 0.482176322624625),
                self.double_mix(0.099520061958437, 0.450239969020782),
                self.double_mix(0.199467521245206, 0.400266239377397),
                self.double_mix(0.495717464058095, 0.252141267970953),
                self.double_mix(0.675905990683077, 0.162047004658461),
                self.double_mix(0.848248235478508, 0.075875882260746),
                self.double_mix(0.968690546064356, 0.015654726967822),
                self.triple_mix(
                    0.010186928826919, 0.334319867363658, 0.655493203809423
                    ),
                self.triple_mix(
                    0.135440871671036, 0.292221537796944, 0.572337590532020
                    ),
                self.triple_mix(
                    0.054423924290583, 0.319574885423190, 0.626001190286228
                    ),
                self.triple_mix(
                    0.012868560833637, 0.190704224192292, 0.796427214974071
                    ),
                self.triple_mix(
                    0.067165782413524, 0.180483211648746, 0.752351005937729
                    ),
                self.triple_mix(
                    0.014663182224828, 0.080711313679564, 0.904625504095608
                    ),
                ])
            self.degree = 17
        # elif index == 18:
        #     self.weights = numpy.concatenate([
        #         0.030809939937647 * numpy.ones(1),
        #         0.009072436679404 * numpy.ones(3),
        #         0.018761316939594 * numpy.ones(3),
        #         0.019441097985477 * numpy.ones(3),
        #         0.027753948610810 * numpy.ones(3),
        #         0.032256225351457 * numpy.ones(3),
        #         0.025074032616922 * numpy.ones(3),
        #         0.015271927971832 * numpy.ones(3),
        #         0.006793922022963 * numpy.ones(3),
        #         -0.002223098729920 * numpy.ones(3),
        #         0.006331914076406 * numpy.ones(6),
        #         0.027257538049138 * numpy.ones(6),
        #         0.017676785649465 * numpy.ones(6),
        #         0.018379484638070 * numpy.ones(6),
        #         0.008104732808192 * numpy.ones(6),
        #         0.007634129070725 * numpy.ones(6),
        #         0.000046187660794 * numpy.ones(6),
        #         ])
        #     bary = numpy.concatenate([
        #         numpy.array([[1.0/3.0, 1.0/3.0, 1.0/3.0]]),
        #         self.double_mix(0.013310382738157, 0.493344808630921),
        #         self.double_mix(0.061578811516086, 0.469210594241957),
        #         self.double_mix(0.127437208225989, 0.436261395887006),
        #         self.double_mix(0.210307658653168, 0.394846170673416),
        #         self.double_mix(0.500410862393686, 0.249794568803157),
        #         self.double_mix(0.677135612512315, 0.161432193743843),
        #         self.double_mix(0.846803545029257, 0.076598227485371),
        #         self.double_mix(0.951495121293100, 0.024252439353450),
        #         self.double_mix(0.913707265566071, 0.043146367216965),
        #         self.triple_mix(
        #             0.008430536202420, 0.358911494940944, 0.632657968856636
        #             ),
        #         self.triple_mix(
        #             0.131186551737188, 0.294402476751957, 0.574410971510855
        #             ),
        #         self.triple_mix(
        #             0.050203151565675, 0.325017801641814, 0.624779046792512
        #             ),
        #         self.triple_mix(
        #             0.066329263810916, 0.184737559666046, 0.748933176523037
        #             ),
        #         self.triple_mix(
        #             0.011996194566236, 0.218796800013321, 0.769207005420443
        #             ),
        #         self.triple_mix(
        #             0.014858100590125, 0.101179597136408, 0.883962302273467
        #             ),
        #         self.triple_mix(
        #             -0.035222015287949, 0.020874755282586, 1.014347260005363
        #             ),
        #         ])

        #     self.degree = 18
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
                numpy.array([[1.0/3.0, 1.0/3.0, 1.0/3.0]]),
                self.double_mix(0.020780025853987, 0.489609987073006),
                self.double_mix(0.090926214604215, 0.454536892697893),
                self.double_mix(0.197166638701138, 0.401416680649431),
                self.double_mix(0.488896691193805, 0.255551654403098),
                self.double_mix(0.645844115695741, 0.177077942152130),
                self.double_mix(0.779877893544096, 0.110061053227952),
                self.double_mix(0.888942751496321, 0.055528624251840),
                self.double_mix(0.974756272445543, 0.012621863777229),
                self.triple_mix(
                    0.003611417848412, 0.395754787356943, 0.600633794794645
                    ),
                self.triple_mix(
                    0.134466754530780, 0.307929983880436, 0.557603261588784
                    ),
                self.triple_mix(
                    0.014446025776115, 0.264566948406520, 0.720987025817365
                    ),
                self.triple_mix(
                    0.046933578838178, 0.358539352205951, 0.594527068955871
                    ),
                self.triple_mix(
                    0.002861120350567, 0.157807405968595, 0.839331473680839
                    ),
                self.triple_mix(
                    0.223861424097916, 0.075050596975911, 0.701087978926173
                    ),
                self.triple_mix(
                    0.034647074816760, 0.142421601113383, 0.822931324069857
                    ),
                self.triple_mix(
                    0.010161119296278, 0.065494628082938, 0.924344252620784
                    ),
                ])
            self.degree = 19
        elif index == 20:
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
                numpy.array([[1.0/3.0, 1.0/3.0, 1.0/3.0]]),
                self.double_mix(-0.001900928704400, 0.500950464352200),
                self.double_mix(0.023574084130543, 0.488212957934729),
                self.double_mix(0.089726636099435, 0.455136681950283),
                self.double_mix(0.196007481363421, 0.401996259318289),
                self.double_mix(0.488214180481157, 0.255892909759421),
                self.double_mix(0.647023488009788, 0.176488255995106),
                self.double_mix(0.791658289326483, 0.104170855336758),
                self.double_mix(0.893862072318140, 0.053068963840930),
                self.double_mix(0.916762569607942, 0.041618715196029),
                self.double_mix(0.976836157186356, 0.011581921406822),
                self.triple_mix(
                    0.048741583664839, 0.344855770229001, 0.606402646106160
                    ),
                self.triple_mix(
                    0.006314115948605, 0.377843269594854, 0.615842614456541
                    ),
                self.triple_mix(
                    0.134316520547348, 0.306635479062357, 0.559048000390295
                    ),
                self.triple_mix(
                    0.013973893962392, 0.249419362774742, 0.736606743262866
                    ),
                self.triple_mix(
                    0.075549132909764, 0.212775724802802, 0.711675142287434
                    ),
                self.triple_mix(
                    -0.008368153208227, 0.146965436053239, 0.861402717154987
                    ),
                self.triple_mix(
                    0.026686063258714, 0.137726978828923, 0.835586957912363
                    ),
                self.triple_mix(
                    0.010547719294141, 0.059696109149007, 0.929756171556853
                    ),
                ])
            self.degree = 20
        else:
            raise ValueError('Illegal Dunavant index')

        # convert barycentric coordinates to reference triangle
        self.points = bary[:, [1, 2]]
        return

    def double_mix(self, a, b):
        return numpy.array([
            [a, b, b],
            [b, a, b],
            [b, b, a],
            ])

    def triple_mix(self, a, b, c):
        return numpy.array([
            [a, b, c],
            [c, a, b],
            [b, c, a],
            [a, c, b],
            [b, a, c],
            [c, b, a],
            ])


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
        elif index == 3:
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
        else:
            raise ValueError('Illegal Zhang index')

        self.points = bary[:, [1, 2]]
        return


class WandzuraXiao(object):
    '''
    S. Wandzurat, H. Xiao,
    Symmetric quadrature rules on a triangle,
    Computers & Mathematics with Applications
    Volume 45, Issue 12, June 2003, Pages 1829-1840,
    doi:10.1016/S0898-1221(03)90004-6.

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
        elif index == 6:
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
        else:
            raise ValueError('Illegal Wandzura index')

        self.points = bary[:, [1, 2]]
        return


class LynessJespersen(object):
    '''
    J.N. Lyness, D. Jespersen,
    Moderate Degree Symmetric Quadrature Rules for the Triangle,
    J. Inst. Maths Applies (1975) 15, 19-32,
    doi: 10.1093/imamat/15.1.19.

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
        elif index == 21:
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
        else:
            raise ValueError('Illegal Lyness-Jespersen index')

        self.points = bary[:, [1, 2]]
        return

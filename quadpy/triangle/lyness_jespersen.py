# -*- coding: utf-8 -*-
#
from .helpers import _s3, _s21, _s111
import numpy


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

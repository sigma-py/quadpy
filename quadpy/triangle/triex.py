# -*- coding: utf-8 -*-
#
from .helpers import _s3, _s21, _s111
import numpy


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

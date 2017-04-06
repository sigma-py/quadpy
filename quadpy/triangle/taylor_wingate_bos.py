# -*- coding: utf-8 -*-
#
from .helpers import _s21, _s111
import numpy


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

# -*- coding: utf-8 -*-
#
from .helpers import _s3, _s21, _s111
import numpy


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

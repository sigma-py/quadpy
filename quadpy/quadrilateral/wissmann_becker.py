# -*- coding: utf-8 -*-
#
import numpy


class WissmannBecker(object):
    '''
    Johannes W. Wissmann and Thomas Becker,
    Partially Symmetric Cubature Formulas for Even Degrees of Exactness,
    SIAM J. Numer. Anal., 23(3), 676–685, 10 pages,
    <https://doi.org/10.1137/0723043>.
    '''
    def __init__(self, index):
        self.name = 'WB({})'.format(index)
        if index == '4-1':
            self.degree = 4
            self.weights = numpy.concatenate([
                numpy.full(1, 1.142857142857143),
                numpy.full(1, 0.439560439560440),
                numpy.full(2, 0.566072207007532),
                numpy.full(2, 0.642719001783677)
                ])
            self.points = numpy.concatenate([
                _z(0.0),
                _z(0.966091783079296),
                _m(0.851914653304601, 0.455603727836193),
                _m(0.630912788976754, -0.731629951573135)
                ])
        elif index == '4-2':
            self.degree = 4
            self.weights = numpy.concatenate([
                numpy.full(1, 1.286412084888852),
                numpy.full(1, 0.491365692888926),
                numpy.full(2, 0.761883709085613),
                numpy.full(2, 0.349227402025498)
                ])
            self.points = numpy.concatenate([
                _z(-0.356822089773090),
                _z(0.934172358962716),
                _m(0.774596669241483, 0.390885162530071),
                _m(0.774596669241483, -0.852765377881771)
                ])
        elif index == '6-1':
            self.degree = 6
            self.weights = numpy.concatenate([
                numpy.full(1, 0.455343245714174),
                numpy.full(1, 0.827395973202966),
                numpy.full(2, 0.144000884599645),
                numpy.full(2, 0.668259104262665),
                numpy.full(2, 0.225474004890679),
                numpy.full(2, 0.320896396788441)
                ])
            self.points = numpy.concatenate([
                _z(0.836405633697626),
                _z(-0.357460165391307),
                _m(0.888764014654765, 0.872101531193131),
                _m(0.604857639464685, 0.305985162155427),
                _m(0.955447506641064, -0.410270899466658),
                _m(0.565459993438754, -0.872869311156879)
                ])
        elif index == '6-2':
            self.degree = 6
            self.weights = numpy.concatenate([
                numpy.full(1, 0.392750590964348),
                numpy.full(1, 0.754762881242610),
                numpy.full(2, 0.206166050588279),
                numpy.full(2, 0.689992138489864),
                numpy.full(2, 0.260517488732317),
                numpy.full(2, 0.269567586086061),
                ])
            self.points = numpy.concatenate([
                _z(0.869833375250059),
                _z(-0.479406351612111),
                _m(0.863742826346154, 0.802837516207657),
                _m(0.518690521392582, 0.262143665508058),
                _m(0.933972544972849, -0.363096583148066),
                _m(0.608977536016356, -0.896608632762453)
                ])
        elif index == '8-1':
            self.degree = 8
            self.weights = numpy.concatenate([
                numpy.full(1, 0.055364705621440),
                numpy.full(1, 0.404389368726076),
                numpy.full(1, 0.533546604952635),
                numpy.full(1, 0.117054188786739),
                numpy.full(2, 0.125614417613747),
                numpy.full(2, 0.136544584733588),
                numpy.full(2, 0.483408479211257),
                numpy.full(2, 0.252528506429544),
                numpy.full(2, 0.361262323882172),
                numpy.full(2, 0.085464254086247)
                ])
            self.points = numpy.concatenate([
                _z(0.0),
                _z(0.757629177660505),
                _z(-0.236871842255702),
                _z(-0.989717929044527),
                _m(0.639091304900370, 0.950520955645667),
                _m(0.937069076924990, 0.663882736885633),
                _m(0.537083530541494, 0.304210681724104),
                _m(0.887188506449625, -0.236496718536120),
                _m(0.494698820670197, -0.698953476086564),
                _m(0.897495818279768, -0.900390774211580),
                ])
        else:
            assert index == '8-2'
            self.degree = 8
            self.weights = numpy.concatenate([
                numpy.full(1, 0.659560131960342),
                numpy.full(1, -0.949142923043125),
                numpy.full(2, 0.765051819557684),
                numpy.full(2, 0.936975981088416),
                numpy.full(2, 0.333656717735747),
                numpy.full(2, -0.079583272377397),
                numpy.full(2, -0.272240080612534),
                numpy.full(2, -0.613735353398028),
                numpy.full(2, -0.888477650535971),
                ])
            self.points = numpy.concatenate([
                _z(0.659560131960342),
                _z(-0.949142923043125),
                _m(0.952509466071562, 0.765051819557684),
                _m(0.532327454074206, 0.936975981088416),
                _m(0.684736297951735, 0.333656717735747),
                _m(0.233143240801405, -0.079583272377397),
                _m(0.927683319306117, -0.272240080612534),
                _m(0.453120687403749, -0.613735353398028),
                _m(0.837503640422812, -0.888477650535971),
                ])
        return


def _z(a):
    return numpy.array([
        [0.0, a]
        ])


def _m(a, b):
    return numpy.array([
        [+a, +b],
        [-a, +b],
        ])

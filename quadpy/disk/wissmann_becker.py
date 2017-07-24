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
        if index == '6-1':
            self.degree = 6
            self.weights = numpy.concatenate([
                numpy.full(1, 0.192106330409167),
                numpy.full(1, 0.642883172453235),
                numpy.full(2, 0.030327815305342),
                numpy.full(2, 0.502358408812314),
                numpy.full(2, 0.340540094500688),
                numpy.full(2, 0.280075256745352),
                ])
            self.points = numpy.concatenate([
                _z(0.914874261544932),
                _z(-0.206355543412215),
                _m(0.978839047064283, 0.662952453661044),
                _m(0.458236515708684, 0.464889889783278),
                _m(0.802017004072343, -0.174207378839556),
                _m(0.373683864304499, -0.770749853148807),
                ])
        elif index == '6-2':
            self.degree = 6
            self.weights = numpy.concatenate([
                numpy.full(1, 0.552687596464871),
                numpy.full(1, 0.076090439302545),
                numpy.full(1, 0.441690572122440),
                numpy.full(2, 0.076090439302545),
                numpy.full(2, 0.076090439302545),
                numpy.full(2, 0.441690572122440),
                numpy.full(2, 0.441690572122440),
                ])
            self.points = numpy.concatenate([
                _z(0.0),
                _z(1.032546063890842),
                _z(-0.726359846042944),
                _m(0.982009662438297, 0.319074281217230),
                _m(0.606915348667678, -0.835347313162651),
                _m(0.426943605361474, 0.587637459480312),
                _m(0.690809264754287, -0.224457536458840),
                ])
        else:
            assert index == '8-1'
            self.degree = 8
            self.weights = numpy.concatenate([
                numpy.full(1, 0.334521439965580),
                numpy.full(1, 0.087911387704639),
                numpy.full(1, 0.167628039106469),
                numpy.full(1, 0.305874815913735),
                numpy.full(2, 0.087911387704639),
                numpy.full(2, 0.087911387704639),
                numpy.full(2, 0.167628039106469),
                numpy.full(2, 0.167628039106469),
                numpy.full(2, 0.305874815913735),
                numpy.full(2, 0.305874815913735),
                ])
            self.points = numpy.concatenate([
                _z(0.000000000000000),
                _z(0.953321175521807),
                _z(-0.882458098900126),
                _z(0.582334188120499),
                _m(0.906662316102171, 0.294592444333741),
                _m(0.560348127669843, -0.771253032094644),
                _m(0.518695856299547, 0.713923598834010),
                _m(0.839267525316398, -0.272694549383947),
                _m(0.553832724273449, 0.179951160534772),
                _m(0.342287447682940, -0.471118254595022),
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

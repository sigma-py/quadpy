# -*- coding: utf-8 -*-
#
import numpy


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
                _symm(0.0, 0.5),
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
                2.0/3.0 - _symm(0.0, 0.5),
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
        #         _symm(0.591506351, 0.158493649),
        #         ])
        #     self.degree = 2
        # elif index == 5:
        #     self.weights = 2.0 * numpy.concatenate([
        #         1.0/18.0 * numpy.ones(1),
        #         2.0/9.0 * numpy.ones(2),
        #         ])
        #     self.points = numpy.concatenate([
        #         2.0/3.0 - numpy.array([[0.0, 0.0]]),
        #         2.0/3.0 - _symm(0.591506351, 0.158493649),
        #         ])
        #     self.degree = 2
        # elif index == 6:
        #     self.weights = 2.0 * numpy.concatenate([
        #         1.0/8.0 * numpy.ones(4),
        #         ])
        #     lmbda = 0.655308609
        #     mu = 0.247060398
        #     self.points = numpy.concatenate([
        #         _symm(lmbda, mu),
        #         2.0/3.0 - _symm(lmbda, mu),
        #         ])
        #     self.degree = 2
        # elif index == 7:
        #     self.weights = 2.0 * numpy.concatenate([
        #         0.159020691 * numpy.ones(2),
        #         0.090979309 * numpy.ones(2),
        #         ])
        #     self.points = numpy.concatenate([
        #         _symm(0.666390246, 0.280019915),
        #         _symm(0.178558728, 0.075031109),
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
        #         _symm(0.0, 0.8),
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
        #         _symm(0.2, 0.6),
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
        #         _symm(0.939332590, 0.0),
        #         _symm(0.0, 0.340667409),
        #         numpy.array([[0.4, 0.4]]),
        #         ])
        #     self.degree = 3

        return

def _symm(a, b):
    return numpy.array([
        [a, b],
        [b, a],
        ])

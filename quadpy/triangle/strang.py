# -*- coding: utf-8 -*-
#
from .helpers import _s3, _s21, _s111
import numpy


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

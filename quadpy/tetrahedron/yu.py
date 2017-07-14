# -*- coding: utf-8 -*-
#
import numpy

from .helpers import _s4, _s31, _s211


class Yu(object):
    '''
    Yu Jinyun,
    Symmetyric Gaussian quadrature formulae for tetrahedronal regions,
    Computer Methods in Applied Mechanics and Engineering, 43 (1984) 349-353,
    <https://dx.doi.org/10.1016/0045-7825(84)90072-0>.

    Abstract:
    Quadrature formulae of degrees 2 to 6 are presented for the numerical
    integration of a function over tetrahedronal regions. The formulae
    presented are of Gaussian type and fully symmetric with respect to the four
    vertices of the tetrahedron.
    '''
    def __init__(self, index):
        if index == 1:
            self.weights = numpy.full(4, 0.25)
            bary = _s31(0.138196601125015)
            self.degree = 2
        elif index == 2:
            self.weights = numpy.concatenate([
                numpy.full(1, -0.8),
                numpy.full(4, 0.45)
                ])
            bary = numpy.concatenate([
                _s4(),
                _s31(1.0/6.0)
                ])
            self.degree = 3
        elif index == 3:
            self.weights = numpy.concatenate([
                numpy.full(4, 0.5037379410012282E-01),
                numpy.full(12, 0.6654206863329239E-01)
                ])
            bary = numpy.concatenate([
                _s31(0.7611903264425430E-01),
                _s211(0.4042339134672644, 0.1197005277978019)
                ])
            self.degree = 4
        elif index == 4:
            self.weights = numpy.concatenate([
                numpy.full(1, 0.1884185567365411),
                numpy.full(4, 0.6703858372604275E-01),
                numpy.full(12, 0.4528559236327399E-01)
                ])
            bary = numpy.concatenate([
                _s4(),
                _s31(0.8945436401412733E-01),
                _s211(0.4214394310662522, 0.1325810999384657),
                ])
            self.degree = 5
        else:
            assert index == 5
            self.weights = numpy.concatenate([
                numpy.full(1, 0.9040129046014750E-01),
                numpy.full(4, 0.1911983427899124E-01),
                numpy.full(12, 0.4361493840666568E-01),
                numpy.full(12, 0.2581167596199161E-01)
                ])
            bary = numpy.concatenate([
                _s4(),
                _s31(0.5742691731735682E-01),
                _s211(0.2312985436519147, 0.5135188412556341E-01),
                _s211(0.4756909881472290E-01, 0.2967538129690260),
                ])
            self.degree = 6

        self.points = bary[:, [1, 2, 3]]
        return

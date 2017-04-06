# -*- coding: utf-8 -*-
#
from .helpers import _s4, _s31, _s22, _s211

import numpy


class ShunnHam(object):
    '''
    Lee Shunn, Frank Ham,
    Symmetric quadrature rules for tetrahedra based on a cubic
    close-packed lattice arrangement,
    Journal of Computational and Applied Mathematics,
    2012,
    <http://dx.doi.org/10.1016/j.cam.2012.03.032>.

    Abstract:
    A family of quadrature rules for integration over tetrahedral volumes is
    developed. The underlying structure of the rules is based on the cubic
    close-packed (CCP) lattice arrangement using 1, 4, 10, 20, 35, and 56
    quadrature points. The rules are characterized by rapid convergence,
    positive weights, and symmetry. Each rule is an optimal approximation in
    the sense that lower-order terms have zero contribution to the truncation
    error and the leading-order error term is minimized. Quadrature formulas up
    to order 9 are presented with relevant numerical examples.
    '''
    def __init__(self, index):
        if index == 1:
            self.weights = numpy.array([
                1.0
                ])
            bary = _s4()
            self.degree = 1
        elif index == 2:
            self.weights = 0.25 * numpy.ones(4)
            bary = _s31(0.1381966011250110)
            self.degree = 2
        elif index == 3:
            self.weights = numpy.concatenate([
                0.0476331348432089 * numpy.ones(4),
                0.1349112434378610 * numpy.ones(6),
                ])
            bary = numpy.concatenate([
                _s31(0.0738349017262234),
                _s22(0.0937556561159491),
                ])
            self.degree = 3
        elif index == 4:
            self.weights = numpy.concatenate([
                0.0070670747944695 * numpy.ones(4),
                0.0469986689718877 * numpy.ones(12),
                0.1019369182898680 * numpy.ones(4),
                ])
            bary = numpy.concatenate([
                _s31(0.0323525947272439),
                _s211(0.0603604415251421, 0.2626825838877790),
                _s31(0.3097693042728620),
                ])
            self.degree = 5
        elif index == 5:
            self.weights = numpy.concatenate([
                0.0021900463965388 * numpy.ones(4),
                0.0143395670177665 * numpy.ones(12),
                0.0250305395686746 * numpy.ones(6),
                0.0479839333057554 * numpy.ones(12),
                0.0931745731195340 * numpy.ones(1)
                ])
            bary = numpy.concatenate([
                _s31(0.0267367755543735),
                _s211(0.0391022406356488, 0.7477598884818090),
                _s22(0.0452454000155172),
                _s211(0.2232010379623150, 0.0504792790607720),
                numpy.array([[0.25, 0.25, 0.25, 0.25]]),
                ])
            self.degree = 6
        else:
            assert index == 6
            self.weights = numpy.concatenate([
                0.0010373112336140 * numpy.ones(4),
                0.0096016645399480 * numpy.ones(12),
                0.0164493976798232 * numpy.ones(12),
                0.0153747766513310 * numpy.ones(12),
                0.0293520118375230 * numpy.ones(12),
                0.0366291366405108 * numpy.ones(4),
                ])
            bary = numpy.concatenate([
                _s31(0.0149520651530592),
                _s211(0.0340960211962615, 0.1518319491659370),
                _s211(0.0462051504150017, 0.5526556431060170),
                _s211(0.2281904610687610, 0.0055147549744775),
                _s211(0.3523052600879940, 0.0992057202494530),
                _s31(0.1344783347929940)
                ])
            self.degree = 8

        self.points = bary[:, 1:]
        return

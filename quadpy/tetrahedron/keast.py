# -*- coding: utf-8 -*-
#
import numpy

from .helpers import _s4, _s31, _s22, _s211


class Keast(object):
    '''
    P. Keast,
    Moderate degree tetrahedral quadrature formulas,
    CMAME 55: 339-348
    1986,
    <http://dx.doi.org/10.1016/0045-7825(86)90059-9>.

    Abstract:
    Quadrature formulas of degrees 4 to 8 for numerical integration over the
    tetrahedron are constructed. The formulas are fully symmetric with respect
    to the tetrahedron, and in some cases are the minimum point rules with this
    symmetry.

    https://people.sc.fsu.edu/~jburkardt/datasets/quadrature_rules_tet/quadrature_rules_tet.html
    '''
    def __init__(self, index):
        if index == 0:
            # Does no appear in Keast's article.
            self.weights = numpy.array([
                1.0
                ])
            bary = _s4()
            self.degree = 1
        elif index == 1:
            # Does no appear in Keast's article.
            self.weights = numpy.full(4, 0.25)
            bary = _s31(0.1381966011250105)
            self.degree = 2
        elif index == 2:
            # Does no appear in Keast's article.
            self.weights = numpy.concatenate([
                numpy.full(1, -0.8),
                numpy.full(4, 0.45),
                ])
            bary = numpy.concatenate([
                _s4(),
                _s31(1.0/6.0),
                ])
            self.degree = 3
        elif index == 3:
            # Does no appear in Keast's article.
            self.weights = numpy.concatenate([
                numpy.full(4, 0.2177650698804054),
                numpy.full(6, 0.0214899534130631),
                ])
            bary = numpy.concatenate([
                _s31(0.1438564719343852),
                _s22(0.5),
                ])
            self.degree = 3
        elif index == 4:
            self.weights = numpy.concatenate([
                -148.0 / numpy.full(1, 1875.0),
                343.0 / numpy.full(4, 7500.0),
                56.0 / numpy.full(6, 375.0),
                ])
            bary = numpy.concatenate([
                _s4(),
                _s31(1.0/14.0),
                _s22(0.3994035761667992),
                ])
            self.degree = 4
        elif index == 5:
            self.weights = numpy.concatenate([
                numpy.full(6, 2.0/105.0),
                numpy.full(4, 0.0885898247429807),
                numpy.full(4, 0.1328387466855907),
                ])
            bary = numpy.concatenate([
                _s22(0.5),
                _s31(0.1005267652252045),
                _s31(0.3143728734931922),
                ])
            self.degree = 4
        elif index == 6:
            self.weights = numpy.concatenate([
                6544.0 / numpy.full(1, 36015.0),
                81.0 / numpy.full(4, 2240.0),
                161051.0 / numpy.full(4, 2304960.0),
                338.0 / numpy.full(6, 5145.0),
                ])
            bary = numpy.concatenate([
                _s4(),
                _s31(1.0/3.0),
                _s31(1.0/11.0),
                _s22(0.0665501535736643),
                ])
            self.degree = 5
        elif index == 7:
            self.weights = numpy.concatenate([
                numpy.full(4, 0.0399227502581679),
                numpy.full(4, 0.0100772110553207),
                numpy.full(4, 0.0553571815436544),
                numpy.full(12, 27.0/560.0),
                ])
            bary = numpy.concatenate([
                _s31(0.2146028712591517),
                _s31(0.0406739585346113),
                _s31(0.3223378901422757),
                _s211(0.0636610018750175, 0.2696723314583159)
                ])
            self.degree = 6
        elif index == 8:
            self.weights = numpy.concatenate([
                numpy.full(1, 0.1095853407966528),
                numpy.full(4, 0.0635996491464850),
                numpy.full(4, -0.3751064406859797),
                numpy.full(4, 0.0293485515784412),
                numpy.full(6, 0.0058201058201058),
                numpy.full(12, 0.1653439153439105),
                ])
            bary = numpy.concatenate([
                _s4(),
                _s31(0.0782131923303186),
                _s31(0.1218432166639044),
                _s31(0.3325391644464206),
                _s22(0.5),
                _s211(0.1, 0.2),
                ])
            self.degree = 7
        elif index == 9:
            self.weights = numpy.concatenate([
                numpy.full(1, -0.2359620398477557),
                numpy.full(4, 0.0244878963560562),
                numpy.full(4, 0.0039485206398261),
                numpy.full(6, 0.0263055529507371),
                numpy.full(6, 0.0829803830550589),
                numpy.full(12, 0.0254426245481023),
                numpy.full(12, 0.0134324384376852),
                ])
            bary = numpy.concatenate([
                _s4(),
                _s31(0.1274709365666390),
                _s31(0.0320788303926323),
                _s22(0.0497770956432810),
                _s22(0.1837304473985499),
                _s211(0.2319010893971509, 0.5132800333608811),
                _s211(0.0379700484718286, 0.1937464752488044),
                ])
            self.degree = 8
        else:
            assert index == 10
            self.weights = 6 * numpy.concatenate([
                # Note: In Keast's article, the first weight is incorrectly
                # given with a positive sign.
                numpy.full(1, -0.393270066412926145e-01),
                numpy.full(4, +0.408131605934270525e-02),
                numpy.full(4, +0.658086773304341943e-03),
                numpy.full(6, +0.438425882512284693e-02),
                numpy.full(6, +0.138300638425098166e-01),
                numpy.full(12, +0.424043742468372453e-02),
                numpy.full(12, +0.223873973961420164e-02),
                ])
            bary = numpy.concatenate([
                _s4(),
                _s31(0.127470936566639015e-00),
                _s31(0.320788303926322960e-01),
                _s22(0.497770956432810185e-01),
                _s22(0.183730447398549945e-00),
                _s211(0.231901089397150906e-00, 0.229177878448171174e-01),
                _s211(0.379700484718286102e-01, 0.730313427807538396e-00),
                ])
            self.degree = 8

        self.points = bary[:, 1:]
        return

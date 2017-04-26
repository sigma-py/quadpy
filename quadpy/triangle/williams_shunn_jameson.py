# -*- coding: utf-8 -*-
#
from .helpers import _s3, _s21, _s111
import numpy


class WilliamsShunnJameson(object):
    '''
    D.M. Williams, L. Shunn, A. Jameson,
    Symmetric quadrature rules for simplexes based on sphere close packed
    lattice arrangements,
    Journal of Computational and Applied Mathematics,
    266 (2014) 18–38,
    <https://doi.org/10.1016/j.cam.2014.01.007>.

    Abstract:
    Sphere close packed (SCP) lattice arrangements of points are well-suited
    for formulating symmetric quadrature rules on simplexes, as they are
    symmetric under affine transformations of the simplex unto itself in 2D and
    3D. As a result, SCP lattice arrangements have been utilized to formulate
    symmetric quadrature rules with Np = 1, 4, 10, 20, 35, and 56 points on the
    3-simplex (Shunn and Ham, 2012). In what follows, the work on the 3-simplex
    is extended, and SCP lattices are employed to identify symmetric quadrature
    rules with Np = 1, 3, 6, 10, 15, 21, 28, 36, 45, 55, and 66 points on the
    2-simplex and Np = 84 points on the 3-simplex. These rules are found to be
    capable of exactly integrating polynomials of up to degree 17 in 2D and up
    to degree 9 in 3D.
    '''
    def __init__(self, index):
        self.name = 'WSJ(%d)' % index
        if index == 1:
            self.weights = [1.0]
            bary = _s3()
            self.degree = 1
        elif index == 2:
            self.weights = 3 * [1.0 / 3.0]
            bary = _s21(1.0 / 6.0)
            self.degree = 2
        elif index == 3:
            self.weights = (
                3 * [0.109951743655333] +
                3 * [0.223381589678000]
                )
            bary = numpy.concatenate([
                _s21(0.091576213509780),
                _s21(0.445948490915964)
                ])
            self.degree = 4
        elif index == 4:
            self.weights = (
                3 * [0.041955512996649] +
                6 * [0.112098412070887] +
                [0.201542988584730]
                )
            bary = numpy.concatenate([
                _s21(0.055564052669793),
                _s111(0.295533711735893, 0.634210747745723),
                _s3(),
                ])
            self.degree = 5
        elif index == 5:
            self.weights = (
                3 * [0.017915455012303] +
                3 * [0.127712195881265] +
                3 * [0.076206062385535] +
                6 * [0.055749810027115]
                )
            bary = numpy.concatenate([
                _s21(0.035870877695734),
                _s21(0.241729395767967),
                _s21(0.474308787777079),
                _s111(0.201503881881800, 0.751183631106484),
                ])
            self.degree = 7
        elif index == 6:
            self.weights = (
                3 * [0.010359374696538] +
                3 * [0.075394884326738] +
                3 * [0.097547802373242] +
                6 * [0.028969269372473] +
                6 * [0.046046366595935]
                )
            bary = numpy.concatenate([
                _s21(0.028112952182664),
                _s21(0.177139098469317),
                _s21(0.405508595867433),
                _s111(0.148565812270887, 0.817900980028499),
                _s111(0.357196298615681, 0.604978911775132),
                ])
            self.degree = 8
        elif index == 7:
            self.weights = (
                3 * [0.005272170280495] +
                3 * [0.044552936679504] +
                1 * [0.083608212215637] +
                3 * [0.033815712804198] +
                6 * [0.015710461340183] +
                6 * [0.028205136280616] +
                6 * [0.066995957127830]
                )
            bary = numpy.concatenate([
                _s21(0.019977187122193),
                _s21(0.131721767529998),
                _s3(),
                _s21(0.485135346793461),
                _s111(0.107951981846011, 0.867911210117951),
                _s111(0.270840772921567, 0.700872570380723),
                _s111(0.316549598844617, 0.536654684206138),
                ])
            self.degree = 10
        else:
            assert index == 8
            self.weights = (
                3 * [0.005639123786910] +
                3 * [0.027148968192278] +
                3 * [0.063100912533359] +
                3 * [0.051752795679899] +
                6 * [0.009866753574646] +
                6 * [0.022008204800147] +
                6 * [0.016644570076736] +
                6 * [0.044326238118914]
                )
            bary = numpy.concatenate([
                _s21(0.021171422779465),
                _s21(0.100584397395888),
                _s21(0.271038307711932),
                _s21(0.440191258403832),
                _s111(0.101763679498021, 0.879979641427232),
                _s111(0.394033271669987, 0.582562022863673),
                _s111(0.226245530909229, 0.751530614542782),
                _s111(0.635737183263105, 0.249079227621332),
                ])
            self.degree = 12

        self.weights = numpy.array(self.weights)
        bary = numpy.array(bary)
        self.points = bary[:, [1, 2]]
        return

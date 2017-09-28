# -*- coding: utf-8 -*-
#
from sympy import Rational as fr

from .helpers import _s3, _s21, _s111 as fs
from ..helpers import untangle


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
            self.degree = 2
            data = [(fr(1, 3), _s21(fr(1, 6)))]
        elif index == 2:
            self.degree = 2
            data = [(fr(1, 3), _s21(fr(1, 2)))]
        elif index == 3:
            self.degree = 3
            data = [
                (-fr(9, 16), _s3()),
                (fr(25, 48), _s21(fr(1, 5))),
                ]
        elif index == 4:
            self.degree = 3
            data = [
                (fr(1, 6), fs(0.659027622374092, 0.231933368553031))
                ]
        elif index == 5:
            self.degree = 4
            data = [
                (0.109951743655322, _s21(0.091576213509771)),
                (0.223381589678011, _s21(0.445948490915965)),
                ]
        elif index == 6:
            self.degree = 4
            data = [
                (fr(3, 8), _s3()),
                (fr(5, 48), fs(0.736712498968435, 0.237932366472434)),
                ]
        elif index == 7:
            self.degree = 5
            data = [
                (fr(9, 40), _s3()),
                (0.12593918054482717, _s21(0.10128650732345633)),
                (0.13239415278850616, _s21(0.47014206410511505)),
                ]
        elif index == 8:
            self.degree = 5
            data = [
                (0.205950504760887, _s21(0.437525248383384)),
                (0.063691414286223, fs(0.797112651860071, 0.165409927389841)),
                ]
        elif index == 9:
            self.degree = 6
            data = [
                (0.050844906370207, _s21(0.063089014491502)),
                (0.116786275726379, _s21(0.249286745170910)),
                (0.082851075618374, fs(0.636502499121399, 0.310352451033785)),
                ]
        else:
            assert index == 10
            self.degree = 7
            data = [
                (-0.149570044467670, _s3()),
                (+0.175615257433204, _s21(0.260345966079038)),
                (+0.053347235608839, _s21(0.065130102902216)),
                (+0.077113760890257, fs(0.638444188569809, 0.312865496004875)),
                ]

        self.bary, self.weights = untangle(data)
        self.points = self.bary[:, 1:]
        return

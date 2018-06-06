# -*- coding: utf-8 -*-
#
from .helpers import _s4, _s22, _s211, _s31

from ..helpers import untangle


class WitherdenVincent(object):
    """
    F.D. Witherden, P.E. Vincent,
    On the identification of symmetric quadrature rules for finite
    element methods,
    Computers & Mathematics with Applications,
    Volume 69, Issue 10, May 2015, Pages 1232–1241,
    <https://doi.org/10.1016/j.camwa.2015.03.017>.

    Abstract:
    In this paper we describe a methodology for the identification of symmetric
    quadrature rules inside of quadrilaterals, triangles, tetrahedra, prisms,
    pyramids, and hexahedra. The methodology is free from manual intervention
    and is capable of identifying a set of rules with a given strength and a
    given number of points. We also present polyquad which is an implementation
    of our methodology. Using polyquad v1.0 we proceed to derive a complete set
    of symmetric rules on the aforementioned domains. All rules possess purely
    positive weights and have all points inside the domain. Many of the rules
    appear to be new, and an improvement over those tabulated in the
    literature.
    """

    def __init__(self, degree):
        self.name = "WV({})".format(degree)
        self.degree = degree

        if degree == 1:
            data = [(1.000000000000000e+00, _s4())]
        elif degree == 2:
            data = [(2.500000000000000e-01, _s31(1.381966011250105e-01))]
        elif degree == 3:
            data = [
                (1.362178425370874e-01, _s31(3.281633025163817e-01)),
                (1.137821574629126e-01, _s31(1.080472498984286e-01)),
            ]
        elif degree == 5:
            data = [
                (1.126879257180159e-01, _s31(3.108859192633006e-01)),
                (7.349304311636196e-02, _s31(9.273525031089125e-02)),
                (4.254602077708147e-02, _s22(4.550370412564964e-02)),
            ]
        elif degree == 6:
            data = [
                (1.007721105532064e-02, _s31(4.067395853461137e-02)),
                (5.535718154365472e-02, _s31(3.223378901422755e-01)),
                (3.992275025816749e-02, _s31(2.146028712591520e-01)),
                (
                    4.821428571428571e-02,
                    _s211(6.366100187501744e-02, 6.030056647916492e-01),
                ),
            ]
        elif degree == 7:
            data = [
                (9.548528946413085e-02, _s4()),
                (4.232958120996703e-02, _s31(3.157011497782028e-01)),
                (3.189692783285758e-02, _s22(5.048982259839635e-02)),
                (
                    3.720713072833462e-02,
                    _s211(1.888338310260010e-01, 5.751716375870000e-01),
                ),
                (
                    8.110770829903342e-03,
                    _s211(2.126547254148314e-02, 8.108302410985486e-01),
                ),
            ]
        elif degree == 8:
            data = [
                (2.642665090840883e-02, _s31(1.079527249622109e-01)),
                (5.203174756373853e-02, _s31(1.851094877825866e-01)),
                (7.525256153540199e-03, _s31(4.231654368476728e-02)),
                (4.176378285693490e-02, _s31(3.141817091240390e-01)),
                (3.628093026130882e-02, _s22(4.355913285838302e-01)),
                (
                    7.156902890844433e-03,
                    _s211(2.143393012713057e-02, 7.174640634263083e-01),
                ),
                (
                    1.545348615096034e-02,
                    _s211(2.041393338760291e-01, 5.837973783021444e-01),
                ),
            ]
        elif degree == 9:
            data = [
                (5.801054891248025e-02, _s4()),
                (6.431928175925639e-05, _s31(6.198169755222693e-10)),
                (2.317333846242546e-02, _s31(1.607745353952616e-01)),
                (2.956291233542929e-02, _s31(3.222765218214210e-01)),
                (8.063979979616182e-03, _s31(4.510891834541358e-02)),
                (3.813408010370246e-02, _s22(1.122965460043761e-01)),
                (
                    8.384422198298552e-03,
                    _s211(4.588714487524592e-01, 2.554579233041310e-03),
                ),
                (
                    1.023455935274533e-02,
                    _s211(3.377587068533860e-02, 7.183503264420745e-01),
                ),
                (
                    2.052491596798814e-02,
                    _s211(1.836413698099279e-01, 3.441591057817528e-02),
                ),
            ]
        else:
            assert degree == 10
            data = [
                (4.739977355602074e-02, _s4()),
                (2.693705999226870e-02, _s31(3.122500686951887e-01)),
                (9.869159716793382e-03, _s31(1.143096538573461e-01)),
                (
                    1.139388122019523e-02,
                    _s211(4.104307392189654e-01, 1.654860256196111e-01),
                ),
                (
                    3.619443443392536e-04,
                    _s211(6.138008824790653e-03, 9.429887673452049e-01),
                ),
                (
                    2.573973198045607e-02,
                    _s211(1.210501811455894e-01, 4.771903799042804e-01),
                ),
                (
                    1.013587167975579e-02,
                    _s211(3.277946821644262e-02, 5.942562694800070e-01),
                ),
                (
                    6.576147277035904e-03,
                    _s211(3.248528156482305e-02, 8.011772846583444e-01),
                ),
                (
                    1.290703579886199e-02,
                    _s211(1.749793421839390e-01, 6.280718454753660e-01),
                ),
            ]

        self.bary, self.weights = untangle(data)
        self.points = self.bary[:, 1:]
        return

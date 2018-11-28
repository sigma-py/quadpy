# -*- coding: utf-8 -*-
#
from .helpers import untangle2


class Triex(object):
    """
    See
    https://people.sc.fsu.edu/~jburkardt/datasets/quadrature_rules_tri/quadrature_rules_tri.html
    and

    E. de Doncker and I. Robinson,
    Algorithm 612: TRIEX: Integration Over a TRIangle Using Nonlinear
    EXtrapolation,
    ACM Trans. Math. Softw.,
    March 1984,
    <http://dl.acm.org/citation.cfm?id=356070&CFID=836775288&CFTOKEN=89206835>.
    """

    def __init__(self, index):
        self.name = "TRIEX({})".format(index)
        if index == 19:
            self.degree = 9
            data = {
                "s3": [[9.71357962827961025e-002]],
                "s2": [
                    [3.13347002271398278e-002, 0.48968251919873701],
                    [7.78275410047754301e-002, 0.43708959149293553],
                    [7.96477389272090969e-002, 0.18820353561903219],
                    [2.55776756586981006e-002, 4.47295133944529688e-002],
                ],
                "s1": [
                    [
                        4.32835393772893970e-002,
                        0.74119859878449801,
                        3.68384120547362581e-002,
                    ]
                ],
            }
        else:
            assert index == 28
            self.degree = 11
            data = {
                "s3": [[8.797730116222190e-2]],
                "s2": [
                    [8.744311553736190e-3, 0.2598914092828833e-01],
                    [3.808157199393533e-2, 0.9428750264792270e-01],
                    [1.885544805613125e-2, 0.4946367750172147],
                    [7.215969754474100e-2, 0.2073433826145142],
                    [6.932913870553720e-2, 0.4389078057004907],
                ],
                "s1": [
                    [4.105631542928860e-2, 0.6779376548825902, 0.04484167758913055],
                    [7.362383783300573e-3, 0.8588702812826364, 0],
                ],
            }

        self.bary, self.weights = untangle2(data)
        self.points = self.bary[:, 1:]
        return

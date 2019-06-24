# -*- coding: utf-8 -*-
#
from ._helpers import concat, s1, s2, s3, TriangleScheme
from ..helpers import article

citation = article(
    authors=["E. de Doncker", "I. Robinson"],
    title="Algorithm 612: TRIEX: Integration Over a TRIangle Using Nonlinear EXtrapolation",
    journal="ACM Trans. Math. Softw.",
    month="mar",
    year="1984",
    url="https://doi.org/10.1145/356068.356070",
)
# See also
# https://people.sc.fsu.edu/~jburkardt/datasets/quadrature_rules_tri/quadrature_rules_tri.html


def triex_19():
    weights, bary = concat(
        s3(9.71357962827961025e-002),
        s2(
            [3.13347002271398278e-002, 0.48968251919873701],
            [7.78275410047754301e-002, 0.43708959149293553],
            [7.96477389272090969e-002, 0.18820353561903219],
            [2.55776756586981006e-002, 4.47295133944529688e-002],
        ),
        s1([4.32835393772893970e-002, 0.74119859878449801, 3.68384120547362581e-002]),
    )
    return TriangleScheme("TRIEX 19", weights, bary, 9, citation)


def triex_28():
    # DUP Lyness-Jespersen 21
    weights, bary = concat(
        s3(8.797730116222190e-2),
        s2(
            [8.744311553736190e-3, 0.2598914092828833e-01],
            [3.808157199393533e-2, 0.9428750264792270e-01],
            [1.885544805613125e-2, 0.4946367750172147],
            [7.215969754474100e-2, 0.2073433826145142],
            [6.932913870553720e-2, 0.4389078057004907],
        ),
        s1(
            [4.105631542928860e-2, 0.6779376548825902, 0.04484167758913055],
            [7.362383783300573e-3, 0.8588702812826364, 0],
        ),
    )
    return TriangleScheme("TRIEX 28", weights, bary, 11, citation)

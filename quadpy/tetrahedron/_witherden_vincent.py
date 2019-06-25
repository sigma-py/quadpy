# -*- coding: utf-8 -*-
#
from ..helpers import article
from ._helpers import TetrahedronScheme, untangle2

citation = article(
    authors=["F.D. Witherden", "P.E. Vincent"],
    title="On the identification of symmetric quadrature rules for finite element methods",
    journal="Computers & Mathematics with Applications",
    volume="69",
    number="10",
    month="may",
    year="2015",
    pages="1232–1241",
    url="https://doi.org/10.1016/j.camwa.2015.03.017",
)


def witherden_vincent_01():
    degree = 1
    data = {"s4": [[1.000000000000000e00]]}
    points, weights = untangle2(data)
    return TetrahedronScheme("Witherden-Vincent 1", weights, points, degree, citation)


def witherden_vincent_02():
    degree = 2
    data = {"s31": [[2.500000000000000e-01, 1.381966011250105e-01]]}
    points, weights = untangle2(data)
    return TetrahedronScheme("Witherden-Vincent 2", weights, points, degree, citation)


def witherden_vincent_03():
    degree = 3
    data = {
        "s31": [
            [1.362178425370874e-01, 3.281633025163817e-01],
            [1.137821574629126e-01, 1.080472498984286e-01],
        ]
    }
    points, weights = untangle2(data)
    return TetrahedronScheme("Witherden-Vincent 3", weights, points, degree, citation)


def witherden_vincent_05():
    degree = 5
    data = {
        "s31": [
            [1.126879257180159e-01, 3.108859192633006e-01],
            [7.349304311636196e-02, 9.273525031089125e-02],
        ],
        "s22": [[4.254602077708147e-02, 4.550370412564964e-02]],
    }
    points, weights = untangle2(data)
    return TetrahedronScheme("Witherden-Vincent 5", weights, points, degree, citation)


def witherden_vincent_06():
    degree = 6
    data = {
        "s31": [
            [1.007721105532064e-02, 4.067395853461137e-02],
            [5.535718154365472e-02, 3.223378901422755e-01],
            [3.992275025816749e-02, 2.146028712591520e-01],
        ],
        "s211": [[4.821428571428571e-02, 6.366100187501744e-02, 6.030056647916492e-01]],
    }
    points, weights = untangle2(data)
    return TetrahedronScheme("Witherden-Vincent 6", weights, points, degree, citation)


def witherden_vincent_07():
    degree = 7
    data = {
        "s4": [[9.548528946413085e-02]],
        "s31": [[4.232958120996703e-02, 3.157011497782028e-01]],
        "s22": [[3.189692783285758e-02, 5.048982259839635e-02]],
        "s211": [
            [3.720713072833462e-02, 1.888338310260010e-01, 5.751716375870000e-01],
            [8.110770829903342e-03, 2.126547254148314e-02, 8.108302410985486e-01],
        ],
    }
    points, weights = untangle2(data)
    return TetrahedronScheme("Witherden-Vincent 7", weights, points, degree, citation)


def witherden_vincent_08():
    degree = 8
    data = {
        "s31": [
            [2.642665090840883e-02, 1.079527249622109e-01],
            [5.203174756373853e-02, 1.851094877825866e-01],
            [7.525256153540199e-03, 4.231654368476728e-02],
            [4.176378285693490e-02, 3.141817091240390e-01],
        ],
        "s22": [[3.628093026130882e-02, 4.355913285838302e-01]],
        "s211": [
            [7.156902890844433e-03, 2.143393012713057e-02, 7.174640634263083e-01],
            [1.545348615096034e-02, 2.041393338760291e-01, 5.837973783021444e-01],
        ],
    }
    points, weights = untangle2(data)
    return TetrahedronScheme("Witherden-Vincent 8", weights, points, degree, citation)


def witherden_vincent_09():
    degree = 9
    data = {
        "s4": [[5.801054891248025e-02]],
        "s31": [
            [6.431928175925639e-05, 6.198169755222693e-10],
            [2.317333846242546e-02, 1.607745353952616e-01],
            [2.956291233542929e-02, 3.222765218214210e-01],
            [8.063979979616182e-03, 4.510891834541358e-02],
        ],
        "s22": [[3.813408010370246e-02, 1.122965460043761e-01]],
        "s211": [
            [8.384422198298552e-03, 4.588714487524592e-01, 2.554579233041310e-03],
            [1.023455935274533e-02, 3.377587068533860e-02, 7.183503264420745e-01],
            [2.052491596798814e-02, 1.836413698099279e-01, 3.441591057817528e-02],
        ],
    }
    points, weights = untangle2(data)
    return TetrahedronScheme("Witherden-Vincent 9", weights, points, degree, citation)


def witherden_vincent_10():
    degree = 10
    data = {
        "s4": [[4.739977355602074e-02]],
        "s31": [
            [2.693705999226870e-02, 3.122500686951887e-01],
            [9.869159716793382e-03, 1.143096538573461e-01],
        ],
        "s211": [
            [1.139388122019523e-02, 4.104307392189654e-01, 1.654860256196111e-01],
            [3.619443443392536e-04, 6.138008824790653e-03, 9.429887673452049e-01],
            [2.573973198045607e-02, 1.210501811455894e-01, 4.771903799042804e-01],
            [1.013587167975579e-02, 3.277946821644262e-02, 5.942562694800070e-01],
            [6.576147277035904e-03, 3.248528156482305e-02, 8.011772846583444e-01],
            [1.290703579886199e-02, 1.749793421839390e-01, 6.280718454753660e-01],
        ],
    }
    points, weights = untangle2(data)
    return TetrahedronScheme("Witherden-Vincent 10", weights, points, degree, citation)

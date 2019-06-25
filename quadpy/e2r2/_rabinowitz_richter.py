# -*- coding: utf-8 -*-
#
from ..helpers import article, untangle
from ._helpers import E2r2Scheme, _s4, _s8, _s40, _z

_citation = article(
    authors=["Philip Rabinowitz", "Nira Richter"],
    title="Perfectly Symmetric Two-Dimensional Integration Formulas with Minimal Numbers of Points",
    journal="Mathematics of Computation",
    volume="23",
    number="108",
    month="oct",
    year="1969",
    pages="765-779",
    url="https://doi.org/10.2307/2004962",
)


def rabinowitz_richter_1():
    data = [
        (0.1237222328857347e00, _s40(1.538189001320852)),
        (0.6544984694978697e-01, _s4(1.224744871391589)),
        (0.5935280476180875e00, _s4(0.4817165220011443)),
        (0.1349017971918148e-02, _s8(2.607349811958554, 0.9663217712794149)),
    ]
    points, weights = untangle(data)
    return E2r2Scheme("RabinowitzRichter 1", weights, points, 9, _citation)


def rabinowitz_richter_2():
    data = [
        (0.8176645817675417e-3, _s40(2.757816396257008)),
        (0.4363323129985824e-1, _s40(1.732050807568877)),
        (0.5373255214498174e0, _s40(0.6280515301597559)),
        (0.3636102608321520e-2, _s8(1.224744871391589, 2.121320343559643)),
        (0.9817477042468103e-1, _s8(0.7071067811865475, 1.224744871391589)),
    ]
    points, weights = untangle(data)
    return E2r2Scheme("RabinowitzRichter 2", weights, points, 11, _citation)


def rabinowitz_richter_3():
    data = [
        (0.4106569066965604e-3, _s40(2.907364117106118)),
        (0.9065690889492120e-1, _s40(1.528230917660483)),
        (0.5266955729327722e0, _s40(0.6178819071436261)),
        (0.9681125175723808e-3, _s4(1.904162039910276)),
        (0.1515812331366514e0, _s4(0.9724173472297303)),
        (0.7542839504417270e-2, _s8(2.061552812808830, 0.8660254037844387)),
    ]
    points, weights = untangle(data)
    return E2r2Scheme("RabinowitzRichter 3", weights, points, 11, _citation)


def rabinowitz_richter_4():
    data = [
        (-0.7482913219380363e0, _z()),
        (+0.3521509661098668e-2, _s40(2.403151765001966)),
        (+0.1650055872539264e0, _s40(1.298479973315986)),
        (+0.8537825937946404e-3, _s4(1.912428205769905)),
        (+0.1326938806789336e0, _s4(0.9478854439698223)),
        (+0.6447719928481539e0, _s4(0.3188824732576547)),
        (+0.1799266413507747e-4, _s8(3.325657829663178, 1.145527285699371)),
        (+0.1279412775888998e-1, _s8(1.882228401823884, 0.8826073082889659)),
    ]
    points, weights = untangle(data)
    return E2r2Scheme("RabinowitzRichter 4", weights, points, 13, _citation)


def rabinowitz_richter_5():
    data = [
        (0.8006483569659628e-5, _s40(3.538388728121807)),
        (0.3604577420838264e-2, _s40(2.359676416877929)),
        (0.1187609330759137e0, _s40(1.312801844620926)),
        (0.4372488543791402e0, _s40(0.5389559482114205)),
        (0.3671735075832989e-4, _s4(2.300279949805658)),
        (0.5654866776461627e-2, _s4(1.581138830084189)),
        (0.1777774268424240e0, _s4(0.8418504335819279)),
        (0.2735449647853290e-3, _s8(2.685533581755341, 1.112384431771456)),
        (0.2087984556938594e-1, _s8(1.740847514397403, 0.7210826504868960)),
    ]
    points, weights = untangle(data)
    return E2r2Scheme("RabinowitzRichter 5", weights, points, 15, _citation)

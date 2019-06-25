# -*- coding: utf-8 -*-
#

import sympy

from ..helpers import article
from ._helpers import TriangleScheme, concat, s1, s2, s3

citation = article(
    authors=["M.E. Laursen", "M. Gellert"],
    title="Some criteria for numerically integrated matrices and quadrature formulas for triangles",
    journal="International Journal for Numerical Methods in Engineering",
    volume="12",
    number="1",
    year="1978",
    pages="67–76",
    url="https://doi.org/10.1002/nme.1620120107",
)


def laursen_gellert_01(symbolic=False):
    weights, points = s3(1)
    return TriangleScheme("Laursen-Gellert 1", weights, points, 1, citation)


def laursen_gellert_02a(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    weights, points = s2([frac(1, 3), frac(1, 6)])
    return TriangleScheme("Laursen-Gellert 2a", weights, points, 2, citation)


def laursen_gellert_02b(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    weights, points = s2([frac(1, 3), frac(1, 2)])
    return TriangleScheme("Laursen-Gellert 2b", weights, points, 2, citation)


def laursen_gellert_03(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    weights, points = concat(s3(-frac(9, 16)), s2([frac(25, 48), frac(1, 5)]))
    return TriangleScheme("Laursen-Gellert 3", weights, points, 3, citation)


def laursen_gellert_04(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    weights, points = s1([frac(1, 6), 0.659027622374092, 0.231933368553031])
    return TriangleScheme("Laursen-Gellert 4", weights, points, 3, citation)


def laursen_gellert_05():
    weights, points = s2(
        [0.109951743655322, 0.091576213509771], [0.223381589678011, 0.445948490915965]
    )
    return TriangleScheme("Laursen-Gellert 5", weights, points, 4, citation)


def laursen_gellert_06(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    weights, points = concat(
        s3(frac(3, 8)), s1([frac(5, 48), 0.736712498968435, 0.237932366472434])
    )
    return TriangleScheme("Laursen-Gellert 6", weights, points, 4, citation)


def laursen_gellert_07(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    weights, points = concat(
        s3(frac(9, 40)),
        s2(
            [0.125939180544827, 0.101286507323456],
            [0.132394152788506, 0.470142064105115],
        ),
    )
    return TriangleScheme("Laursen-Gellert 7", weights, points, 5, citation)


def laursen_gellert_08():
    weights, points = concat(
        s2([0.205950504760887, 0.437525248383384]),
        s1([0.063691414286223, 0.797112651860071, 0.165409927389841]),
    )
    return TriangleScheme("Laursen-Gellert 8", weights, points, 5, citation)


def laursen_gellert_09():
    weights, points = concat(
        s2(
            [0.050844906370207, 0.063089014491502],
            [0.116786275726379, 0.249286745170910],
        ),
        s1([0.082851075618374, 0.636502499121399, 0.310352451033785]),
    )
    return TriangleScheme("Laursen-Gellert 9", weights, points, 6, citation)


def laursen_gellert_10():
    weights, points = concat(
        s3(-0.149570044467670),
        s2(
            [+0.175615257433204, 0.260345966079038],
            [+0.053347235608839, 0.065130102902216],
        ),
        s1([+0.077113760890257, 0.638444188569809, 0.312865496004875]),
    )
    return TriangleScheme("Laursen-Gellert 10", weights, points, 7, citation)


def laursen_gellert_11():
    weights, points = concat(
        s2([0.053077801790233, 0.064930513159165]),
        s1(
            [0.070853083692136, 0.284575584249173, 0.517039939069325],
            [0.069274682079415, 0.313559184384932, 0.043863471792371],
        ),
    )
    return TriangleScheme("Laursen-Gellert 11", weights, points, 7, citation)


def laursen_gellert_12():
    weights, points = concat(
        s3(0.144315607677787),
        s2(
            [0.103217370534718, 0.170569307751761],
            [0.032458497623198, 0.050547228317031],
            [0.095091634267284, 0.459292588292723],
        ),
        s1([0.027230314174435, 0.008394777409958, 0.263112829634638]),
    )
    return TriangleScheme("Laursen-Gellert 12", weights, points, 8, citation)


def laursen_gellert_13():
    weights, points = concat(
        s3(0.097135796282799),
        s2(
            [0.031334700227139, 0.489682519198738],
            [0.077827541004774, 0.437089591492937],
            [0.079647738927210, 0.188203535619033],
            [0.025577675658698, 0.044729513394453],
        ),
        s1([0.043283539377289, 0.036838412054736, 0.221962989160766]),
    )
    return TriangleScheme("Laursen-Gellert 13", weights, points, 9, citation)


def laursen_gellert_14():
    weights, points = concat(
        s2(
            [0.051617202569021, 0.481519834783311],
            [0.094080073458356, 0.403603979817940],
            [0.025993571032320, 0.045189009784377],
        ),
        s1(
            [0.045469538047619, 0.136991201264904, 0.218290070971381],
            [0.035351705089199, 0.030424361728820, 0.222063165537318],
        ),
    )
    return TriangleScheme("Laursen-Gellert 14", weights, points, 9, citation)


def laursen_gellert_15a():
    weights, points = concat(
        s3(0.079894504741240),
        s2(
            [0.071123802232377, 0.425086210602091],
            [0.008223818690464, 0.023308867510000],
        ),
        s1(
            [0.045430592296170, 0.147925626209534, 0.223766973576973],
            [0.037359856234305, 0.029946031954171, 0.358740141864431],
            [0.030886656884564, 0.035632559587504, 0.143295370426867],
        ),
    )
    return TriangleScheme("Laursen-Gellert 15a", weights, points, 10, citation)


def laursen_gellert_15b():
    weights, points = concat(
        s3(0.081743329146286),
        s2(
            [0.045957963604745, 0.142161101056564],
            [0.013352968813150, 0.032055373216944],
        ),
        s1(
            [0.063904906396424, 0.148132885783821, 0.321812995288835],
            [0.034184648162959, 0.029619889488730, 0.369146781827811],
            [0.025297757707288, 0.028367665339938, 0.163701733737182],
        ),
    )
    return TriangleScheme("Laursen-Gellert 15b", weights, points, 10, citation)

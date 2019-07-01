# -*- coding: utf-8 -*-
#
from sympy import Rational as frac
from sympy import sqrt

from ..helpers import article
from ._helpers import QuadrilateralScheme, concat, symm_r0, symm_s, symm_s_t, zero

citation = article(
    authors=["D.A. Dunavant"],
    title="Economical symmetrical quadrature rules for complete polynomials over a square domain",
    journal="Numerical Methods in Engineering",
    volume="21",
    number="10",
    month="oct",
    year="1985",
    pages="1777–1784",
    url="https://doi.org/10.1002/nme.1620211004",
)


def dunavant_00():
    weights, points = zero(4)
    return QuadrilateralScheme("Dunavant 0", weights, points, 1, citation)


def dunavant_01():
    weights, points = symm_s([1, sqrt(frac(1, 3))])
    return QuadrilateralScheme("Dunavant 1", weights, points, 3, citation)


def dunavant_02():
    weights, points = concat(
        symm_r0([frac(40, 49), sqrt(frac(7, 15))]),
        symm_s([frac(9, 49), sqrt(frac(7, 9))]),
    )
    return QuadrilateralScheme("Dunavant 2", weights, points, 5, citation)


def dunavant_03():
    weights, points = concat(
        symm_r0([frac(98, 405), sqrt(frac(6, 7))]),
        symm_s(
            [0.237431774690630, 0.805979782918599],
            [0.520592916667394, 0.380554433208316],
        ),
    )
    return QuadrilateralScheme("Dunavant 3", weights, points, 7, citation)


def dunavant_04():
    weights, points = concat(
        symm_r0(
            [0.018475842507491, 1.121225763866564],
            [0.390052939160735, 0.451773049920657],
        ),
        symm_s([0.083095178026482, 0.891849420851512]),
        symm_s_t([0.254188020152646, 0.824396370749276, 0.411623426336542]),
    )
    return QuadrilateralScheme("Dunavant 4", weights, points, 9, citation)


def dunavant_05():
    weights, points = concat(
        zero(0.365379525585903),
        symm_r0(
            [0.027756165564204, 1.044402915409813],
            [0.244272057751754, 0.769799068396649],
        ),
        symm_s(
            [0.034265103851229, 0.935787012440540],
            [0.308993036133713, 0.413491953449114],
        ),
        symm_s_t([0.146684377651312, 0.883025508525690, 0.575653595840465]),
    )
    return QuadrilateralScheme("Dunavant 5", weights, points, 11, citation)


def dunavant_06():
    weights, points = concat(
        symm_r0(
            [0.005656169693764, 1.086056158573971],
            [0.192443867470396, 0.658208197042585],
        ),
        symm_s(
            [0.005166832979773, 1.001300602991729],
            [0.200302559622138, 0.584636168775946],
            [0.228125175912536, 0.246795612720261],
        ),
        symm_s_t(
            [0.117496926974491, 0.900258815287201, 0.304720678579870],
            [0.066655770186205, 0.929866705560780, 0.745052720131169],
        ),
    )
    return QuadrilateralScheme("Dunavant 6", weights, points, 13, citation)


def dunavant_07():
    weights, points = concat(
        zero(-0.001768979827207),
        symm_r0(
            [+0.012816726617512, 1.027314357719367],
            [+0.119897873101347, 0.856766776147643],
            [+0.210885452208801, 0.327332998189723],
        ),
        symm_s(
            [+0.006392720128215, 0.967223740028505],
            [+0.104415680788580, 0.732168901749711],
        ),
        symm_s_t(
            [0.168053047203816, 0.621974427996805, 0.321696694921009],
            [0.076169694452294, 0.928618480068352, 0.455124178121179],
            [0.028794154400064, 0.960457474887516, 0.809863684081217],
        ),
    )
    return QuadrilateralScheme("Dunavant 7", weights, points, 15, citation)


def dunavant_08():
    weights, points = concat(
        symm_r0(
            [0.020614915919991, 0.989353074512600],
            [0.128025716179910, 0.376285207157973],
        ),
        symm_s(
            [0.005511739534032, 0.978848279262233],
            [0.039207712457142, 0.885794729164116],
            [0.076396945079863, 0.171756123838348],
        ),
        symm_s_t(
            [0.141513729949972, 0.590499273806002, 0.319505036634574],
            [0.083903279363798, 0.799079131916863, 0.597972451929457],
            [0.060394163649685, 0.803743962958745, 0.058344481776551],
            [0.057387752969213, 0.936506276127495, 0.347386316166203],
            [0.021922559481864, 0.981321179805452, 0.706000287798646],
        ),
    )
    return QuadrilateralScheme("Dunavant 8", weights, points, 17, citation)


def dunavant_09():
    weights, points = concat(
        symm_r0(
            [0.038205406871462, 0.943962831808239],
            [0.135368502976521, 0.536918434376013],
        ),
        symm_s(
            [0.005773708558664, 0.973981076394170],
            [0.067460759759473, 0.742995535327609],
            [0.140899115227892, 0.285010052188916],
            [0.047466627685662, 0.068354569272491],
        ),
        symm_s_t(
            [0.078619467342982, 0.802952004398543, 0.203345534163332],
            [0.094979169511394, 0.634244672807882, 0.426572172992877],
            [0.022331162356015, 0.978350706908227, 0.295830776620995],
            [0.055594877793785, 0.901672714410389, 0.541983037327871],
            [0.006049054506376, 1.007018449383116, 0.669414798783936],
            [0.024839207949609, 0.945161453573471, 0.829501421477824],
        ),
    )
    # TODO ERR the article claims 19
    return QuadrilateralScheme("Dunavant 9", weights, points, 16, citation)


def dunavant_10():
    weights, points = concat(
        symm_r0(
            [0.019503841092684, 0.980883148832881],
            [0.089012127744268, 0.678152700336576],
            [0.114568584702749, 0.240599282275864],
        ),
        symm_s(
            [0.007463627359106, 0.965176994929162],
            [0.050585943594705, 0.749698539312765],
            [0.074613865184212, 0.568983925500818],
        ),
        symm_s_t(
            [0.023501091310143, 0.971086142843168, 0.355832132274584],
            [0.011588562644144, 0.983453947854968, 0.645588139196562],
            [0.023073245798171, 0.933927707027213, 0.821920249234369],
            [0.001570221774472, 1.014086498915039, 0.862185099566557],
            [0.049102258016277, 0.877914842155496, 0.168914072450263],
            [0.042512352239126, 0.882246882640128, 0.568113580166780],
            [0.067270936863160, 0.741324453314596, 0.371360260002223],
            [0.103507336515645, 0.469570217710647, 0.237333359193547],
        ),
    )
    return QuadrilateralScheme("Dunavant 10", weights, points, 21, citation)

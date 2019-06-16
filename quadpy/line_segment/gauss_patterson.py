# -*- coding: utf-8 -*-
#
import numpy

from .gauss_legendre import GaussLegendre
from .tools import integrate
from .helpers import LineSegmentScheme


def GaussPatterson(index):
    """
    Gauss-Patterson quadrature.
    <https://people.sc.fsu.edu/~jburkardt/datasets/quadrature_rules_patterson/quadrature_rules_patterson.html>

    The optimum addition of points to quadrature formulae,
    T.N.L. Patterson,
    Math. Comp. 22 (1968), 847-856,
    <https://doi.org/10.1090/S0025-5718-68-99866-9>.
    """
    degree = 3 * 2 ** index - 1 if index > 0 else 1

    points = numpy.array(_get_points(index))

    # weights
    if index < 6:
        weights = _get_weights(points)
    else:
        # _get_weights is flawed with round-off for index == 6. Use
        # explicit values.
        assert index == 6
        weights = numpy.array(
            [0.2818881418019236e-01]
            + 2 * [0.1680193857410386e-01]
            + 2 * [0.2507856965294977e-01]
            + 2 * [0.6451900050175737e-02]
            + 2 * [0.2739460526398143e-01]
            + 2 * [0.2143898001250387e-01]
            + 2 * [0.1161572331995513e-01]
            + 2 * [0.2108815245726633e-02]
            + 2 * [0.2798921825523816e-01]
            + 2 * [0.2641747339505826e-01]
            + 2 * [0.2340677749531401e-01]
            + 2 * [0.1921990512472777e-01]
            + 2 * [0.1424487737291678e-01]
            + 2 * [0.8989275784064136e-02]
            + 2 * [0.4111503978654693e-02]
            + 2 * [0.6326073193626335e-03]
            + 2 * [0.2813884991562715e-01]
            + 2 * [0.2774070217827968e-01]
            + 2 * [0.2695274966763303e-01]
            + 2 * [0.2579162697602423e-01]
            + 2 * [0.2428216520333660e-01]
            + 2 * [0.2245726582681610e-01]
            + 2 * [0.2035775505847216e-01]
            + 2 * [0.1803221639039129e-01]
            + 2 * [0.1553677555584398e-01]
            + 2 * [0.1293483966360737e-01]
            + 2 * [0.1029711695795636e-01]
            + 2 * [0.7703375233279742e-02]
            + 2 * [0.5249123454808859e-02]
            + 2 * [0.3057753410175531e-02]
            + 2 * [0.1289524082610417e-02]
            + 2 * [0.1807395644453884e-03]
            + 2 * [0.2817631903301660e-01]
            + 2 * [0.2807645579381725e-01]
            + 2 * [0.2787725147661370e-01]
            + 2 * [0.2757974956648187e-01]
            + 2 * [0.2718551322962479e-01]
            + 2 * [0.2669662292745036e-01]
            + 2 * [0.2611567337670610e-01]
            + 2 * [0.2544576996546477e-01]
            + 2 * [0.2469052474448768e-01]
            + 2 * [0.2385405210603854e-01]
            + 2 * [0.2294096422938775e-01]
            + 2 * [0.2195636630531782e-01]
            + 2 * [0.2090585144581202e-01]
            + 2 * [0.1979549504809750e-01]
            + 2 * [0.1863184825613879e-01]
            + 2 * [0.1742193015946417e-01]
            + 2 * [0.1617321872957772e-01]
            + 2 * [0.1489364166481518e-01]
            + 2 * [0.1359157100976555e-01]
            + 2 * [0.1227583056008277e-01]
            + 2 * [0.1095573338783790e-01]
            + 2 * [0.9641177729702537e-02]
            + 2 * [0.8342838753968157e-02]
            + 2 * [0.7072489995433555e-02]
            + 2 * [0.5843449875835640e-02]
            + 2 * [0.4671050372114322e-02]
            + 2 * [0.3572892783517299e-02]
            + 2 * [0.2568764943794020e-02]
            + 2 * [0.1681142865421470e-02]
            + 2 * [0.9383698485423815e-03]
            + 2 * [0.3777466463269846e-03]
            + 2 * [0.5053609520786252e-04]
        )

    return LineSegmentScheme("Gauss-Patterson", degree, weights, points)


def _get_points(index):
    def _pm(a):
        return [a, -a]

    if index == 0:
        return [0.0]

    if index == 1:
        new_points = _pm(0.7745966692414834)
    elif index == 2:
        new_points = _pm(0.4342437493468025) + _pm(0.9604912687080203)
    elif index == 3:
        new_points = (
            _pm(0.2233866864289669)
            + _pm(0.6211029467372264)
            + _pm(0.8884592328722570)
            + _pm(0.9938319632127550)
        )
    elif index == 4:
        new_points = (
            _pm(0.1124889431331866)
            + _pm(0.3311353932579768)
            + _pm(0.5313197436443756)
            + _pm(0.7024962064915271)
            + _pm(0.8367259381688688)
            + _pm(0.9296548574297401)
            + _pm(0.9815311495537401)
            + _pm(0.9990981249676676)
        )
    elif index == 5:
        new_points = (
            _pm(0.5634431304659279e-01)
            + _pm(0.1682352515522075)
            + _pm(0.2777498220218243)
            + _pm(0.3833593241987304)
            + _pm(0.4836180269458411)
            + _pm(0.5771957100520458)
            + _pm(0.6629096600247806)
            + _pm(0.7397560443526947)
            + _pm(0.8069405319502176)
            + _pm(0.8639079381936905)
            + _pm(0.9103711569570043)
            + _pm(0.9463428583734029)
            + _pm(0.9721828747485818)
            + _pm(0.9886847575474295)
            + _pm(0.9972062593722220)
            + _pm(0.9998728881203576)
        )
    else:
        assert index == 6
        new_points = (
            _pm(0.2818464894974569e-01)
            + _pm(0.8445404008371088e-01)
            + _pm(0.1404242331525602)
            + _pm(0.1958975027111002)
            + _pm(0.2506787303034832)
            + _pm(0.3045764415567140)
            + _pm(0.3574038378315322)
            + _pm(0.4089798212298887)
            + _pm(0.4591300119898323)
            + _pm(0.5076877575337166)
            + _pm(0.5544951326319325)
            + _pm(0.5994039302422429)
            + _pm(0.6422766425097595)
            + _pm(0.6829874310910792)
            + _pm(0.7214230853700989)
            + _pm(0.7574839663805136)
            + _pm(0.7910849337998483)
            + _pm(0.8221562543649804)
            + _pm(0.8506444947683502)
            + _pm(0.8765134144847053)
            + _pm(0.8997448997769401)
            + _pm(0.9203400254700124)
            + _pm(0.9383203977795929)
            + _pm(0.9537300064257611)
            + _pm(0.9666378515584165)
            + _pm(0.9771415146397057)
            + _pm(0.9853714995985203)
            + _pm(0.9914957211781061)
            + _pm(0.9957241046984072)
            + _pm(0.9983166353184074)
            + _pm(0.9995987996719107)
            + _pm(0.9999824303548916)
        )

    return _get_points(index - 1) + new_points


def _get_weights(pts):
    """Given a number of points in [-1, 1], according to

        On some Gauss and Lobatto based integration formulae,
        T. N. L. Patterson,
        Math. Comp. 22 (1968), 877-881,

    one can compute the corresponding weights. One reads there:

    > Thus the weights of an n-point integration formula [...] are given by
    >
    >     omega_i = int_{-1}^{1} L_i(x) dx,
    >
    > (where L_i is the Lagrange polynomial for the point x_i).
    > These weights can be evaluated exactly in a numerically stable fashion
    > using a Gauss formula with n/2 points when n is even and (n + 1)/2 points
    > when n is odd.
    """
    n = len(pts)

    # Unnormalized Lagrange polynomial: Degree n, 0 at all x_j except x_i.
    def L(i, x):
        return numpy.prod([(x - pts[j]) for j in range(n) if j != i], axis=0)

    # Gauss-Legendre of order k integrates polynomials of degree 2*k-1 exactly.
    # L has degree n-1, so k needs to be n/2 if n is even, and (n+1)/2 if n is
    # odd.
    k = (n // 2) - 1 if n % 2 == 0 else (n + 1) // 2
    out = numpy.array(
        [
            integrate(
                lambda x, i=i: L(i, x[0]),
                numpy.array([[-1.0], [1.0]]),
                GaussLegendre(k),
            )
            / numpy.prod([(pts[i] - pts[j]) for j in range(n) if j != i])
            for i in range(n)
        ]
    )
    return out

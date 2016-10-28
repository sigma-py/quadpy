# -*- coding: utf-8 -*-
#
import math
import numpy
import sympy


def integrate(f, a, b, scheme):
    out = math.fsum(
        scheme.weights * f(0.5 * (scheme.points.T + 1) * (b-a) + a)
        )
    return 0.5 * (b - a) * out


def show(a, b, scheme):
    from matplotlib import pyplot as plt
    pts = 0.5 * (scheme.points + 1) * (b-a) + a
    plt.plot([0.0, 1.0], [0.0, 0.0], '-k')
    plt.bar(
        pts, scheme.weights,
        color='b',
        alpha=0.5,
        width=(b-a) / len(scheme.weights)
        )
    return


class Midpoint(object):
    def __init__(self):
        self.weights = [2.0]
        self.points = numpy.array([
            0.0
            ])
        self.degree = 1
        return


class Trapezoidal(object):
    def __init__(self):
        self.weights = [1.0, 1.0]
        self.points = numpy.array([
            -1.0,
            1.0
            ])
        self.degree = 1
        return


class ChebyshevGauss1(object):
    '''
    Chebyshev-Gauß quadrature for \int_{-1}^1 f(x) / sqrt(1+x^2) dx.
    '''
    def __init__(self, n):
        self.degree = n if n % 2 == 1 else n+1
        self.points = numpy.cos(
                (2*numpy.arange(1, n+1) - 1.0) / (2*n) * numpy.pi
                )
        self.weights = numpy.pi / n * numpy.ones(n)
        return


class ChebyshevGauss2(object):
    '''
    Chebyshev-Gauß quadrature for \int_{-1}^1 f(x) * sqrt(1+x^2) dx.
    '''
    def __init__(self, n):
        self.degree = n if n % 2 == 1 else n+1
        self.points = numpy.cos(
                numpy.pi * numpy.arange(1, n+1) / (n+1)
                )
        self.weights = numpy.pi / (n+1) \
            * (numpy.sin(numpy.pi * numpy.arange(1, n+1) / (n+1)))**2
        return


class GaussHermite(object):
    '''
    Gauß-Hermite quadrature.
    '''
    def __init__(self, n):
        self.degree = 2*n - 1
        self.points, self.weights = numpy.polynomial.hermite.hermgauss(n)
        return


class GaussLaguerre(object):
    '''
    Gauß-Laguerre quadrature.
    '''
    def __init__(self, n):
        self.degree = 2*n - 1
        self.points, self.weights = numpy.polynomial.laguerre.laggauss(n)
        return


class GaussLegendre(object):
    '''
    Gauß-Legendre quadrature.
    '''
    def __init__(self, n):
        self.degree = 2*n - 1
        self.points, self.weights = numpy.polynomial.legendre.leggauss(n)
        return


def _jacobi_recursion_coefficients(n, a, b):
    '''
    Generate the recursion coefficients alpha_k, beta_k

    P_{k+1}(x) = (x-alpha_k)*P_{k}(x) - beta_k P_{k-1}(x)

    for the Jacobi polynomials which are orthogonal on [-1,1]
    with respect to the weight w(x)=[(1-x)^a]*[(1+x)^b].

    Adapted from the MATLAB code by Dirk Laurie and Walter Gautschi
    http://www.cs.purdue.edu/archives/2002/wxg/codes/r_jacobi.m
    and from Greg van Winckel's
    https://github.com/gregvw/orthopoly-quadrature/blob/master/rec_jacobi.pyx
    '''
    assert a > -1.0 or b > -1.0
    assert n >= 1

    mu = 2.0**(a+b+1.0) \
        * numpy.exp(
            math.lgamma(a+1.0) + math.lgamma(b+1.0) - math.lgamma(a+b+2.0)
            )
    nu = (b-a) / (a+b+2.0)

    if n == 1:
        return nu, mu

    N = numpy.arange(1, n)

    nab = 2.0*N + a + b
    alpha = numpy.hstack((nu, (b**2 - a**2) / (nab * (nab + 2.0))))
    N = N[1:]
    nab = nab[1:]
    B1 = 4.0 * (a+1.0) * (b+1.0) / ((a+b+2.0)**2.0 * (a+b+3.0))
    B = 4.0 * (N+a) * (N+b) * N * (N+a+b) / (nab**2.0 * (nab+1.0) * (nab-1.0))
    beta = numpy.hstack((mu, B1, B))
    return alpha, beta


def _gauss(alpha, beta):
    '''
    Compute the Gauss nodes and weights from the recursion
    coefficients associated with a set of orthogonal polynomials

    Adapted from the MATLAB code by Walter Gautschi
    http://www.cs.purdue.edu/archives/2002/wxg/codes/gauss.m

    and

    http://www.scientificpython.net/pyblog/radau-quadrature
    '''
    from scipy.linalg import eig_banded
    import scipy as sp
    A = numpy.vstack((numpy.sqrt(beta), alpha))
    x, V = eig_banded(A, lower=False)
    w = beta[0]*sp.real(sp.power(V[0, :], 2))
    return x, w


def _lobatto(alpha, beta, xl1, xl2):
    ''' Compute the Lobatto nodes and weights with the preassigned
        node xl1,xl2

        Based on the section 7 of the paper 'Some modified matrix eigenvalue
        problems' by Gene Golub, SIAM Review Vol 15, No. 2, April 1973,
        pp.318--334

        and

        http://www.scientificpython.net/pyblog/radau-quadrature
    '''
    from scipy.linalg import solve_banded, solve
    n = len(alpha)-1
    en = numpy.zeros(n)
    en[-1] = 1
    A1 = numpy.vstack((numpy.sqrt(beta), alpha-xl1))
    J1 = numpy.vstack((A1[:, 0:-1], A1[0, 1:]))
    A2 = numpy.vstack((numpy.sqrt(beta), alpha-xl2))
    J2 = numpy.vstack((A2[:, 0:-1], A2[0, 1:]))
    g1 = solve_banded((1, 1), J1, en)
    g2 = solve_banded((1, 1), J2, en)
    C = numpy.array(((1, -g1[-1]), (1, -g2[-1])))
    xl = numpy.array((xl1, xl2))
    ab = solve(C, xl)

    alphal = alpha
    alphal[-1] = ab[0]
    betal = beta
    betal[-1] = ab[1]
    x, w = _gauss(alphal, betal)
    return x, w


class GaussLobatto(object):
    '''
    Gauß-Lobatto quadrature.
    '''
    def __init__(self, n, a=0.0, b=0.0):
        assert n >= 2
        self.degree = 2*n - 3
        alpha, beta = _jacobi_recursion_coefficients(n, a, b)
        self.points, self.weights = _lobatto(alpha, beta, -1.0, 1.0)
        return


class GaussPatterson(object):
    '''
    Gauß-Patterson quadrature.
    <https://people.sc.fsu.edu/~jburkardt/datasets/quadrature_rules_patterson/quadrature_rules_patterson.html>
    '''
    def __init__(self, index):
        self.degree = 3*2**index - 1
        if index == 0:
            self.degree = 1  # override degree
            self.weights = [2.0]
            self.points = numpy.array([
                0.0
                ])
        elif index == 1:
            self.weights = numpy.array([
                5.0 / 9.0,
                8.0 / 9.0,
                5.0 / 9.0,
                ])
            self.points = numpy.array([
                -0.7745966692414834,
                0.0,
                0.7745966692414834,
                ])
        elif index == 2:
            self.weights = numpy.array([
                0.1046562260264673,
                0.2684880898683334,
                0.4013974147759622,
                0.4509165386584741,
                0.4013974147759622,
                0.2684880898683334,
                0.1046562260264673,
                ])
            self.points = numpy.array([
                -0.9604912687080203,
                -0.7745966692414834,
                -0.4342437493468025,
                0.0,
                0.4342437493468025,
                0.7745966692414834,
                0.9604912687080203,
                ])
        elif index == 3:
            self.weights = numpy.array([
                0.1700171962994026E-01,
                0.5160328299707974E-01,
                0.9292719531512454E-01,
                0.1344152552437842,
                0.1715119091363914,
                0.2006285293769890,
                0.2191568584015875,
                0.2255104997982067,
                0.2191568584015875,
                0.2006285293769890,
                0.1715119091363914,
                0.1344152552437842,
                0.9292719531512454E-01,
                0.5160328299707974E-01,
                0.1700171962994026E-01,
                ])
            self.points = numpy.array([
                -0.9938319632127550,
                -0.9604912687080203,
                -0.8884592328722570,
                -0.7745966692414834,
                -0.6211029467372264,
                -0.4342437493468025,
                -0.2233866864289669,
                0.0,
                0.2233866864289669,
                0.4342437493468025,
                0.6211029467372264,
                0.7745966692414834,
                0.8884592328722570,
                0.9604912687080203,
                0.9938319632127550,
                ])
        elif index == 4:
            self.weights = numpy.array([
                0.2544780791561875E-02,
                0.8434565739321106E-02,
                0.1644604985438781E-01,
                0.2580759809617665E-01,
                0.3595710330712932E-01,
                0.4646289326175799E-01,
                0.5697950949412336E-01,
                0.6720775429599070E-01,
                0.7687962049900353E-01,
                0.8575592004999034E-01,
                0.9362710998126447E-01,
                0.1003142786117956,
                0.1056698935802348,
                0.1095784210559246,
                0.1119568730209535,
                0.1127552567207687,
                0.1119568730209535,
                0.1095784210559246,
                0.1056698935802348,
                0.1003142786117956,
                0.9362710998126447E-01,
                0.8575592004999034E-01,
                0.7687962049900353E-01,
                0.6720775429599070E-01,
                0.5697950949412336E-01,
                0.4646289326175799E-01,
                0.3595710330712932E-01,
                0.2580759809617665E-01,
                0.1644604985438781E-01,
                0.8434565739321106E-02,
                0.2544780791561875E-02,
                ])
            self.points = numpy.array([
                -0.9990981249676676,
                -0.9938319632127550,
                -0.9815311495537401,
                -0.9604912687080203,
                -0.9296548574297401,
                -0.8884592328722570,
                -0.8367259381688688,
                -0.7745966692414834,
                -0.7024962064915271,
                -0.6211029467372264,
                -0.5313197436443756,
                -0.4342437493468025,
                -0.3311353932579768,
                -0.2233866864289669,
                -0.1124889431331866,
                0.0,
                0.1124889431331866,
                0.2233866864289669,
                0.3311353932579768,
                0.4342437493468025,
                0.5313197436443756,
                0.6211029467372264,
                0.7024962064915271,
                0.7745966692414834,
                0.8367259381688688,
                0.8884592328722570,
                0.9296548574297401,
                0.9604912687080203,
                0.9815311495537401,
                0.9938319632127550,
                0.9990981249676676,
                ])
        elif index == 5:
            self.weights = numpy.array([
                0.3632214818455306E-03,
                0.1265156556230068E-02,
                0.2579049794685688E-02,
                0.4217630441558855E-02,
                0.6115506822117246E-02,
                0.8223007957235930E-02,
                0.1049824690962132E-01,
                0.1290380010035127E-01,
                0.1540675046655950E-01,
                0.1797855156812827E-01,
                0.2059423391591271E-01,
                0.2323144663991027E-01,
                0.2586967932721475E-01,
                0.2848975474583355E-01,
                0.3107355111168797E-01,
                0.3360387714820773E-01,
                0.3606443278078257E-01,
                0.3843981024945553E-01,
                0.4071551011694432E-01,
                0.4287796002500773E-01,
                0.4491453165363220E-01,
                0.4681355499062801E-01,
                0.4856433040667320E-01,
                0.5015713930589954E-01,
                0.5158325395204846E-01,
                0.5283494679011652E-01,
                0.5390549933526606E-01,
                0.5478921052796287E-01,
                0.5548140435655936E-01,
                0.5597843651047632E-01,
                0.5627769983125430E-01,
                0.5637762836038471E-01,
                0.5627769983125430E-01,
                0.5597843651047632E-01,
                0.5548140435655936E-01,
                0.5478921052796287E-01,
                0.5390549933526606E-01,
                0.5283494679011652E-01,
                0.5158325395204846E-01,
                0.5015713930589954E-01,
                0.4856433040667320E-01,
                0.4681355499062801E-01,
                0.4491453165363220E-01,
                0.4287796002500773E-01,
                0.4071551011694432E-01,
                0.3843981024945553E-01,
                0.3606443278078257E-01,
                0.3360387714820773E-01,
                0.3107355111168797E-01,
                0.2848975474583355E-01,
                0.2586967932721475E-01,
                0.2323144663991027E-01,
                0.2059423391591271E-01,
                0.1797855156812827E-01,
                0.1540675046655950E-01,
                0.1290380010035127E-01,
                0.1049824690962132E-01,
                0.8223007957235930E-02,
                0.6115506822117246E-02,
                0.4217630441558855E-02,
                0.2579049794685688E-02,
                0.1265156556230068E-02,
                0.3632214818455306E-03,
                ])
            self.points = numpy.array([
                -0.9998728881203576,
                -0.9990981249676676,
                -0.9972062593722220,
                -0.9938319632127550,
                -0.9886847575474295,
                -0.9815311495537401,
                -0.9721828747485818,
                -0.9604912687080203,
                -0.9463428583734029,
                -0.9296548574297401,
                -0.9103711569570043,
                -0.8884592328722570,
                -0.8639079381936905,
                -0.8367259381688688,
                -0.8069405319502176,
                -0.7745966692414834,
                -0.7397560443526947,
                -0.7024962064915271,
                -0.6629096600247806,
                -0.6211029467372264,
                -0.5771957100520458,
                -0.5313197436443756,
                -0.4836180269458411,
                -0.4342437493468025,
                -0.3833593241987304,
                -0.3311353932579768,
                -0.2777498220218243,
                -0.2233866864289669,
                -0.1682352515522075,
                -0.1124889431331866,
                -0.5634431304659279E-01,
                0.0,
                0.5634431304659279E-01,
                0.1124889431331866,
                0.1682352515522075,
                0.2233866864289669,
                0.2777498220218243,
                0.3311353932579768,
                0.3833593241987304,
                0.4342437493468025,
                0.4836180269458411,
                0.5313197436443756,
                0.5771957100520458,
                0.6211029467372264,
                0.6629096600247806,
                0.7024962064915271,
                0.7397560443526947,
                0.7745966692414834,
                0.8069405319502176,
                0.8367259381688688,
                0.8639079381936905,
                0.8884592328722570,
                0.9103711569570043,
                0.9296548574297401,
                0.9463428583734029,
                0.9604912687080203,
                0.9721828747485818,
                0.9815311495537401,
                0.9886847575474295,
                0.9938319632127550,
                0.9972062593722220,
                0.9990981249676676,
                0.9998728881203576,
                ])
        elif index == 6:
            self.weights = numpy.array([
                0.5053609520786252E-04,
                0.1807395644453884E-03,
                0.3777466463269846E-03,
                0.6326073193626335E-03,
                0.9383698485423815E-03,
                0.1289524082610417E-02,
                0.1681142865421470E-02,
                0.2108815245726633E-02,
                0.2568764943794020E-02,
                0.3057753410175531E-02,
                0.3572892783517299E-02,
                0.4111503978654693E-02,
                0.4671050372114322E-02,
                0.5249123454808859E-02,
                0.5843449875835640E-02,
                0.6451900050175737E-02,
                0.7072489995433555E-02,
                0.7703375233279742E-02,
                0.8342838753968157E-02,
                0.8989275784064136E-02,
                0.9641177729702537E-02,
                0.1029711695795636E-01,
                0.1095573338783790E-01,
                0.1161572331995513E-01,
                0.1227583056008277E-01,
                0.1293483966360737E-01,
                0.1359157100976555E-01,
                0.1424487737291678E-01,
                0.1489364166481518E-01,
                0.1553677555584398E-01,
                0.1617321872957772E-01,
                0.1680193857410386E-01,
                0.1742193015946417E-01,
                0.1803221639039129E-01,
                0.1863184825613879E-01,
                0.1921990512472777E-01,
                0.1979549504809750E-01,
                0.2035775505847216E-01,
                0.2090585144581202E-01,
                0.2143898001250387E-01,
                0.2195636630531782E-01,
                0.2245726582681610E-01,
                0.2294096422938775E-01,
                0.2340677749531401E-01,
                0.2385405210603854E-01,
                0.2428216520333660E-01,
                0.2469052474448768E-01,
                0.2507856965294977E-01,
                0.2544576996546477E-01,
                0.2579162697602423E-01,
                0.2611567337670610E-01,
                0.2641747339505826E-01,
                0.2669662292745036E-01,
                0.2695274966763303E-01,
                0.2718551322962479E-01,
                0.2739460526398143E-01,
                0.2757974956648187E-01,
                0.2774070217827968E-01,
                0.2787725147661370E-01,
                0.2798921825523816E-01,
                0.2807645579381725E-01,
                0.2813884991562715E-01,
                0.2817631903301660E-01,
                0.2818881418019236E-01,
                0.2817631903301660E-01,
                0.2813884991562715E-01,
                0.2807645579381725E-01,
                0.2798921825523816E-01,
                0.2787725147661370E-01,
                0.2774070217827968E-01,
                0.2757974956648187E-01,
                0.2739460526398143E-01,
                0.2718551322962479E-01,
                0.2695274966763303E-01,
                0.2669662292745036E-01,
                0.2641747339505826E-01,
                0.2611567337670610E-01,
                0.2579162697602423E-01,
                0.2544576996546477E-01,
                0.2507856965294977E-01,
                0.2469052474448768E-01,
                0.2428216520333660E-01,
                0.2385405210603854E-01,
                0.2340677749531401E-01,
                0.2294096422938775E-01,
                0.2245726582681610E-01,
                0.2195636630531782E-01,
                0.2143898001250387E-01,
                0.2090585144581202E-01,
                0.2035775505847216E-01,
                0.1979549504809750E-01,
                0.1921990512472777E-01,
                0.1863184825613879E-01,
                0.1803221639039129E-01,
                0.1742193015946417E-01,
                0.1680193857410386E-01,
                0.1617321872957772E-01,
                0.1553677555584398E-01,
                0.1489364166481518E-01,
                0.1424487737291678E-01,
                0.1359157100976555E-01,
                0.1293483966360737E-01,
                0.1227583056008277E-01,
                0.1161572331995513E-01,
                0.1095573338783790E-01,
                0.1029711695795636E-01,
                0.9641177729702537E-02,
                0.8989275784064136E-02,
                0.8342838753968157E-02,
                0.7703375233279742E-02,
                0.7072489995433555E-02,
                0.6451900050175737E-02,
                0.5843449875835640E-02,
                0.5249123454808859E-02,
                0.4671050372114322E-02,
                0.4111503978654693E-02,
                0.3572892783517299E-02,
                0.3057753410175531E-02,
                0.2568764943794020E-02,
                0.2108815245726633E-02,
                0.1681142865421470E-02,
                0.1289524082610417E-02,
                0.9383698485423815E-03,
                0.6326073193626335E-03,
                0.3777466463269846E-03,
                0.1807395644453884E-03,
                0.5053609520786252E-04,
                ])
            self.points = numpy.array([
                -0.9999824303548916,
                -0.9998728881203576,
                -0.9995987996719107,
                -0.9990981249676676,
                -0.9983166353184074,
                -0.9972062593722220,
                -0.9957241046984072,
                -0.9938319632127550,
                -0.9914957211781061,
                -0.9886847575474295,
                -0.9853714995985203,
                -0.9815311495537401,
                -0.9771415146397057,
                -0.9721828747485818,
                -0.9666378515584165,
                -0.9604912687080203,
                -0.9537300064257611,
                -0.9463428583734029,
                -0.9383203977795929,
                -0.9296548574297401,
                -0.9203400254700124,
                -0.9103711569570043,
                -0.8997448997769401,
                -0.8884592328722570,
                -0.8765134144847053,
                -0.8639079381936905,
                -0.8506444947683502,
                -0.8367259381688688,
                -0.8221562543649804,
                -0.8069405319502176,
                -0.7910849337998483,
                -0.7745966692414834,
                -0.7574839663805136,
                -0.7397560443526947,
                -0.7214230853700989,
                -0.7024962064915271,
                -0.6829874310910792,
                -0.6629096600247806,
                -0.6422766425097595,
                -0.6211029467372264,
                -0.5994039302422429,
                -0.5771957100520458,
                -0.5544951326319325,
                -0.5313197436443756,
                -0.5076877575337166,
                -0.4836180269458411,
                -0.4591300119898323,
                -0.4342437493468025,
                -0.4089798212298887,
                -0.3833593241987304,
                -0.3574038378315322,
                -0.3311353932579768,
                -0.3045764415567140,
                -0.2777498220218243,
                -0.2506787303034832,
                -0.2233866864289669,
                -0.1958975027111002,
                -0.1682352515522075,
                -0.1404242331525602,
                -0.1124889431331866,
                -0.8445404008371088E-01,
                -0.5634431304659279E-01,
                -0.2818464894974569E-01,
                0.0,
                0.2818464894974569E-01,
                0.5634431304659279E-01,
                0.8445404008371088E-01,
                0.1124889431331866,
                0.1404242331525602,
                0.1682352515522075,
                0.1958975027111002,
                0.2233866864289669,
                0.2506787303034832,
                0.2777498220218243,
                0.3045764415567140,
                0.3311353932579768,
                0.3574038378315322,
                0.3833593241987304,
                0.4089798212298887,
                0.4342437493468025,
                0.4591300119898323,
                0.4836180269458411,
                0.5076877575337166,
                0.5313197436443756,
                0.5544951326319325,
                0.5771957100520458,
                0.5994039302422429,
                0.6211029467372264,
                0.6422766425097595,
                0.6629096600247806,
                0.6829874310910792,
                0.7024962064915271,
                0.7214230853700989,
                0.7397560443526947,
                0.7574839663805136,
                0.7745966692414834,
                0.7910849337998483,
                0.8069405319502176,
                0.8221562543649804,
                0.8367259381688688,
                0.8506444947683502,
                0.8639079381936905,
                0.8765134144847053,
                0.8884592328722570,
                0.8997448997769401,
                0.9103711569570043,
                0.9203400254700124,
                0.9296548574297401,
                0.9383203977795929,
                0.9463428583734029,
                0.9537300064257611,
                0.9604912687080203,
                0.9666378515584165,
                0.9721828747485818,
                0.9771415146397057,
                0.9815311495537401,
                0.9853714995985203,
                0.9886847575474295,
                0.9914957211781061,
                0.9938319632127550,
                0.9957241046984072,
                0.9972062593722220,
                0.9983166353184074,
                0.9990981249676676,
                0.9995987996719107,
                0.9998728881203576,
                0.9999824303548916,
                ])
        else:
            raise ValueError('Illegal Gauss-Patterson order')


class ClenshawCurtis(object):
    '''
    Clenshaw-Curtis quadrature.

    Weights are constructed after

    J. Waldvogel,
    Fast Construction of the Fejér and Clenshaw–Curtis Quadrature Rules,
    BIT Numerical Mathematics, March 2006, Volume 46, Issue 1, pp 195–202,
    DOI: 10.1007/s10543-006-0045-4,
    <https://dx.doi.org/10.1007/s10543-006-0045-4>.
    '''
    def __init__(self, n):
        self.degree = n

        self.points = -numpy.cos((numpy.pi * numpy.arange(n)) / (n-1))

        if n == 2:
            self.weights = numpy.array([1.0, 1.0])
            return

        n -= 1
        N = numpy.arange(1, n, 2)
        l = len(N)
        m = n - l
        v0 = numpy.concatenate([
            2.0 / N / (N-2),
            numpy.array([1.0 / N[-1]]),
            numpy.zeros(m),
            ])
        v2 = - v0[:-1] - v0[:0:-1]
        g0 = -numpy.ones(n)
        g0[l] += n
        g0[m] += n
        g = g0 / (n**2 - 1 + (n % 2))

        w = numpy.fft.ihfft(v2 + g)
        assert max(w.imag) < 1.0e-15
        w = w.real

        if n % 2 == 1:
            self.weights = numpy.concatenate([
                w,
                w[::-1]
                ])
        else:
            self.weights = numpy.concatenate([
                w,
                w[len(w)-2::-1]
                ])

        return


class Fejer1(object):
    '''
    Fejér-type-1 quadrature.

    Weights are constructed after

    J. Waldvogel,
    Fast Construction of the Fejér and Clenshaw–Curtis Quadrature Rules,
    BIT Numerical Mathematics, March 2006, Volume 46, Issue 1, pp 195–202,
    DOI: 10.1007/s10543-006-0045-4,
    <https://dx.doi.org/10.1007/s10543-006-0045-4>.
    '''
    def __init__(self, n):
        self.degree = n

        self.points = -numpy.cos(numpy.pi * (numpy.arange(n) + 0.5) / n)

        # n -= 1
        N = numpy.arange(1, n, 2)
        l = len(N)
        m = n - l
        K = numpy.arange(m)

        v0 = numpy.concatenate([
            2 * numpy.exp(1j*numpy.pi*K/n) / (1 - 4*K**2),
            numpy.zeros(l+1)
            ])
        v1 = v0[:-1] + numpy.conjugate(v0[:0:-1])

        w = numpy.fft.ifft(v1)
        assert max(w.imag) < 1.0e-15
        self.weights = w.real

        return


class Fejer2(object):
    '''
    Fejér-type-2 quadrature.

    Weights are constructed after

    J. Waldvogel,
    Fast Construction of the Fejér and Clenshaw–Curtis Quadrature Rules,
    BIT Numerical Mathematics, March 2006, Volume 46, Issue 1, pp 195–202,
    DOI: 10.1007/s10543-006-0045-4,
    <https://dx.doi.org/10.1007/s10543-006-0045-4>.
    '''
    def __init__(self, n):
        self.degree = n

        self.points = -numpy.cos((numpy.pi * numpy.arange(1, n+1)) / (n+1))

        n += 1
        N = numpy.arange(1, n, 2)
        l = len(N)
        m = n - l
        v0 = numpy.concatenate([
            2.0 / N / (N-2),
            numpy.array([1.0 / N[-1]]),
            numpy.zeros(m),
            ])
        v2 = - v0[:-1] - v0[:0:-1]

        w = numpy.fft.ihfft(v2)
        assert max(w.imag) < 1.0e-15
        w = w.real

        if n % 2 == 1:
            self.weights = numpy.concatenate([w, w[::-1]])
        else:
            self.weights = numpy.concatenate([w, w[len(w)-2::-1]])

        # cut off first and last
        self.weights = self.weights[1:-1]

        return


class NewtonCotesClosed(object):
    '''
    Closed Newton-Cotes formulae.
    <https://en.wikipedia.org/wiki/Newton%E2%80%93Cotes_formulas#Closed_Newton.E2.80.93Cotes_formulae>,
    <http://mathworld.wolfram.com/Newton-CotesFormulas.html>.
    '''
    def __init__(self, index):
        self.points = numpy.linspace(-1.0, 1.0, index+1)
        self.degree = index + 1 if index % 2 == 0 else index

        # Formula (26) from
        # <http://mathworld.wolfram.com/Newton-CotesFormulas.html>.
        # Note that Sympy carries out all operations in rationals, i.e.,
        # _exactly_. Only at the end, the rational is converted into a float.
        n = index
        self.weights = numpy.empty(n+1)
        for r in range(n+1):
            t = sympy.Symbol('t')
            f = 1
            for i in range(n+1):
                if i != r:
                    f *= (t - i)
            alpha = 2 * \
                (-1)**(n-r) * sympy.integrate(f, (t, 0, n)) \
                / (math.factorial(r) * math.factorial(n-r)) \
                / index
            self.weights[r] = alpha

        return


class NewtonCotesOpen(object):
    '''
    Open Newton-Cotes formulae.
    <http://math.stackexchange.com/a/1959071/36678>
    '''
    def __init__(self, index):
        self.points = numpy.linspace(-1.0, 1.0, index+2)[1:-1]
        self.degree = index if (index+1) % 2 == 0 else index - 1
        #
        n = index+1
        self.weights = numpy.empty(n-1)
        for r in range(1, n):
            t = sympy.Symbol('t')
            f = 1
            for i in range(1, n):
                if i != r:
                    f *= (t - i)
            alpha = 2 * \
                (-1)**(n-r+1) * sympy.integrate(f, (t, 0, n)) \
                / (math.factorial(r-1) * math.factorial(n-1-r)) \
                / n
            self.weights[r-1] = alpha
        return

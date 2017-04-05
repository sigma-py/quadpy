# -*- coding: utf-8 -*-
#
from . import helpers

import math
import numpy
import sympy


def integrate(f, interval, scheme, sumfun=helpers.kahan_sum):
    alpha = 0.5 * (interval[1] - interval[0])
    beta = 0.5 * (interval[0] + interval[1])
    out = sumfun(
        # use multiply.outer instead of outer here. The handles floats in the
        # second argument such that no extra dimension is produces.
        scheme.weights *
        f((numpy.multiply.outer(scheme.points, alpha) + beta).T),
        axis=-1
        )
    return alpha * out


def adaptive_integrate(f, interval, eps, sumfun=helpers.kahan_sum):
    # Use Gauss-Kronrod 7/15 schemes for error estimation and adaptivity.
    # TODO merge
    scheme1 = GaussLegendre(7)
    scheme2 = GaussKronrod(7)
    val1 = integrate(f, interval, scheme1, sumfun=sumfun)
    val2 = integrate(f, interval, scheme2, sumfun=sumfun)

    # TODO better error estimation
    error_estimation = abs(val1 - val2)

    if error_estimation < eps:
        return val1

    # mark the only interval as bad
    quad_sum = 0.0
    intervals = numpy.array([interval])
    is_bad = [True]

    original_interval_length = abs(interval[1] - interval[0])

    while any(is_bad):
        # split the intervals in half
        midpoints = \
            0.5 * (intervals[is_bad][:, 0] + intervals[is_bad][:, 1])
        intervals = numpy.concatenate([
            numpy.column_stack([intervals[is_bad][:, 0], midpoints]),
            numpy.column_stack([midpoints, intervals[is_bad][:, 1]]),
            ])
        # compute values and error estimates for the new intervals
        val1 = integrate(f, intervals.T, scheme1, sumfun=sumfun)
        val2 = integrate(f, intervals.T, scheme2, sumfun=sumfun)
        # mark bad intervals
        error_estimates = abs(val1 - val2)
        lengths = abs(intervals[:, 1] - intervals[:, 0])
        is_bad = error_estimates > eps * lengths / original_interval_length
        # add values from good intervals to sum
        is_good = numpy.logical_not(is_bad)
        quad_sum += sumfun(val1[is_good])

    return quad_sum


def show(scheme, interval=numpy.array([-1.0, 1.0]), show_axes=False):
    from matplotlib import pyplot as plt
    # change default range so that new disks will work
    plt.axis('equal')
    # ax.set_xlim((-1.5, 1.5))
    # ax.set_ylim((-1.5, 1.5))

    if not show_axes:
        plt.gca().set_axis_off()

    plt.plot(interval, [0, 0], color='k')

    pts = numpy.column_stack([scheme.points, numpy.zeros(len(scheme.points))])

    # The total area is used to gauge the disk radii. This is only meaningful
    # for 2D manifolds, not for the circle. What we do instead is choose the
    # total_area such that the sum of the disk radii equals b-a.
    length = interval[1] - interval[0]
    total_area = 0.25 * length**2 * numpy.pi * sum(scheme.weights) \
        / sum(numpy.sqrt(abs(scheme.weights)))**2

    helpers.plot_disks(
        plt, pts, scheme.weights, total_area
        )
    plt.show()
    return


class Midpoint(object):
    def __init__(self):
        self.weights = numpy.array([2.0])
        self.points = numpy.array([
            0.0
            ])
        self.degree = 1
        return


class Trapezoidal(object):
    def __init__(self):
        self.weights = numpy.array([1.0, 1.0])
        self.points = numpy.array([
            -1.0,
            1.0
            ])
        self.degree = 1
        return


class ChebyshevGauss1(object):
    '''
    Chebyshev-Gauss quadrature for \int_{-1}^1 f(x) / sqrt(1+x^2) dx.
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
    Chebyshev-Gauss quadrature for \int_{-1}^1 f(x) * sqrt(1+x^2) dx.
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
    Gauss-Hermite quadrature.
    '''
    def __init__(self, n):
        self.degree = 2*n - 1
        self.points, self.weights = numpy.polynomial.hermite.hermgauss(n)
        return


class GaussLaguerre(object):
    '''
    Gauss-Laguerre quadrature.
    '''
    def __init__(self, n):
        self.degree = 2*n - 1
        self.points, self.weights = numpy.polynomial.laguerre.laggauss(n)
        return


class GaussLegendre(object):
    '''
    Gauss-Legendre quadrature.
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
    alpha = numpy.hstack([nu, (b**2 - a**2) / (nab * (nab + 2.0))])
    N = N[1:]
    nab = nab[1:]
    B1 = 4.0 * (a+1.0) * (b+1.0) / ((a+b+2.0)**2.0 * (a+b+3.0))
    B = 4.0 * (N+a) * (N+b) * N * (N+a+b) / (nab**2.0 * (nab+1.0) * (nab-1.0))
    beta = numpy.hstack((mu, B1, B))
    return alpha, beta


def _gauss(alpha, beta):
    '''
    Compute the Gauss nodes and weights from the recursion coefficients
    associated with a set of orthogonal polynomials.

    Algorithm 726: ORTHPOL–a package of routines for generating orthogonal
    polynomials and Gauss-type quadrature rules,
    W. Gautschi,
    ACM Transactions on Mathematical Software (TOMS),
    Volume 20, Issue 1, March 1994,
    Pages 21-62,
    <http://doi.org/10.1145/174603.174605>,
    <http://www.cs.purdue.edu/archives/2002/wxg/codes/gauss.m>,

    and

    http://www.scientificpython.net/pyblog/radau-quadrature
    '''
    from scipy.linalg import eig_banded
    import scipy
    A = numpy.vstack((numpy.sqrt(beta), alpha))
    x, V = eig_banded(A, lower=False)
    w = beta[0]*scipy.real(scipy.power(V[0, :], 2))
    return x, w


def _lobatto(alpha, beta, xl1, xl2):
    '''Compute the Lobatto nodes and weights with the preassigned node xl1, xl2.
    Based on the section 7 of the paper

        Some modified matrix eigenvalue problems,
        Gene Golub,
        SIAM Review Vol 15, No. 2, April 1973, pp.318--334,

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


def _radau(alpha, beta, xr):
    '''From <http://www.scientificpython.net/pyblog/radau-quadrature>:
    Compute the Radau nodes and weights with the preassigned node xr.

    Based on the section 7 of the paper

        Some modified matrix eigenvalue problems,
        Gene Golub,
        SIAM Review Vol 15, No. 2, April 1973, pp.318--334.
    '''
    from scipy.linalg import solve_banded

    n = len(alpha)-1
    f = numpy.zeros(n)
    f[-1] = beta[-1]
    A = numpy.vstack((numpy.sqrt(beta), alpha-xr))
    J = numpy.vstack((A[:, 0:-1], A[0, 1:]))
    delta = solve_banded((1, 1), J, f)
    alphar = alpha.copy()
    alphar[-1] = xr + delta[-1]
    x, w = _gauss(alphar, beta)
    return x, w


class GaussLobatto(object):
    '''
    Gauss-Lobatto quadrature.
    '''
    def __init__(self, n, a=0.0, b=0.0):
        assert n >= 2
        self.degree = 2*n - 3
        alpha, beta = _jacobi_recursion_coefficients(n, a, b)
        self.points, self.weights = _lobatto(alpha, beta, -1.0, 1.0)
        return


class GaussRadau(object):
    '''
    Gauss-Radau quadrature.
    '''
    def __init__(self, n, a=0.0, b=0.0):
        assert n >= 2
        self.degree = 2*n - 1
        alpha, beta = _jacobi_recursion_coefficients(n, a, b)
        self.points, self.weights = _radau(alpha, beta, -1.0)
        return


def _get_weights(pts):
    '''Given a number of points in [-1, 1], according to

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
    '''
    n = len(pts)

    # Unnormalized Lagrange polynomial: Degree n, 0 at all x_j except x_i.
    def L(i, x):
        return numpy.prod([(x - pts[j]) for j in range(n) if j != i], axis=0)

    # Gauss-Legendre of order k integrates polynomials of degree 2*k-1 exactly.
    # L has degree n-1, so k needs to be n/2 if n is even, and (n+1)/2 if n is
    # odd.
    k = (n // 2) - 1 if n % 2 == 0 else (n+1) // 2
    return numpy.array([
        integrate(
            lambda x, i=i: L(i, x), [-1.0, 1.0], GaussLegendre(k),
            sumfun=lambda a, axis: numpy.array([math.fsum(a)])
            )[0]
        /
        numpy.prod([(pts[i] - pts[j]) for j in range(n) if j != i])
        for i in range(n)
        ])


class GaussPatterson(object):
    '''
    Gauss-Patterson quadrature.
    <https://people.sc.fsu.edu/~jburkardt/datasets/quadrature_rules_patterson/quadrature_rules_patterson.html>

    The optimum addition of points to quadrature formulae,
    T.N.L. Patterson,
    Math. Comp. 22 (1968), 847-856,
    <https://doi.org/10.1090/S0025-5718-68-99866-9>.
    '''
    def __init__(self, index):
        self.degree = 3*2**index - 1 if index > 0 else 1

        self.points = numpy.array(self._get_points(index))

        # weights
        if index < 6:
            self.weights = _get_weights(self.points)
        else:
            # _get_weights is flawed with round-off for index == 6. Use
            # explicit values.
            assert index == 6
            self.weights = numpy.array(
                [0.2818881418019236E-01] +
                #
                2 * [0.1680193857410386E-01] +
                #
                2 * [0.2507856965294977E-01] +
                2 * [0.6451900050175737E-02] +
                #
                2 * [0.2739460526398143E-01] +
                2 * [0.2143898001250387E-01] +
                2 * [0.1161572331995513E-01] +
                2 * [0.2108815245726633E-02] +
                #
                2 * [0.2798921825523816E-01] +
                2 * [0.2641747339505826E-01] +
                2 * [0.2340677749531401E-01] +
                2 * [0.1921990512472777E-01] +
                2 * [0.1424487737291678E-01] +
                2 * [0.8989275784064136E-02] +
                2 * [0.4111503978654693E-02] +
                2 * [0.6326073193626335E-03] +
                #
                2 * [0.2813884991562715E-01] +
                2 * [0.2774070217827968E-01] +
                2 * [0.2695274966763303E-01] +
                2 * [0.2579162697602423E-01] +
                2 * [0.2428216520333660E-01] +
                2 * [0.2245726582681610E-01] +
                2 * [0.2035775505847216E-01] +
                2 * [0.1803221639039129E-01] +
                2 * [0.1553677555584398E-01] +
                2 * [0.1293483966360737E-01] +
                2 * [0.1029711695795636E-01] +
                2 * [0.7703375233279742E-02] +
                2 * [0.5249123454808859E-02] +
                2 * [0.3057753410175531E-02] +
                2 * [0.1289524082610417E-02] +
                2 * [0.1807395644453884E-03] +
                #
                2 * [0.2817631903301660E-01] +
                2 * [0.2807645579381725E-01] +
                2 * [0.2787725147661370E-01] +
                2 * [0.2757974956648187E-01] +
                2 * [0.2718551322962479E-01] +
                2 * [0.2669662292745036E-01] +
                2 * [0.2611567337670610E-01] +
                2 * [0.2544576996546477E-01] +
                2 * [0.2469052474448768E-01] +
                2 * [0.2385405210603854E-01] +
                2 * [0.2294096422938775E-01] +
                2 * [0.2195636630531782E-01] +
                2 * [0.2090585144581202E-01] +
                2 * [0.1979549504809750E-01] +
                2 * [0.1863184825613879E-01] +
                2 * [0.1742193015946417E-01] +
                2 * [0.1617321872957772E-01] +
                2 * [0.1489364166481518E-01] +
                2 * [0.1359157100976555E-01] +
                2 * [0.1227583056008277E-01] +
                2 * [0.1095573338783790E-01] +
                2 * [0.9641177729702537E-02] +
                2 * [0.8342838753968157E-02] +
                2 * [0.7072489995433555E-02] +
                2 * [0.5843449875835640E-02] +
                2 * [0.4671050372114322E-02] +
                2 * [0.3572892783517299E-02] +
                2 * [0.2568764943794020E-02] +
                2 * [0.1681142865421470E-02] +
                2 * [0.9383698485423815E-03] +
                2 * [0.3777466463269846E-03] +
                2 * [0.5053609520786252E-04]
                )

        return

    def _get_points(self, index):
        def _pm(a):
            return [a, -a]

        if index == 0:
            return [0.0]

        if index == 1:
            new_points = (
                _pm(0.7745966692414834)
                )
        elif index == 2:
            new_points = (
                _pm(0.4342437493468025) +
                _pm(0.9604912687080203)
                )
        elif index == 3:
            new_points = (
                _pm(0.2233866864289669) +
                _pm(0.6211029467372264) +
                _pm(0.8884592328722570) +
                _pm(0.9938319632127550)
                )
        elif index == 4:
            new_points = (
                _pm(0.1124889431331866) +
                _pm(0.3311353932579768) +
                _pm(0.5313197436443756) +
                _pm(0.7024962064915271) +
                _pm(0.8367259381688688) +
                _pm(0.9296548574297401) +
                _pm(0.9815311495537401) +
                _pm(0.9990981249676676)
                )
        elif index == 5:
            new_points = (
                _pm(0.5634431304659279E-01) +
                _pm(0.1682352515522075) +
                _pm(0.2777498220218243) +
                _pm(0.3833593241987304) +
                _pm(0.4836180269458411) +
                _pm(0.5771957100520458) +
                _pm(0.6629096600247806) +
                _pm(0.7397560443526947) +
                _pm(0.8069405319502176) +
                _pm(0.8639079381936905) +
                _pm(0.9103711569570043) +
                _pm(0.9463428583734029) +
                _pm(0.9721828747485818) +
                _pm(0.9886847575474295) +
                _pm(0.9972062593722220) +
                _pm(0.9998728881203576)
                )
        else:
            assert index == 6
            new_points = (
                _pm(0.2818464894974569E-01) +
                _pm(0.8445404008371088E-01) +
                _pm(0.1404242331525602) +
                _pm(0.1958975027111002) +
                _pm(0.2506787303034832) +
                _pm(0.3045764415567140) +
                _pm(0.3574038378315322) +
                _pm(0.4089798212298887) +
                _pm(0.4591300119898323) +
                _pm(0.5076877575337166) +
                _pm(0.5544951326319325) +
                _pm(0.5994039302422429) +
                _pm(0.6422766425097595) +
                _pm(0.6829874310910792) +
                _pm(0.7214230853700989) +
                _pm(0.7574839663805136) +
                _pm(0.7910849337998483) +
                _pm(0.8221562543649804) +
                _pm(0.8506444947683502) +
                _pm(0.8765134144847053) +
                _pm(0.8997448997769401) +
                _pm(0.9203400254700124) +
                _pm(0.9383203977795929) +
                _pm(0.9537300064257611) +
                _pm(0.9666378515584165) +
                _pm(0.9771415146397057) +
                _pm(0.9853714995985203) +
                _pm(0.9914957211781061) +
                _pm(0.9957241046984072) +
                _pm(0.9983166353184074) +
                _pm(0.9995987996719107) +
                _pm(0.9999824303548916)
                )

        return self._get_points(index-1) + new_points


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
        length = len(N)
        m = n - length
        v0 = numpy.concatenate([
            2.0 / N / (N-2),
            numpy.array([1.0 / N[-1]]),
            numpy.zeros(m),
            ])
        v2 = - v0[:-1] - v0[:0:-1]
        g0 = -numpy.ones(n)
        g0[length] += n
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
        length = len(N)
        m = n - length
        K = numpy.arange(m)

        v0 = numpy.concatenate([
            2 * numpy.exp(1j*numpy.pi*K/n) / (1 - 4*K**2),
            numpy.zeros(length+1)
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
        length = len(N)
        m = n - length
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
        t = sympy.Symbol('t')
        for r in range(n+1):
            # Compare with get_weights().
            f = sympy.prod([(t - i) for i in range(n+1) if i != r])
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
        t = sympy.Symbol('t')
        for r in range(1, n):
            # Compare with get_weights().
            f = sympy.prod([(t - i) for i in range(1, n) if i != r])
            alpha = 2 * \
                (-1)**(n-r+1) * sympy.integrate(f, (t, 0, n)) \
                / (math.factorial(r-1) * math.factorial(n-1-r)) \
                / n
            self.weights[r-1] = alpha
        return


class GaussKronrod(object):
    '''
    Gauss-Kronrod quadrature; see
    <https://en.wikipedia.org/wiki/Gauss%E2%80%93Kronrod_quadrature_formula>.

    Besides points and weights, this class provides the weights of the
    corresponding Gauss-Legendre scheme in self.gauss_weights.

    Code adapted from
    <https://www.cs.purdue.edu/archives/2002/wxg/codes/r_kronrod.m>,
    <https://www.cs.purdue.edu/archives/2002/wxg/codes/kronrod.m>.
    See

    Calculation of Gauss-Kronrod quadrature rules,
    Dirk P. Laurie,
    Math. Comp. 66 (1997), 1133-1145,
    <https://doi.org/10.1090/S0025-5718-97-00861-2>

    Abstract:
    The Jacobi matrix of the $(2n+1)$-point Gauss-Kronrod quadrature rule for a
    given measure is calculated efficiently by a five-term recurrence relation.
    The algorithm uses only rational operations and is therefore also useful
    for obtaining the Jacobi-Kronrod matrix analytically. The nodes and weights
    can then be computed directly by standard software for Gaussian quadrature
    formulas.
    '''
    def __init__(self, n, a=0.0, b=0.0):
        # The general scheme is:
        # Get the Jacobi recursion coefficients, get the Kronrod vectors alpha
        # and beta, and hand those off to _gauss. There, the eigenproblem for a
        # tridiagonal matrix with alpha and beta is solved to retrieve the
        # points and weights.
        length = int(math.ceil(3*n/2.0)) + 1
        self.degree = 2*length + 1
        alpha, beta = _jacobi_recursion_coefficients(length, a, b)
        a, b = self.r_kronrod(n, alpha, beta)
        x, w = _gauss(a, b)
        # sort by x
        i = numpy.argsort(x)
        self.points = x[i]
        self.weights = w[i]
        return

    def r_kronrod(self, n, a0, b0):
        assert len(a0) == int(math.ceil(3*n/2.0)) + 1
        assert len(b0) == int(math.ceil(3*n/2.0)) + 1

        a = numpy.zeros(2*n+1)
        b = numpy.zeros(2*n+1)

        k = int(math.floor(3*n/2.0)) + 1
        a[:k] = a0[:k]
        k = int(math.ceil(3*n/2.0)) + 1
        b[:k] = b0[:k]
        s = numpy.zeros(int(math.floor(n/2.0)) + 2)
        t = numpy.zeros(int(math.floor(n/2.0)) + 2)
        t[1] = b[n+1]
        for m in range(n-1):
            k0 = int(math.floor((m+1)/2.0))
            k = numpy.arange(k0, -1, -1)
            l = m - k
            s[k+1] = numpy.cumsum(
                (a[k+n+1] - a[l])*t[k+1] + b[k+n+1]*s[k] - b[l]*s[k+1]
                )
            s, t = t, s

        j = int(math.floor(n/2.0)) + 1
        s[1:j+1] = s[:j]
        for m in range(n-1, 2*n-2):
            k0 = m+1-n
            k1 = int(math.floor((m-1)/2.0))
            k = numpy.arange(k0, k1 + 1)
            l = m - k
            j = n-1-l
            s[j+1] = numpy.cumsum(
                -(a[k+n+1] - a[l])*t[j+1] - b[k+n+1]*s[j+1] + b[l]*s[j+2]
                )
            j = j[-1]
            k = int(math.floor((m+1)/2.0))
            if m % 2 == 0:
                a[k+n+1] = a[k] + (s[j+1] - b[k+n+1]*s[j+2]) / t[j+2]
            else:
                b[k+n+1] = s[j+1] / s[j+2]
            s, t = t, s

        a[2*n] = a[n-1] - b[2*n] * s[1]/t[1]
        return a, b

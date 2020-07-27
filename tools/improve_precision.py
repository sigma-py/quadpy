import numpy
import orthopy
from scipy.optimize import least_squares, lsq_linear, minimize

import quadpy
from quadpy.enr._helpers import integrate_monomial_over_enr
from quadpy.helpers import untangle


def partition(boxes, balls):
    """Create all nonnegative tuples of length d which sum up to n.
    """
    # <https://stackoverflow.com/a/36748940/353337>
    # See <https://stackoverflow.com/a/45348441/353337> for an alterative
    # solution.
    def rec(boxes, balls, parent=tuple()):
        if boxes > 1:
            for i in range(balls + 1):
                yield from rec(boxes - 1, i, parent + (balls - i,))
        else:
            yield parent + (balls,)

    return list(rec(boxes, balls))


def one():
    def weights_from_points(azimuthal, polar):
        out = orthopy.sphere.tree_sph(
            polar, azimuthal, scheme.degree, standardization="quantum mechanic"
        )

        A = numpy.array([row for level in out for row in level])

        b = numpy.zeros(A.shape[0])
        # b[0] = numpy.sqrt(4 * numpy.pi)
        b[0] = 1.0 / (2 * numpy.sqrt(numpy.pi))

        # solve linear least-squares problem for the weights
        res = lsq_linear(A, b, tol=1.0e-15)
        return res.x, res.fun

    def f(x):
        azimuthal, polar = x.reshape(2, -1)
        _, err = weights_from_points(azimuthal, polar)
        v = numpy.sqrt(err.real ** 2 + err.imag ** 2)
        return v

    scheme = quadpy.sphere.heo_xu_13()
    print(scheme.points)

    # x0 = numpy.column_stack([scheme.weights, scheme.points]).T.reshape(-1)
    x0 = scheme.azimuthal_polar.T.reshape(-1)
    # x0 += 1.0e-10 * numpy.random.rand(*x0.shape)

    out = least_squares(f, x0, gtol=1.0e-16, xtol=1.0e-16)
    # print(out.x)
    # print(out.status, out.nfev, out.njev)
    # print(out.message)
    assert out.success

    azimuthal, polar = out.x.reshape(2, -1)
    w, _ = weights_from_points(azimuthal, polar)
    assert numpy.all(numpy.imag(w) < 1.0e-14)
    w = numpy.real(w)
    print("weights:")
    for item in w:
        print(f"{item:.15e}")

    x = numpy.sin(polar) * numpy.cos(azimuthal)
    y = numpy.sin(polar) * numpy.sin(azimuthal)
    z = numpy.cos(polar)
    X = numpy.column_stack([x, y, z])
    print()
    print("points:")
    for item in X:
        print("{:.15e}  {:.15e}  {:.15e}".format(*item))
    return


def stenger_7a_5():
    from quadpy.helpers import fsd, get_all_exponents, untangle, z
    from quadpy.nball._helpers import integrate_monomial_over_unit_nball

    def f(x):
        degree = 7
        n = 5
        u = x[0]
        v = x[1]
        B = x[2:]
        data = [
            (B[0], z(n)),
            (B[1], fsd(n, (u, 1))),
            (B[2], fsd(n, (v, 1))),
            (B[3], fsd(n, (u, 2))),
            (B[4], fsd(n, (v, 2))),
            (B[5], fsd(n, (u, 3))),
        ]
        points, weights = untangle(data)

        exponents = get_all_exponents(n, degree)
        # flatten list
        exponents = numpy.array([item for sublist in exponents for item in sublist])

        def evaluate_all_monomials(x):
            return numpy.prod(x[..., None] ** exponents.T[:, None], axis=0).T

        flt = numpy.vectorize(float)
        exact_vals = flt([integrate_monomial_over_unit_nball(k) for k in exponents])

        A = evaluate_all_monomials(points.T)

        out = numpy.dot(A, weights)
        out -= exact_vals

        norm_v = numpy.sqrt(numpy.dot(out, out))
        print(norm_v)
        return norm_v

    x0 = [
        0.250_562_808_085_732,
        0.694_746_590_606_866,
        -0.220_221_371_883_822e03,
        +0.730_167_125_339_176e02,
        +0.143_281_369_027_706,
        -0.203_714_128_400_494e02,
        +0.448_293_291_677_155e-01,
        +0.383_685_702_879_441e01,
    ]

    out = minimize(f, x0, method="Powell", tol=1.0e-12, options={"maxiter": 10000})
    print(out.status, out.nfev)
    print(out.message)
    assert out.success
    print()
    for x in out.x:
        print(f"{x:.15e}")
    return


def stenger_11a_5():
    from quadpy.helpers import fsd, get_all_exponents, untangle, z
    from quadpy.nball._helpers import integrate_monomial_over_unit_nball

    def f(x):
        degree = 11
        n = 5
        u = x[0]
        v = x[1]
        w = x[2]
        B = x[3:]
        data = [
            (B[0], z(n)),
            (B[1], fsd(n, (u, 1))),
            (B[2], fsd(n, (v, 1))),
            (B[3], fsd(n, (w, 1))),
            (B[4], fsd(n, (u, 2))),
            (B[5], fsd(n, (v, 2))),
            (B[6], fsd(n, (w, 2))),
            (B[7], fsd(n, (u, 1), (v, 1))),
            (B[8], fsd(n, (u, 1), (w, 1))),
            (B[9], fsd(n, (u, 3))),
            (B[10], fsd(n, (v, 3))),
            (B[11], fsd(n, (w, 3))),
            (B[12], fsd(n, (u, 2), (v, 1))),
        ]
        if n > 3:
            data += [(B[13], fsd(n, (u, 4))), (B[14], fsd(n, (v, 4)))]
        if n > 4:
            data += [(B[15], fsd(n, (u, 5)))]
        points, weights = untangle(data)

        exponents = get_all_exponents(n, degree)
        # flatten list
        exponents = numpy.array([item for sublist in exponents for item in sublist])

        def evaluate_all_monomials(x):
            return numpy.prod(x[..., None] ** exponents.T[:, None], axis=0).T

        flt = numpy.vectorize(float)
        exact_vals = flt([integrate_monomial_over_unit_nball(k) for k in exponents])

        A = evaluate_all_monomials(points.T)

        out = numpy.dot(A, weights)
        out -= exact_vals

        norm_v = numpy.sqrt(numpy.dot(out, out))
        print()
        print(norm_v)
        print()
        for xx in x:
            print(f"{xx:.15e}")
        return norm_v

    x0 = [
        0.819_845_995_463_488,
        0.540_604_637_387_361,
        0.188_677_422_490_785,
        -0.455_346_412_352_218e03,
        -0.643_198_057_179_463,
        -0.313_723_910_937_508,
        +0.142_863_899_851_242e03,
        -0.115_044_196_304_602e-02,
        +0.149_484_688_898_586,
        -0.373_831_314_185_824e02,
        -0.141_469_445_049_282e-01,
        +0.104_863_970_436_266,
        -0.973_070_178_977_534e-03,
        -0.862_234_640_073_899e-02,
        +0.654_476_925_250_512e01,
        +0.153_312_717_360_660e-02,
        -0.611_928_443_128_898e-04,
        +0.622_126_777_947_090e-02,
        +0.887_274_197_302_674e-05,
    ]

    out = minimize(f, x0, method="Powell", tol=1.0e-12, options={"maxiter": 10000})
    print(out.status, out.nfev)
    print(out.message)
    assert out.success
    print()
    for x in out.x:
        print(f"{x:.15e}")
    return


def stenger_11b_3():
    from quadpy.helpers import fsd, get_all_exponents, untangle, z
    from quadpy.nball._helpers import integrate_monomial_over_unit_nball

    def f(x):
        degree = 11
        n = 3

        u = 0.871_740_148_509_601
        v = 0.591_700_181_433_148
        w = 0.209_299_217_902_484

        # u = x[0]
        # v = x[1]
        # w = x[2]
        # B = x[3:]
        B = x

        data = [
            (B[0], z(n)),
            (B[1], fsd(n, (u, 1))),
            (B[2], fsd(n, (v, 1))),
            (B[3], fsd(n, (w, 1))),
            (B[4], fsd(n, (u, 2))),
            (B[5], fsd(n, (v, 2))),
            (B[6], fsd(n, (w, 2))),
            (B[7], fsd(n, (u, 1), (v, 1))),
            (B[8], fsd(n, (u, 1), (w, 1))),
            (B[9], fsd(n, (u, 3))),
            (B[10], fsd(n, (v, 3))),
            (B[11], fsd(n, (w, 3))),
            (B[12], fsd(n, (u, 2), (v, 1))),
        ]

        points, weights = untangle(data)

        exponents = get_all_exponents(n, degree)
        # flatten list
        exponents = numpy.array([item for sublist in exponents for item in sublist])

        def evaluate_all_monomials(x):
            return numpy.prod(x[..., None] ** exponents.T[:, None], axis=0).T

        flt = numpy.vectorize(float)
        exact_vals = flt([integrate_monomial_over_unit_nball(k) for k in exponents])

        A = evaluate_all_monomials(points.T)

        out = numpy.dot(A, weights)
        out -= exact_vals

        norm_v = numpy.sqrt(numpy.dot(out, out))
        # print()
        print(norm_v)
        # print()
        # for xx in x:
        #     print(f"{xx:.15e}")
        return norm_v

    x0 = [
        # 0.871_740_148_509_601,
        # 0.591_700_181_433_148,
        # 0.209_299_217_902_484,
        -0.369_608_422_804_220e02,
        -0.292_281_498_761_429,
        +0.981_554_121_166_593e-01,
        +0.200_417_005_906_264e02,
        +0.915_649_460_852_702e-03,
        +0.140_881_585_655_056,
        -0.109_681_663_185_532e02,
        -0.867_910_432_951_898e-02,
        +0.118_181_679_110_079,
        -0.172_723_379_038_583e-02,
        +0.153_215_767_717_862e-01,
        +0.614_931_311_140_891e01,
        +0.205_381_636_291_801e-02,
    ]

    out = minimize(f, x0, method="Powell", tol=1.0e-15, options={"maxiter": 20000})
    print(out.status, out.nfev)
    print(out.message)
    print()
    for x in out.x:
        print(f"{x:.15e}")
    assert out.success
    return


def stroud_e2r2_gauss():
    from quadpy.helpers import fsd

    def f(x):
        degree = 7
        data = [
            (x[0], fsd(2, (x[3], 1))),
            (x[1], fsd(2, (x[4], 1))),
            (x[2], fsd(2, (x[3], 1), (x[4], 1))),
        ]

        points = numpy.array(
            [
                [0.0, +x[3]],
                [0.0, -x[3]],
                [+x[3], 0.0],
                [-x[3], 0.0],
                #
                [0.0, +x[4]],
                [0.0, -x[4]],
                [+x[4], 0.0],
                [-x[4], 0.0],
                #
                #
                [+x[3], +x[4]],
                [+x[3], -x[4]],
                [-x[3], +x[4]],
                [-x[3], -x[4]],
            ]
        )

        points, weights = untangle(data)

        A = numpy.concatenate(orthopy.e2r2.tree(points.T, degree, symbolic=False))

        out = numpy.dot(A, weights)
        out[0] -= numpy.sqrt(numpy.pi)

        norm_v = numpy.sqrt(numpy.vdot(out, out))
        # print(norm_v)
        return norm_v

    # Cartesian product Gauss formula
    # sqrt6 = numpy.sqrt(6)
    # r, s = [numpy.sqrt((3 + p_m * sqrt6) / 2) for p_m in [+1, -1]]
    # A, B = [(5 - p_m * 2 * sqrt6) / 48 for p_m in [+1, -1]]
    # C = 1.0 / 48.0
    # A *= math.pi
    # B *= math.pi
    # C *= math.pi
    # x0 = [A, B, C, r, s]

    while True:
        x0 = numpy.random.rand(5) * 10 - 5

        print()
        print("x0", x0)

        out = minimize(
            f, x0, method="Nelder-Mead", tol=1.0e-15, options={"maxiter": 20000}
        )
        print(out.status, out.nfev, out.message, "Function value", out.fun)
        # assert out.success
        if abs(out.fun) < 1.0e-10:
            print()
            for x in out.x:
                print(f"{x:.15e}")
            break
    return


def rabinowitz_richter_4():
    from quadpy.e2r._helpers import _s4, _s8, _s40

    def f(x):
        degree = 13
        data = [
            (x[0], [[x[8], x[9]]]),
            (x[1], _s40(x[10])),
            (x[2], _s40(x[11])),
            (x[3], _s4(x[12])),
            (x[4], _s4(x[13])),
            (x[5], _s4(x[14])),
            (x[6], _s8(x[15], x[16])),
            (x[7], _s8(x[17], x[18])),
        ]
        points, weights = untangle(data)

        exponents = numpy.concatenate([partition(2, d) for d in range(degree + 1)])
        exact_vals = numpy.array([integrate_monomial_over_enr(k) for k in exponents])

        def fun(x):
            k = exponents.T
            # <https://stackoverflow.com/a/46689653/353337>
            s = x.shape[1:] + k.shape[1:]
            return (
                (x.reshape(x.shape[0], -1, 1) ** k.reshape(k.shape[0], 1, -1))
                .prod(0)
                .reshape(s)
            )

        A = fun(points.T).T
        print(A)
        print(exact_vals)
        print(sum(weights))
        out = numpy.dot(A, weights) - exact_vals
        nrm = numpy.sqrt(numpy.dot(out, out))
        print(nrm)
        exit(1)
        return nrm

    x0 = [
        +0.349_777_602_241_248_0e1,
        +0.442_580_256_591_559_0e-6,
        +0.455_340_971_239_599_4e-2,
        +0.277_530_326_587_565_2e-4,
        +0.331_277_792_488_418_2e1,
        -0.101_044_092_999_506_7e1,
        +0.112_721_370_308_653_4e-3,
        +0.492_114_301_738_741_9e2,
        #
        0.0,
        0.0,
        19.676_381_860_412_46,
        8.770_037_945_037_203,
        10.205_685_192_384_36,
        3.591_105_603_680_783,
        3.242_171_893_025_389,
        11.941_693_015_408_18,
        4.911_904_665_577_694,
        3.287_383_483_530_638,
        3.162_277_660_168_379,
    ]

    out = minimize(f, x0, method="Nelder-Mead", tol=1.0e-17)
    print(out.status, out.nfev)
    print(out.message)
    assert out.success
    print()
    for x in out.x:
        print(f"{x:.15e}")
    return


def stroud_1967_5():
    from quadpy.helpers import rd

    def f(x):
        degree = 5
        n = 7

        lmbda, xi, mu, gamma = x
        eta = 0
        A = 1 / 9
        B = 1 / 72
        C = B

        # data = [
        #     (B, rd(n, [(+lmbda, 1), (+xi, n - 1)])),
        #     (B, rd(n, [(-lmbda, 1), (-xi, n - 1)])),
        #     (C, rd(n, [(+mu, 2), (+gamma, n - 2)])),
        #     (C, rd(n, [(-mu, 2), (-gamma, n - 2)])),
        #     (2 * A, numpy.full((1, n), eta)),
        # ]
        # points, weights = untangle(data)
        # weights *= numpy.sqrt(numpy.pi) ** n

        data = [
            (B, rd(n, [(+lmbda, 1), (+xi, n - 1)])),
            (B, rd(n, [(-lmbda, 1), (-xi, n - 1)])),
            (C, rd(n, [(+mu, 2), (+gamma, n - 2)])),
            (C, rd(n, [(-mu, 2), (-gamma, n - 2)])),
            (2 * A, numpy.full((1, n), eta)),
        ]

        points, weights = untangle(data)
        weights *= numpy.sqrt(numpy.pi) ** n

        A = numpy.concatenate(orthopy.enr2.tree(points.T, degree, symbolic=False))

        out = numpy.dot(A, weights)
        out[0] -= numpy.sqrt(numpy.sqrt(numpy.pi)) ** n

        norm_v = numpy.sqrt(numpy.vdot(out, out))
        return norm_v

    x0 = [
        2.009_505_637_083_749e00,
        2.774_548_295_173_737e-01,
        -1.062_215_595_206_724e00,
        6.698_352_123_613_097e-01,
        # 2.009_505_6,
        # 0.277_454_83,
        # -1.062_215_60,
        # 0.669_835_21,
    ]

    out = minimize(f, x0, method="Powell", tol=1.0e-20, options={"maxiter": 20000})
    print(out.status, out.nfev, out.message, "Function value", out.fun)
    assert out.success
    print()
    for x in out.x:
        print(f"{x:.15e}")
    return


if __name__ == "__main__":
    stroud_1967_5()

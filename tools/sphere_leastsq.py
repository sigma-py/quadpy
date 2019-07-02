import numpy
import orthopy
import quadpy
from scipy.optimize import least_squares, lsq_linear, minimize


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
        print("{:.15e}".format(item))

    x = numpy.sin(polar) * numpy.cos(azimuthal)
    y = numpy.sin(polar) * numpy.sin(azimuthal)
    z = numpy.cos(polar)
    X = numpy.column_stack([x, y, z])
    print()
    print("points:")
    for item in X:
        print("{:.15e}  {:.15e}  {:.15e}".format(*item))
    return


def two():
    from quadpy.sphere._heo_xu import _f, _f1, _f2
    from quadpy.helpers import untangle
    from quadpy.sphere._helpers import cartesian_to_spherical

    def f(x):
        degree = 13
        data = [
            (x[0], _f((1.0, 1))),
            (x[1], _f2(x[4])),
            (x[2], _f2(x[5])),
            (x[3], _f1(x[6])),
        ]
        points, weights = untangle(data)
        # print(sum(weights))
        azimuthal, polar = cartesian_to_spherical(points).T

        out = orthopy.sphere.tree_sph(
            polar, azimuthal, degree, standardization="quantum mechanic"
        )

        A = numpy.array([row for level in out for row in level])
        out = numpy.dot(A, weights)
        out[0] -= 1.0 / (2 * numpy.sqrt(numpy.pi))
        v = numpy.sqrt(out.real ** 2 + out.imag ** 2)
        # return v
        norm_v = numpy.sqrt(numpy.vdot(v, v))
        print(norm_v)
        return norm_v

    x0 = numpy.array(
        [
            0.013_866_592_105,
            0.013_050_931_863,
            0.013_206_423_223,
            0.011_942_663_555,
            0.286_640_146_767,
            0.659_905_001_656,
            0.539_490_098_706,
        ]
    )
    # out = least_squares(f, x0, jac="3-point", gtol=1.0e-15, xtol=1.0e-50, ftol=1.0e-15)
    out = minimize(f, x0, method="Nelder-Mead", tol=1.0e-17)
    assert out.success
    print(out.status, out.nfev)
    print(out.message)
    print()
    for x in out.x:
        print(f"{x:.15e}")

    return


if __name__ == "__main__":
    two()

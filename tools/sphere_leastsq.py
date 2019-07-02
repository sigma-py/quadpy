import numpy
import orthopy
import quadpy
from scipy.optimize import least_squares, lsq_linear


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

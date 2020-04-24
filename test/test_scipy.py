import numpy

import quadpy


def test_simple():
    val, err = quadpy.quad(lambda x: numpy.sin(x) - x, 0.0, 6.0)

    ref = -17 - numpy.cos(6)
    assert abs(val - ref) < 1.0e-8 * abs(ref)


def test_ln():
    val, err = quadpy.quad(lambda x: numpy.log(x), 0.5, 5.0)

    ref = (numpy.log(5) * 5 - 5) - (numpy.log(0.5) * 0.5 - 0.5)
    assert abs(val - ref) < 1.0e-8 * abs(ref)


def test_vector():
    val, err = quadpy.quad(lambda x: [numpy.sin(x) - x, x], 0.0, 6.0)

    ref = [-17 - numpy.cos(6), 18]
    assert numpy.all(numpy.abs(val - ref) < 1.0e-8 * numpy.abs(ref))


def test_args():
    val, err = quadpy.quad(lambda x, a: numpy.sin(a * x) - x, 0.0, 6.0, args=(1,))

    ref = -17 - numpy.cos(6)
    assert abs(val - ref) < 1.0e-8 * abs(ref)


def test_gh255():
    # https://github.com/nschloe/quadpy/issues/255
    T = 2 * numpy.pi

    def ex1(t):
        return numpy.where(
            numpy.logical_and((t % T >= 0), (t % T < numpy.pi)), t % T, numpy.pi
        )

    a = 6.3
    val, err = quadpy.quad(ex1, 0.0, a, epsabs=1.0e-13)

    # ref = 14.804547973620716
    ref = 1.5 * numpy.pi ** 2 + (a - 2 * numpy.pi) ** 2 / 2
    assert abs(val - ref) < 1.0e-8 * abs(ref), "\n" + "\n".join(
        [f"reference = {ref}", f"computed  = {val}", f"error     = {abs(val-ref)}"]
    )


if __name__ == "__main__":
    # test_ln()
    test_gh255()

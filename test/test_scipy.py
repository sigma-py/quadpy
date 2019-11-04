import numpy
import quadpy


def test_simple():
    val, err = quadpy.quad(lambda x: numpy.sin(x) - x, 0.0, 6.0)

    ref = -17.960170286650353
    assert abs(val - ref) < 1.0e-10 * abs(ref)


def test_vector():
    val, err = quadpy.quad(lambda x: [numpy.sin(x) - x, x], 0.0, 6.0)

    ref = [-17.960170286650353, 18]
    assert numpy.all(numpy.abs(val - ref) < 1.0e-10 * numpy.abs(ref))


def test_args():
    val, err = quadpy.quad(lambda x, a: numpy.sin(a * x) - x, 0.0, 6.0, args=(1,))

    ref = -17.960170286650353
    assert abs(val - ref) < 1.0e-10 * abs(ref)

import numpy
import quadpy


def test_scipy():
    val, err = quadpy.quad(lambda x: numpy.sin(x) - x, [0.0, 6.0])

    ref = -17.960170286650353
    assert abs(val - ref) < 1.0e-10 * abs(ref)

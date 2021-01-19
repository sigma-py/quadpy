import numpy as np

import quadpy


def test_x():
    val, _ = quadpy.quad(lambda x: x, -1.0, 1.0)
    exact = 0.0
    assert abs(exact - val) < 1.0e-10


def test_sin_x():
    val, err = quadpy.quad(lambda x: np.sin(x) - x, 0.0, 6.0)

    ref = -17 - np.cos(6)
    assert abs(val - ref) < 1.0e-8 * abs(ref)


def test_sin_x_other_way():
    val, err = quadpy.quad(lambda x: np.sin(x) - x, 1.0, 0.0)
    # import scipy.integrate
    # val, err = scipy.integrate.quad(lambda x: np.sin(x) - x, 1.0, 0.0)
    ref = -0.5 + np.cos(1)
    assert abs(val - ref) < 1.0e-8 * abs(ref)


def test_ln():
    val, err = quadpy.quad(lambda x: np.log(x), 0.5, 5.0)

    ref = (np.log(5) * 5 - 5) - (np.log(0.5) * 0.5 - 0.5)
    assert abs(val - ref) < 1.0e-8 * abs(ref)


def test_vector():
    val, err = quadpy.quad(lambda x: [np.sin(x) - x, x], 0.0, 6.0)

    ref = [-17 - np.cos(6), 18]
    assert np.all(np.abs(val - ref) < 1.0e-8 * np.abs(ref))


def test_args():
    val, err = quadpy.quad(lambda x, a: np.sin(a * x) - x, 0.0, 6.0, args=(1,))

    ref = -17 - np.cos(6)
    assert abs(val - ref) < 1.0e-8 * abs(ref)


def test_gh255a():
    # https://github.com/nschloe/quadpy/issues/255
    T = 2 * np.pi

    def ex1(t):
        return np.where(np.logical_and((t % T >= 0), (t % T < np.pi)), t % T, np.pi)

    a = 3.15
    val, err = quadpy.quad(ex1, 0.0, a, epsabs=1.0e-13)

    ref = 0.5 * np.pi ** 2 + (a - np.pi) * np.pi

    assert abs(val - ref) < 1.0e-8 * abs(ref), "\n" + "\n".join(
        [f"reference = {ref}", f"computed  = {val}", f"error     = {abs(val-ref)}"]
    )


def test_gh255b():
    # https://github.com/nschloe/quadpy/issues/255
    T = 2 * np.pi

    def ex1(t):
        return np.where(np.logical_and((t % T >= 0), (t % T < np.pi)), t % T, np.pi)

    a = 6.3
    val, err = quadpy.quad(ex1, 0.0, a, epsabs=1.0e-13, limit=100)

    # ref = 14.804547973620716
    ref = 1.5 * np.pi ** 2 + (a - 2 * np.pi) ** 2 / 2
    assert abs(val - ref) < 1.0e-8 * abs(ref), "\n" + "\n".join(
        [f"reference = {ref}", f"computed  = {val}", f"error     = {abs(val-ref)}"]
    )


def test_gh295():
    def f(x):
        return 1e-20 * np.sin(x)

    # import scipy.integrate
    # out = scipy.integrate.quad(f, 0.0, 1.0, epsabs=0.0, epsrel=1.0e-10)
    quadpy.quad(f, 0.0, 1.0, epsabs=1.0e-8, epsrel=1.0e-8)


def test_complex_valued():
    def f(x):
        return np.exp(1j * x)

    # import scipy.integrate
    # out = scipy.integrate.quad(f, 0.0, 1.0, epsabs=0.0, epsrel=1.0e-10)
    val, _ = quadpy.quad(f, 0.0, 1.0, epsabs=1.0e-8, epsrel=1.0e-8)
    exact = np.sin(1.0) - 1j * (np.cos(1.0) - 1.0)
    assert np.abs(val - exact) < 1.0e-10


if __name__ == "__main__":
    test_x()

<p align="center">
  <a href="https://github.com/nschloe/quadpy"><img alt="quadpy" src="https://raw.githubusercontent.com/sigma-py/quadpy/assets/logo-with-text.svg" width="60%"></a>
  <p align="center">Your one-stop shop for numerical integration in Python.</p>
</p>

[![PyPi Version](https://img.shields.io/pypi/v/quadpy.svg?style=flat-square)](https://pypi.org/project/quadpy/)
[![PyPI pyversions](https://img.shields.io/pypi/pyversions/quadpy.svg?style=flat-square)](https://pypi.org/project/quadpy/)
[![GitHub stars](https://img.shields.io/github/stars/nschloe/quadpy.svg?style=flat-square&logo=github&label=Stars&logoColor=white)](https://github.com/nschloe/quadpy)
[![PyPi downloads](https://img.shields.io/pypi/dm/quadpy.svg?style=flat-square)](https://pypistats.org/packages/quadpy)

[![Discord](https://img.shields.io/static/v1?logo=discord&label=chat&message=on%20discord&color=7289da&style=flat-square)](https://discord.gg/hnTJ5MRX2Y)
[![awesome](https://img.shields.io/badge/awesome-yes-brightgreen.svg?style=flat-square)](https://github.com/nschloe/quadpy)

More than 1500 numerical integration schemes for
[line segments](#line-segment-c1),
[circles](#circle-u2),
[disks](#disk-s2),
[triangles](#triangle-t2),
[quadrilaterals](#quadrilateral-c2),
[spheres](#sphere-u3),
[balls](#ball-s3),
[tetrahedra](#tetrahedron-t3),
[hexahedra](#hexahedra-c3),
[wedges](#wedge-w3),
[pyramids](#pyramid-p3),
[n-spheres](#n-sphere-un),
[n-balls](#n-ball-sn),
[n-cubes](#n-cube-cn),
[n-simplices](#n-simplex-tn),
[the 1D half-space with weight functions exp(-r)](#1d-half-space-with-weight-function-exp-r-e1r),
[the 2D space with weight functions exp(-r)](#2d-space-with-weight-function-exp-r-e2r),
[the 3D space with weight functions exp(-r)](#3d-space-with-weight-function-exp-r-e3r),
[the nD space with weight functions exp(-r)](#nd-space-with-weight-function-exp-r-enr),
[the 1D space with weight functions exp(-r<sup>2</sup>)](#1d-space-with-weight-function-exp-r2-e1r2),
[the 2D space with weight functions exp(-r<sup>2</sup>)](#2d-space-with-weight-function-exp-r2-e2r2),
[the 3D space with weight functions exp(-r<sup>2</sup>)](#3d-space-with-weight-function-exp-r2-e3r2),
and
[the nD space with weight functions exp(-r<sup>2</sup>)](#nd-space-with-weight-function-exp-r2-enr2),
for fast integration of real-, complex-, and vector-valued functions.

For example, to numerically integrate any function over any given interval, install
quadpy [from the Python Package Index](https://pypi.org/project/quadpy/) with

```
pip install quadpy
```

and do

```python
import numpy as np
import quadpy


def f(x):
    return np.sin(x) - x


val, err = quadpy.quad(f, 0.0, 6.0)
```

This is like
[scipy](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.quad.html)
with the addition that quadpy handles complex-, vector-, matrix-valued integrands,
and "intervals" in spaces of arbitrary dimension.

To integrate over a _triangle_, do

```python
import numpy as np
import quadpy


def f(x):
    return np.sin(x[0]) * np.sin(x[1])


triangle = np.array([[0.0, 0.0], [1.0, 0.0], [0.7, 0.5]])

# get a "good" scheme of degree 10
scheme = quadpy.t2.get_good_scheme(10)
val = scheme.integrate(f, triangle)
```

Most domains have `get_good_scheme(degree)`. If you would like to use a particular
scheme, you can pick one from the dictionary `quadpy.t2.schemes`.

All schemes have

<!--pytest-codeblocks:skip-->

```python
scheme.points
scheme.weights
scheme.degree
scheme.source
scheme.test_tolerance

scheme.show()
scheme.integrate(
    # ...
)
```

and many have

<!--pytest-codeblocks:skip-->

```python
scheme.points_symbolic
scheme.weights_symbolic
```

quadpy is fully vectorized, so if you like to compute the integral of a function on many
domains at once, you can provide them all in one `integrate()` call, e.g.,

<!--pytest-codeblocks:skip-->

```python
# shape (3, 5, 2), i.e., (corners, num_triangles, xy_coords)
triangles = np.stack(
    [
        [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]],
        [[1.2, 0.6], [1.3, 0.7], [1.4, 0.8]],
        [[26.0, 31.0], [24.0, 27.0], [33.0, 28]],
        [[0.1, 0.3], [0.4, 0.4], [0.7, 0.1]],
        [[8.6, 6.0], [9.4, 5.6], [7.5, 7.4]],
    ],
    axis=-2,
)
```

The same goes for functions with vectorized output, e.g.,

```python
def f(x):
    return [np.sin(x[0]), np.sin(x[1])]
```

More examples under [test/examples_test.py](test/examples_test.py).

Read more about the dimensionality of the input/output arrays [in the
wiki](https://github.com/nschloe/quadpy/wiki#dimensionality-of-input-and-output-arrays).

Advanced topics:

- [Adaptive quadrature](https://github.com/nschloe/quadpy/wiki/Adaptive-quadrature)
- [Creating your own Gauss scheme](https://github.com/nschloe/quadpy/wiki/Creating-your-own-Gauss-quadrature-in-two-simple-steps)
- [tanh-sinh quadrature](https://github.com/nschloe/tanh_sinh/)

## Schemes

### Line segment (_C<sub>1</sub>_)

<img src="https://raw.githubusercontent.com/sigma-py/quadpy/assets/line-segment-gauss-legendre-20.svg" width="50%">

- [Chebyshev-Gauss](src/quadpy/c1/_chebyshev_gauss.py) (type 1 and 2, arbitrary degree)
- [Clenshaw-Curtis](src/quadpy/c1/_clenshaw_curtis.py) (arbitrary degree)
- [Fejér](src/quadpy/c1/_fejer.py) (type 1 and 2, arbitrary degree)
- [Gauss-Jacobi](src/quadpy/c1/_gauss_jacobi.py) (arbitrary degree)
- [Gauss-Legendre](src/quadpy/c1/_gauss_legendre.py) (arbitrary degree)
- [Gauss-Lobatto](src/quadpy/c1/_gauss_lobatto.py) (arbitrary degree)
- [Gauss-Kronrod](src/quadpy/c1/_gauss_kronrod.py) (arbitrary degree)
- [Gauss-Patterson](src/quadpy/c1/_gauss_patterson.py) (9 nested schemes up to degree 767)
- [Gauss-Radau](src/quadpy/c1/_gauss_radau.py) (arbitrary degree)
- [Newton-Cotes](src/quadpy/c1/_newton_cotes.py) (open and closed, arbitrary degree)

[See
here](https://github.com/nschloe/quadpy/wiki/Creating-your-own-Gauss-quadrature-in-two-simple-steps)
for how to generate Gauss formulas for your own weight functions.

Example:

```python
import numpy as np
import quadpy

scheme = quadpy.c1.gauss_patterson(5)
scheme.show()
val = scheme.integrate(lambda x: np.exp(x), [0.0, 1.0])
```

### 1D half-space with weight function exp(-r) (_E<sub>1</sub><sup>r</sup>_)

<img src="https://raw.githubusercontent.com/sigma-py/quadpy/assets/e1r-gauss-laguerre-3.svg" width="50%">

- [Generalized Gauss-Laguerre](src/quadpy/e1r/_gauss_laguerre.py)

Example:

```python
import quadpy

scheme = quadpy.e1r.gauss_laguerre(5, alpha=0)
scheme.show()
val = scheme.integrate(lambda x: x**2)
```

### 1D space with weight function exp(-r<sup>2</sup>) (_E<sub>1</sub><sup>r<sup>2</sup></sup>_)

<img src="https://raw.githubusercontent.com/sigma-py/quadpy/assets/e1r2-gauss-hermite-8.svg" width="50%">

- [Gauss-Hermite](src/quadpy/e1r2/_gauss_hermite.py) (arbitrary degree)
- [Genz-Keister](src/quadpy/e1r2/_genz_keister.py) (1996, 8 nested schemes up to degree 67)

Example:

```python
import quadpy

scheme = quadpy.e1r2.gauss_hermite(5)
scheme.show()
val = scheme.integrate(lambda x: x**2)
```

### Circle (_U<sub>2</sub>_)

<img src="https://raw.githubusercontent.com/sigma-py/quadpy/assets/circle-krylov-30.svg" width="25%">

- [Krylov](src/quadpy/u2/_krylov.py) (1959, arbitrary degree)

Example:

```python
import numpy as np
import quadpy

scheme = quadpy.u2.get_good_scheme(7)
scheme.show()
val = scheme.integrate(lambda x: np.exp(x[0]), [0.0, 0.0], 1.0)
```

### Triangle (_T<sub>2</sub>_)

<img src="https://raw.githubusercontent.com/sigma-py/quadpy/assets/triangle-dunavant-15.svg" width="25%">

Apart from the classical centroid, vertex, and seven-point schemes we have

- [Hammer-Marlowe-Stroud](src/quadpy/t2/_hammer_marlowe_stroud.py) (1956, 5 schemes up to
  degree 5)
- [Albrecht-Collatz](src/quadpy/t2/_albrecht_collatz.py) (1958, degree 3)
- [Stroud](src/quadpy/t2/_stroud.py) (1971, conical product scheme of degree 7)
- [Franke](src/quadpy/t2/_franke.py) (1971, 2 schemes of degree 7)
- [Strang-Fix/Cowper](src/quadpy/t2/_strang_fix_cowper) (1973, 10 schemes up to degree
  7),
- [Lyness-Jespersen](src/quadpy/t2/_lyness_jespersen.py) (1975, 21 schemes up to degree 11,
  two of which are used in [TRIEX](https://doi.org/10.1145/356068.356070)),
- [Lether](src/quadpy/t2/_lether.py) (1976, degree 2n-2, arbitrary n, not symmetric),
- [Hillion](src/quadpy/t2/_hillion.py) (1977, 10 schemes up to degree 3),
- [Laursen-Gellert](src/quadpy/t2/_laursen_gellert) (1978, 17 schemes up to degree 10),
- [CUBTRI](src/quadpy/t2/_cubtri) (Laurie, 1982, degree 8),
- [Dunavant](src/quadpy/t2/_dunavant) (1985, 20 schemes up to degree 20),
- [Cools-Haegemans](src/quadpy/t2/_cools_haegemans) (1987, degrees 8 and 11),
- [Gatermann](src/quadpy/t2/_gatermann) (1988, degree 7)
- [Berntsen-Espelid](src/quadpy/t2/_berntsen_espelid) (1990, 4 schemes of degree 13, the
  first one being [DCUTRI](https://doi.org/10.1145/131766.131772)),
- [Liu-Vinokur](src/quadpy/t2/_liu_vinokur.py) (1998, 13 schemes up to degree 5),
- [Griener-Schmid](src/quadpy/t2/_griener_schmid), (1999, 2 schemes of degree 6),
- [Walkington](src/quadpy/t2/_walkington.py) (2000, 5 schemes up to degree 5),
- [Wandzura-Xiao](src/quadpy/t2/_wandzura_xiao) (2003, 6 schemes up to degree 30),
- [Taylor-Wingate-Bos](src/quadpy/t2/_taylor_wingate_bos) (2005, 5 schemes up to degree
  14),
- [Zhang-Cui-Liu](src/quadpy/t2/_zhang_cui_liu) (2009, 3 schemes up to degree 20),
- [Xiao-Gimbutas](src/quadpy/t2/_xiao_gimbutas) (2010, 50 schemes up to degree 50),
- [Vioreanu-Rokhlin](src/quadpy/t2/_vioreanu_rokhlin) (2014, 20 schemes up to degree 62),
- [Williams-Shunn-Jameson](src/quadpy/t2/_williams_shunn_jameson) (2014, 8 schemes up to
  degree 12),
- [Witherden-Vincent](src/quadpy/t2/_witherden_vincent) (2015, 19 schemes up to degree 20),
- [Papanicolopulos](src/quadpy/t2/_papanicolopulos) (2016, 27 schemes up to degree 25),
- [all schemes for the n-simplex](#n-simplex-tn).

Example:

```python
import numpy as np
import quadpy

scheme = quadpy.t2.get_good_scheme(12)
scheme.show()
val = scheme.integrate(lambda x: np.exp(x[0]), [[0.0, 0.0], [1.0, 0.0], [0.5, 0.7]])
```

### Disk (_S<sub>2</sub>_)

<img src="https://raw.githubusercontent.com/sigma-py/quadpy/assets/disk-hammer-stroud-20.svg" width="25%">

- [Radon](src/quadpy/s2/_radon.py) (1948, degree 5)
- [Peirce](src/quadpy/s2/_peirce_1956.py) (1956, 3 schemes up to degree 11)
- [Peirce](src/quadpy/s2/_peirce_1957.py) (1957, arbitrary degree)
- [Albrecht-Collatz](src/quadpy/s2/_albrecht_collatz.py) (1958, degree 3)
- [Hammer-Stroud](src/quadpy/s2/_hammer_stroud.py) (1958, 8 schemes up to degree 15)
- [Albrecht](src/quadpy/s2/_albrecht.py) (1960, 8 schemes up to degree 17)
- [Mysovskih](src/quadpy/s2/_mysovskih.py) (1964, 3 schemes up to degree 15)
- [Rabinowitz-Richter](src/quadpy/s2/_rabinowitz_richter) (1969, 6 schemes up to degree 15)
- [Lether](src/quadpy/s2/_lether.py) (1971, arbitrary degree)
- [Piessens-Haegemans](src/quadpy/s2/_piessens_haegemans/) (1975, 1 scheme of degree 9)
- [Haegemans-Piessens](src/quadpy/s2/_haegemans_piessens/) (1977, degree 9)
- [Cools-Haegemans](src/quadpy/s2/_cools_haegemans/) (1985, 4 schemes up to degree 13)
- [Wissmann-Becker](src/quadpy/s2/_wissmann_becker.py) (1986, 3 schemes up to degree 8)
- [Kim-Song](src/quadpy/s2/_kim_song/) (1997, 15 schemes up to degree 17)
- [Cools-Kim](src/quadpy/s2/_cools_kim/) (2000, 3 schemes up to degree 21)
- [Luo-Meng](src/quadpy/s2/_luo_meng/) (2007, 6 schemes up to degree 17)
- [all schemes from the n-ball](#n-ball-sn)

Example:

```python
import numpy as np
import quadpy

scheme = quadpy.s2.get_good_scheme(6)
scheme.show()
val = scheme.integrate(lambda x: np.exp(x[0]), [0.0, 0.0], 1.0)
```

### Quadrilateral (_C<sub>2</sub>_)

<img src="https://raw.githubusercontent.com/sigma-py/quadpy/assets/quad-maxwell.svg" width="25%">

- [Maxwell](src/quadpy/c2/_maxwell.py) (1890, degree 7)
- [Burnside](src/quadpy/c2/_burnside.py) (1908, degree 5)
- [Irwin](src/quadpy/c2/_irwin.py) (1923, 3 schemes up to degree 5)
- [Tyler](src/quadpy/c2/_tyler.py) (1953, 3 schemes up to degree 7)
- [Hammer-Stroud](src/quadpy/c2/_hammer_stroud.py) (1958, 3 schemes up to degree 7)
- [Albrecht-Collatz](src/quadpy/c2/_albrecht_collatz.py) (1958, 4 schemes up to degree 5)
- [Miller](src/quadpy/c2/_miller.py) (1960, degree 1, degree 11 for harmonic integrands)
- [Meister](src/quadpy/c2/_meister.py) (1966, degree 7)
- [Phillips](src/quadpy/c2/_phillips.py) (1967, degree 7)
- [Rabinowitz-Richter](src/quadpy/c2/_rabinowitz_richter) (1969, 6 schemes up to degree 15)
- [Franke](src/quadpy/c2/_franke.py) (1971, 10 schemes up to degree 9)
- [Piessens-Haegemans](src/quadpy/c2/_piessens_haegemans) (1975, 2 schemes of degree 9)
- [Haegemans-Piessens](src/quadpy/c2/_haegemans_piessens) (1977, degree 7)
- [Schmid](src/quadpy/c2/_schmid) (1978, 3 schemes up to degree 6)
- [Cools-Haegemans](src/quadpy/c2/_cools_haegemans_1985/) (1985, 6 schemes up to degree 17)
- [Dunavant](src/quadpy/c2/_dunavant) (1985, 11 schemes up to degree 19)
- [Morrow-Patterson](src/quadpy/c2/_morrow_patterson) (1985, 2 schemes up to degree 20, single precision)
- [Cohen-Gismalla](src/quadpy/c2/_cohen_gismalla.py), (1986, 2 schemes up to degree 3)
- [Wissmann-Becker](src/quadpy/c2/_wissmann_becker) (1986, 6 schemes up to degree 8)
- [Cools-Haegemans](src/quadpy/c2/_cools_haegemans_1988) (1988, 2 schemes up to degree 13)
- [Waldron](src/quadpy/c2/_waldron.py) (1994, infinitely many schemes of degree 3)
- [Sommariva](src/quadpy/c2/_sommariva) (2012, 55 schemes up to degree 55)
- [Witherden-Vincent](src/quadpy/c2/_witherden_vincent) (2015, 11 schemes up to degree 21)
- products of line segment schemes
- [all schemes from the n-cube](#n-cube-cn)

Example:

```python
import numpy as np
import quadpy

scheme = quadpy.c2.get_good_scheme(7)
val = scheme.integrate(
    lambda x: np.exp(x[0]),
    [[[0.0, 0.0], [1.0, 0.0]], [[0.0, 1.0], [1.0, 1.0]]],
)
```

The points are specified in an array of shape (2, 2, ...) such that `arr[0][0]`
is the lower left corner, `arr[1][1]` the upper right. If your c2
has its sides aligned with the coordinate axes, you can use the convenience
function

<!--pytest-codeblocks:skip-->

```python
quadpy.c2.rectangle_points([x0, x1], [y0, y1])
```

to generate the array.

### 2D space with weight function exp(-r) (_E<sub>2</sub><sup>r</sup>_)

<img src="https://raw.githubusercontent.com/sigma-py/quadpy/assets/e2r-rabinowitz-richter-5.svg" width="25%">

- [Stroud-Secrest](src/quadpy/e2r/_stroud_secrest.py) (1963, 2 schemes up to degree 7)
- [Rabinowitz-Richter](src/quadpy/e2r/_rabinowitz_richter) (1969, 4 schemes up to degree 15)
- [Stroud](src/quadpy/e2r/_stroud.py) (1971, degree 4)
- [Haegemans-Piessens](src/quadpy/e2r/_haegemans_piessens/) (1977, 2 schemes up to degree 9)
- [Cools-Haegemans](src/quadpy/e2r/_cools_haegemans/) (1985, 3 schemes up to degree 13)
- [all schemes from the nD space with weight function exp(-r)](#nd-space-with-weight-function-exp-r-enr)

Example:

```python
import quadpy

scheme = quadpy.e2r.get_good_scheme(5)
scheme.show()
val = scheme.integrate(lambda x: x[0] ** 2)
```

### 2D space with weight function exp(-r<sup>2</sup>) (_E<sub>2</sub><sup>r<sup>2</sup></sup>_)

<img src="https://raw.githubusercontent.com/sigma-py/quadpy/assets/e2r2-rabinowitz-richter-3.svg" width="25%">

- [Stroud-Secrest](src/quadpy/e2r2/_stroud_secrest.py) (1963, 2 schemes up to degree 7)
- [Rabinowitz-Richter](src/quadpy/e2r2/_rabinowitz_richter/) (1969, 5 schemes up to degree 15)
- [Stroud](src/quadpy/e2r2/_stroud.py) (1971, 3 schemes up to degree 7)
- [Haegemans-Piessens](src/quadpy/e2r2/_haegemans_piessens/) (1977, 2 schemes of degree 9)
- [Cools-Haegemans](src/quadpy/e2r2/_cools_haegemans/) (1985, 3 schemes up to degree 13)
- [all schemes from the nD space with weight function exp(-r<sup>2</sup>)](#nd-space-with-weight-function-exp-r2-enr2)

Example:

```python
import quadpy

scheme = quadpy.e2r2.get_good_scheme(3)
scheme.show()
val = scheme.integrate(lambda x: x[0] ** 2)
```

### Sphere (_U<sub>3</sub>_)

<img src="https://raw.githubusercontent.com/sigma-py/quadpy/assets/sphere.png" width="25%">

- [Albrecht-Collatz](src/quadpy/u3/_albrecht_collatz.py) (1958, 5 schemes up to degree 7)
- [McLaren](src/quadpy/u3/_mclaren.py) (1963, 10 schemes up to degree 14)
- [Lebedev](src/quadpy/u3/_lebedev/) (1976, 34 schemes up to degree 131)
- [Bažant-Oh](src/quadpy/u3/_bazant_oh/) (1986, 3 schemes up to degree 13)
- [Heo-Xu](src/quadpy/u3/_heo_xu/) (2001, 27 schemes up to degree 39)
- [Fliege-Maier](src/quadpy/u3/_fliege_maier/) (2007, 4 schemes up to degree 4,
  single-precision)
- [all schemes from the n-sphere](#n-sphere-un)

Example:

```python
import numpy as np
import quadpy

scheme = quadpy.u3.get_good_scheme(19)
# scheme.show()
val = scheme.integrate(lambda x: np.exp(x[0]), [0.0, 0.0, 0.0], 1.0)
```

Integration on the sphere can also be done for functions defined in spherical
coordinates:

```python
import numpy as np
import quadpy


def f(theta_phi):
    theta, phi = theta_phi
    return np.sin(phi) ** 2 * np.sin(theta)


scheme = quadpy.u3.get_good_scheme(19)
val = scheme.integrate_spherical(f)
```

### Ball (_S<sub>3</sub>_)

<img src="https://raw.githubusercontent.com/sigma-py/quadpy/assets/ball.png" width="25%">

- [Ditkin](src/quadpy/s3/_ditkin.py) (1948, 3 schemes up to degree 7)
- [Hammer-Stroud](src/quadpy/s3/_hammer_stroud.py) (1958, 6 schemes up to degree 7)
- [Mysovskih](src/quadpy/s3/_mysovskih.py) (1964, degree 7)
- [Stroud](src/quadpy/s3/_stroud.py) (1971, 2 schemes up to degree 14)
- [all schemes from the n-ball](#n-ball-sn)

Example:

```python
import numpy as np
import quadpy

scheme = quadpy.s3.get_good_scheme(4)
# scheme.show()
val = scheme.integrate(lambda x: np.exp(x[0]), [0.0, 0.0, 0.0], 1.0)
```

### Tetrahedron (_T<sub>3</sub>_)

<img src="https://raw.githubusercontent.com/sigma-py/quadpy/assets/tet.png" width="25%">

- [Hammer-Marlowe-Stroud](src/quadpy/t3/_hammer_marlowe_stroud.py)
  (1956, 3 schemes up to degree 3, also appearing in [Hammer-Stroud](https://doi.org/10.1090/S0025-5718-1958-0102176-6))
- [Stroud](src/quadpy/t3/_stroud.py) (1971, degree 7)
- [Yu](src/quadpy/t3/_yu) (1984, 5 schemes up to degree 6)
- [Keast](src/quadpy/t3/_keast) (1986, 10 schemes up to degree 8)
- [Beckers-Haegemans](src/quadpy/t3/_beckers_haegemans) (1990, degrees 8 and 9)
- [Gatermann](src/quadpy/t3/_gatermann) (1992, degree 5)
- [Liu-Vinokur](src/quadpy/t3/_liu_vinokur.py) (1998, 14 schemes up to degree 5)
- [Walkington](src/quadpy/t3/_walkington/) (2000, 6 schemes up to degree 7)
- [Zhang-Cui-Liu](src/quadpy/t3/_zhang_cui_liu/) (2009, 2 schemes up to degree 14)
- [Xiao-Gimbutas](src/quadpy/t3/_xiao_gimbutas/) (2010, 15 schemes up to degree 15)
- [Shunn-Ham](src/quadpy/t3/_shunn_ham/) (2012, 6 schemes up to degree 7)
- [Vioreanu-Rokhlin](src/quadpy/t3/_vioreanu_rokhlin/) (2014, 10 schemes up to degree 13)
- [Williams-Shunn-Jameson](src/quadpy/t3/_williams_shunn_jameson/) (2014, 1 scheme with
  degree 9)
- [Witherden-Vincent](src/quadpy/t3/_witherden_vincent/) (2015, 9 schemes up to degree 10)
- [Jaśkowiec-Sukumar](src/quadpy/t3/_jaskowiec_sukumar/) (2020, 21 schemes up to degree 20)
- [all schemes for the n-simplex](#n-simplex-tn).

Example:

```python
import numpy as np
import quadpy

scheme = quadpy.t3.get_good_scheme(5)
# scheme.show()
val = scheme.integrate(
    lambda x: np.exp(x[0]),
    [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 0.7, 0.0], [0.3, 0.9, 1.0]],
)
```

### Hexahedron (_C<sub>3</sub>_)

<img src="https://raw.githubusercontent.com/sigma-py/quadpy/assets/hexa.png" width="25%">

- [Sadowsky](src/quadpy/c3/_sadowsky.py) (1940, degree 5)
- [Tyler](src/quadpy/c3/_tyler.py) (1953, 2 schemes up to degree 5)
- [Hammer-Wymore](src/quadpy/c3/_hammer_wymore.py) (1957, degree 7)
- [Albrecht-Collatz](src/quadpy/c3/_albrecht_collatz.py) (1958, degree 3)
- [Hammer-Stroud](src/quadpy/c3/_hammer_stroud.py) (1958, 6 schemes up to degree 7)
- [Mustard-Lyness-Blatt](src/quadpy/c3/_mustard_lyness_blatt.py) (1963, 6 schemes up to degree 5)
- [Stroud](src/quadpy/c3/_stroud_1967.py) (1967, degree 5)
- [Sarma-Stroud](src/quadpy/c3/_sarma_stroud.py) (1969, degree 7)
- [Witherden-Vincent](src/quadpy/c3/_witherden_vincent/) (2015, 7 schemes up to degree degree 11)
- [all schemes from the n-cube](#n-cube-cn)
- Product schemes derived from line segment schemes

Example:

```python
import numpy as np
import quadpy

scheme = quadpy.c3.product(quadpy.c1.newton_cotes_closed(3))
# scheme.show()
val = scheme.integrate(
    lambda x: np.exp(x[0]),
    quadpy.c3.cube_points([0.0, 1.0], [-0.3, 0.4], [1.0, 2.1]),
)
```

### Pyramid (_P<sub>3</sub>_)

<img src="https://raw.githubusercontent.com/sigma-py/quadpy/assets/pyra.png" width="25%">

- [Felippa](src/quadpy/p3/_felippa.py) (2004, 9 schemes up to degree 5)

Example:

```python
import numpy as np
import quadpy

scheme = quadpy.p3.felippa_5()

val = scheme.integrate(
    lambda x: np.exp(x[0]),
    [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.5, 0.7, 0.0],
        [0.3, 0.9, 0.0],
        [0.0, 0.1, 1.0],
    ],
)
```

### Wedge (_W<sub>3</sub>_)

<img src="https://raw.githubusercontent.com/sigma-py/quadpy/assets/wedge.png" width="15%">

- [Felippa](src/quadpy/w3/_felippa.py) (2004, 6 schemes up to degree 6)
- [Kubatko-Yeager-Maggi](src/quadpy/w3/_kubatko_yeager_maggi.py) (2013, 21 schemes up to
  degree 9)

Example:

```python
import numpy as np
import quadpy

scheme = quadpy.w3.felippa_3()
val = scheme.integrate(
    lambda x: np.exp(x[0]),
    [
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 0.7, 0.0]],
        [[0.0, 0.0, 1.0], [1.0, 0.0, 1.0], [0.5, 0.7, 1.0]],
    ],
)
```

### 3D space with weight function exp(-r) (_E<sub>3</sub><sup>r</sup>_)

<img src="https://raw.githubusercontent.com/sigma-py/quadpy/assets/e3r.png" width="25%">

- [Stroud-Secrest](src/quadpy/e3r/_stroud_secrest.py) (1963, 5 schemes up to degree 7)
- [all schemes from the nD space with weight function
  exp(-r)](#nd-space-with-weight-function-exp-r-enr)

Example:

```python
import quadpy

scheme = quadpy.e3r.get_good_scheme(5)
# scheme.show()
val = scheme.integrate(lambda x: x[0] ** 2)
```

### 3D space with weight function exp(-r<sup>2</sup>) (_E<sub>3</sub><sup>r<sup>2</sup></sup>_)

<img src="https://raw.githubusercontent.com/sigma-py/quadpy/assets/e3r2.png" width="25%">

- [Stroud-Secrest](src/quadpy/e3r/_stroud_secrest.py) (1963, 7 schemes up to degree 7)
- [Stroud](src/quadpy/e3r/_stroud.py) (1971, scheme of degree 14)
- [all schemes from the nD space with weight function
  exp(-r<sup>2</sup>)](#nd-space-with-weight-function-exp-r2-enr2)

Example:

```python
import quadpy

scheme = quadpy.e3r2.get_good_scheme(6)
# scheme.show()
val = scheme.integrate(lambda x: x[0] ** 2)
```

### n-Simplex (_T<sub>n</sub>_)

- [Lauffer](src/quadpy/tn/_lauffer.py) (1955, 5 schemes up to degree 5)
- [Hammer-Stroud](src/quadpy/tn/_hammer_stroud.py) (1956, 3 schemes up to degree 3)
- [Stroud](src/quadpy/tn/_stroud_1964.py) (1964, degree 3)
- [Stroud](src/quadpy/tn/_stroud_1966.py) (1966, 7 schemes of degree 3)
- [Stroud](src/quadpy/tn/_stroud_1969.py) (1969, degree 5)
- [Silvester](src/quadpy/tn/_silvester.py) (1970, arbitrary degree),
- [Grundmann-Möller](src/quadpy/tn/_grundmann_moeller.py) (1978, arbitrary degree)
- [Walkington](src/quadpy/tn/_walkington.py) (2000, 5 schemes up to degree 7)

Example:

```python
import numpy as np
import quadpy

dim = 4
scheme = quadpy.tn.grundmann_moeller(dim, 3)
val = scheme.integrate(
    lambda x: np.exp(x[0]),
    np.array(
        [
            [0.0, 0.0, 0.0, 0.0],
            [1.0, 2.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 3.0, 1.0, 0.0],
            [0.0, 0.0, 4.0, 1.0],
        ]
    ),
)
```

### n-Sphere (_U<sub>n</sub>_)

- [Stroud](src/quadpy/un/_stroud_1967.py) (1967, degree 7)
- [Stroud](src/quadpy/un/_stroud_1969.py) (1969, 3 <= n <= 16, degree 11)
- [Stroud](src/quadpy/un/_stroud.py) (1971, 6 schemes up to degree 5)
- [Dobrodeev](src/quadpy/un/_dobrodeev_1978.py) (1978, n >= 2, degree 5)
- [Mysovskikh](src/quadpy/un/_mysovskikh.py) (1980, 2 schemes up to degree 5)

Example:

```python
import numpy as np
import quadpy

dim = 4
scheme = quadpy.un.dobrodeev_1978(dim)
val = scheme.integrate(lambda x: np.exp(x[0]), np.zeros(dim), 1.0)
```

### n-Ball (_S<sub>n</sub>_)

- [Stroud](src/quadpy/sn/_stroud_1957.py) (1957, degree 2)
- [Hammer-Stroud](src/quadpy/sn/_hammer_stroud.py) (1958, 2 schemes up to degree 5)
- [Stroud](src/quadpy/sn/_stroud_1966.py) (1966, 4 schemes of degree 5)
- [Stroud](src/quadpy/sn/_stroud_1967_5.py) (1967, 4 <= n <= 7, 2 schemes of degree 5)
- [Stroud](src/quadpy/sn/_stroud_1967_7.py) (1967, n >= 3, 3 schemes of degree 7)
- [Stenger](src/quadpy/sn/_stenger.py) (1967, 6 schemes up to degree 11)
- [McNamee-Stenger](src/quadpy/sn/_mcnamee_stenger.py) (1967, 6 schemes up to degree 9)
- [Dobrodeev](src/quadpy/sn/_dobrodeev_1970.py) (1970, n >= 3, degree 7)
- [Dobrodeev](src/quadpy/sn/_dobrodeev_1978.py) (1978, 2 <= n <= 20, degree 5)
- [Stoyanova](src/quadpy/sn/_stoyanova.py) (1997, n >= 5, degree 7)

Example:

```python
import numpy as np
import quadpy

dim = 4
scheme = quadpy.sn.dobrodeev_1970(dim)
val = scheme.integrate(lambda x: np.exp(x[0]), np.zeros(dim), 1.0)
```

### n-Cube (_C<sub>n</sub>_)

- [Ewing](src/quadpy/cn/_ewing.py) (1941, degree 3)
- [Tyler](src/quadpy/cn/_tyler.py) (1953, degree 3)
- [Stroud](src/quadpy/cn/_stroud_1957.py) (1957, 2 schemes up to degree 3)
- [Hammer-Stroud](src/quadpy/cn/_hammer_stroud.py) (1958, degree 5)
- [Mustard-Lyness-Blatt](src/quadpy/cn/_mustard_lyness_blatt.py) (1963, degree 5)
- [Thacher](src/quadpy/cn/_thacher.py) (1964, degree 2)
- [Stroud](src/quadpy/cn/_stroud_1966.py) (1966, 4 schemes of degree 5)
- [Phillips](src/quadpy/cn/_phillips.py) (1967, degree 7)
- [McNamee-Stenger](src/quadpy/cn/_mcnamee_stenger.py) (1967, 6 schemes up to degree 9)
- [Stroud](src/quadpy/cn/_stroud_1968.py) (1968, degree 5)
- [Dobrodeev](src/quadpy/cn/_dobrodeev_1970.py) (1970, n >= 5, degree 7)
- [Dobrodeev](src/quadpy/cn/_dobrodeev_1978.py) (1978, n >= 2, degree 5)
- [Cools-Haegemans](src/quadpy/cn/_cools_haegemans.py) (1994, 2 schemes up to degree 5)

Example:

```python
import numpy as np
import quadpy

dim = 4
scheme = quadpy.cn.stroud_cn_3_3(dim)
val = scheme.integrate(
    lambda x: np.exp(x[0]),
    quadpy.cn.ncube_points([0.0, 1.0], [0.1, 0.9], [-1.0, 1.0], [-1.0, -0.5]),
)
```

### nD space with weight function exp(-r) (_E<sub>n</sub><sup>r</sup>_)

- [Stroud-Secrest](src/quadpy/enr/_stroud_secrest.py) (1963, 4 schemes up to degree 5)
- [McNamee-Stenger](src/quadpy/enr/_mcnamee_stenger.py) (1967, 6 schemes up to degree 9)
- [Stroud](src/quadpy/enr/_stroud.py) (1971, 2 schemes up to degree 5)

Example:

```python
import quadpy

dim = 4
scheme = quadpy.enr.stroud_enr_5_4(dim)
val = scheme.integrate(lambda x: x[0] ** 2)
```

### nD space with weight function exp(-r<sup>2</sup>) (_E<sub>n</sub><sup>r<sup>2</sup></sup>_)

- [Stroud-Secrest](src/quadpy/enr2/_stroud_secrest.py) (1963, 4 schemes up to degree 5)
- [McNamee-Stenger](src/quadpy/enr2/_mcnamee_stenger.py) (1967, 6 schemes up to degree 9)
- [Stroud](src/quadpy/enr2/_stroud_1967_5.py) (1967, 2 schemes of degree 5)
- [Stroud](src/quadpy/enr2/_stroud_1967_7.py) (1967, 3 schemes of degree 7)
- [Stenger](src/quadpy/enr2/_stenger.py) (1971, 6 schemes up to degree 11, varying dimensionality restrictions)
- [Stroud](src/quadpy/enr2/_stroud.py) (1971, 5 schemes up to degree 5)
- [Phillips](src/quadpy/enr2/_phillips.py) (1980, degree 5)
- [Cools-Haegemans](src/quadpy/enr2/_cools_haegemans.py) (1994, 3 schemes up to degree 7)
- [Lu-Darmofal](src/quadpy/enr2/_lu_darmofal.py) (2004, degree 5)
- [Xiu](src/quadpy/enr2/_xiu.py) (2008, degree 2)

Example:

```python
import quadpy

dim = 4
scheme = quadpy.enr2.stroud_enr2_5_2(dim)
val = scheme.integrate(lambda x: x[0] ** 2)
```

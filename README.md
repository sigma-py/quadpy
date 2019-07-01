# quadpy

Your one-stop shop for numerical integration in Python.

[![CircleCI](https://img.shields.io/circleci/project/github/nschloe/quadpy/master.svg?style=flat-square)](https://circleci.com/gh/nschloe/quadpy/tree/master)
[![codecov](https://img.shields.io/codecov/c/github/nschloe/quadpy.svg?style=flat-square)](https://codecov.io/gh/nschloe/quadpy)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg?style=flat-square)](https://github.com/ambv/black)
[![awesome](https://img.shields.io/badge/awesome-yes-brightgreen.svg?style=flat-square)](https://github.com/nschloe/quadpy)
[![PyPi Version](https://img.shields.io/pypi/v/quadpy.svg?style=flat-square)](https://pypi.org/project/quadpy)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1173132.svg?style=flat-square)](https://doi.org/10.5281/zenodo.1173132)
[![GitHub stars](https://img.shields.io/github/stars/nschloe/quadpy.svg?style=flat-square&logo=github&label=Stars&logoColor=white)](https://github.com/nschloe/quadpy)
[![PyPi downloads](https://img.shields.io/pypi/dd/quadpy.svg?style=flat-square)](https://pypistats.org/packages/quadpy)

<p align="center">
  <img src="https://nschloe.github.io/quadpy/quad.png" width="20%">
</p>

More than 1500 numerical integration schemes for
[line segments](#line-segment),
[circles](#circle),
[disks](#disk),
[triangles](#triangle),
[quadrilaterals](#quadrilateral),
[spheres](#sphere),
[balls](#ball),
[tetrahedra](#tetrahedron),
[hexahedra](#hexahedron),
[wedges](#wedge),
[pyramids](#pyramid),
[n-spheres](#n-sphere),
[n-balls](#n-ball),
[n-cubes](#n-cube),
[n-simplices](#n-simplex), and the
1D/2D/3D/nD spaces with weight functions exp(-r) and exp(-r<sup>2</sup>)
for fast integration of real-, complex-, and vector-valued functions.

To numerically integrate any function over any given triangle, install quadpy [from the
Python Package Index](https://pypi.org/project/quadpy/) with
```
pip3 install quadpy --user
```
and do
```python
import numpy
import quadpy

def f(x):
    return numpy.sin(x[0]) * numpy.sin(x[1])

triangle = numpy.array([[0.0, 0.0], [1.0, 0.0], [0.7, 0.5]])

val = quadpy.triangle.strang_fix_cowper_09().integrate(f, triangle)
```
This uses [Strang's rule](https://bookstore.siam.org/wc08/) of degree 6.

quadpy is fully vectorized, so if you like to compute the integral of a function on many
domains at once, you can provide them all in one `integrate()` call, e.g.,
```python
# shape (3, 5, 2), i.e., (corners, num_triangles, xy_coords)
triangles = numpy.stack([
    [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]],
    [[1.2, 0.6], [1.3, 0.7], [1.4, 0.8]],
    [[26.0, 31.0], [24.0, 27.0], [33.0, 28]],
    [[0.1, 0.3], [0.4, 0.4], [0.7, 0.1]],
    [[8.6, 6.0], [9.4, 5.6], [7.5, 7.4]]
    ], axis=-2)
```
The same goes for functions with vectorized output, e.g.,
```python
def f(x):
    return [numpy.sin(x[0]), numpy.sin(x[1])]
```

More examples under [test/examples_test.py](https://github.com/nschloe/quadpy/blob/master/test/examples_test.py).

Read more about the dimensionality of the input/output arrays [in the
wiki](https://github.com/nschloe/quadpy/wiki#dimensionality-of-input-and-output-arrays).

## Schemes

### Line segment
<img src="https://nschloe.github.io/quadpy/line_segment.png" width="50%">

 * Chebyshev-Gauss (both variants, arbitrary degree)
 * Clenshaw-Curtis (after
   [Waldvogel](https://doi.org/10.1007/s10543-006-0045-4), arbitrary degree)
 * Fejér-type-1 (after
   [Waldvogel](https://doi.org/10.1007/s10543-006-0045-4), arbitrary degree)
 * Fejér-type-2 (after
   [Waldvogel](https://doi.org/10.1007/s10543-006-0045-4), arbitrary degree)
 * Gauss-Jacobi
 * Gauss-Legendre (via
   [NumPy](https://docs.scipy.org/doc/numpy/reference/generated/numpy.polynomial.legendre.leggauss.html), arbitrary degree)
 * Gauss-Lobatto (arbitrary degree)
 * Gauss-Kronrod (after [Laurie](https://doi.org/10.1090/S0025-5718-97-00861-2), arbitrary degree)
 * [Gauss-Patterson](https://doi.org/10.1090/S0025-5718-68-99866-9) (9 nested schemes up to degree 767)
 * Gauss-Radau (arbitrary degree)
 * closed Newton-Cotes (arbitrary degree)
 * open Newton-Cotes (arbitrary degree)
 * [tanh-sinh quadrature](https://en.wikipedia.org/wiki/Tanh-sinh_quadrature)
   (see above)

[See below](#generating-your-own-gauss-quadrature-in-three-simple-steps) for how to
generate Gauss formulas for your own weight functions.

Example:
```python
import numpy
import quadpy

scheme = quadpy.line_segment.gauss_patterson(5)
scheme.show()
val = scheme.integrate(lambda x: numpy.exp(x), [0.0, 1.0])
```

### 1D half-space with weight function exp(-r)
<img src="https://nschloe.github.io/quadpy/e1r.png" width="50%">

 * Generalized Gauss-Laguerre

Example:
```python
import quadpy

scheme = quadpy.e1r.gauss_laguerre(5, alpha=0)
scheme.show()
val = scheme.integrate(lambda x: x**2)
```


### 1D space with weight function exp(-r<sup>2</sup>)
<img src="https://nschloe.github.io/quadpy/e1r2.png" width="50%">

 * [Gauss-Hermite](https://en.wikipedia.org/wiki/Gauss%E2%80%93Hermite_quadrature) (via
   [NumPy](https://docs.scipy.org/doc/numpy/reference/generated/numpy.polynomial.hermite.hermgauss.html), arbitrary degree)
 * [Genz-Keister](https://doi.org/10.1016/0377-0427(95)00232-4) (1996, 8 nested schemes up to degree 67)

Example:
```python
import quadpy

scheme = quadpy.e1r2.gauss_hermite(5)
scheme.show()
val = scheme.integrate(lambda x: x**2)
```

### Circle
<img src="https://nschloe.github.io/quadpy/circle.png" width="25%">

 * [Krylov](https://books.google.de/books/about/Approximate_Calculation_of_Integrals.html?id=ELeRwR27IRIC&redir_esc=y) (1959, arbitrary degree)

Example:
```python
import quadpy

scheme = quadpy.circle.krylov(7)
scheme.show()
val = scheme.integrate(lambda x: numpy.exp(x[0]), [0.0, 0.0], 1.0)
```

### Triangle
<img src="https://nschloe.github.io/quadpy/triangle.png" width="25%">

Apart from the classical centroid, vertex, and seven-point schemes we have

 * [Hammer-Marlowe-Stroud](https://doi.org/10.1090/S0025-5718-1956-0086389-6)
   (1956, 5 schemes up to degree 5, also appearing in [Hammer-Stroud](https://doi.org/10.1090/S0025-5718-1958-0102176-6))
 * open and closed Newton-Cotes schemes (1970, after [Silvester](https://doi.org/10.1090/S0025-5718-1970-0258283-6), arbitrary degree),
 * via [Stroud](https://books.google.de/books/about/Approximate_calculation_of_multiple_inte.html?id=L_tQAAAAMAAJ&redir_esc=y) (1971):
   - [Albrecht-Collatz](https://doi.org/10.1002/zamm.19580380102) (1958, degree 3)
   - conical product scheme (degree 7)
 * [Franke](https://doi.org/10.1090/S0025-5718-1971-0300440-5) (1971, 2 schemes of degree 7)
 * [Strang-Fix](https://bookstore.siam.org/wc08/)/[Cowper](https://doi.org/10.1002/nme.1620070316) (1973, 10 schemes up to
   degree 7),
 * [Lyness-Jespersen](https://doi.org/10.1093/imamat/15.1.19) (1975, 21
   schemes up to degree 11, two of which are used in [TRIEX](https://doi.org/10.1145/356068.356070)),
 * [Lether](https://doi.org/10.1016/0771-050X(76)90008-5) (1976, degree 2n-2, arbitrary
   n, not symmetric; reproduced in
   [Rathod-Nagaraja-Venkatesudu](https://doi.org/10.1016/j.amc.2006.10.041), 2007),
 * [Hillion](https://doi.org/10.1002/nme.1620110504) (1977, 10 schemes up to
   degree 3),
 * [Grundmann-Möller](https://doi.org/10.1137/0715019) (1978, arbitrary degree),
 * [Laursen-Gellert](https://doi.org/10.1002/nme.1620120107) (1978, 17
   schemes up to degree 10),
 * [CUBTRI](https://dl.acm.org/citation.cfm?id=356001) (Laurie, 1982, degree 8),
 * [Dunavant](https://doi.org/10.1002/nme.1620210612) (1985, 20 schemes up
   to degree 20),
 * [Cools-Haegemans](https://lirias.kuleuven.be/handle/123456789/131869) (1987,
   degrees 8 and 11),
 * [Gatermann](https://doi.org/10.1007/BF02251251) (1988, degree 7)
 * Berntsen-Espelid (1990, 4 schemes of degree 13, the first one being
   [DCUTRI](https://dl.acm.org/citation.cfm?id=131772)),
 * [Liu-Vinokur](https://doi.org/10.1006/jcph.1998.5884) (1998, 13 schemes
   up to degree 5),
 * [Griener-Schmid](https://doi.org/10.1016/S0377-0427(99)00215-0), (1999, 2 schemes of degree 6),
 * [Walkington](https://www.math.cmu.edu/~nw0z/publications/00-CNA-023/023abs/)
   (2000, 5 schemes up to degree 5),
 * [Wandzura-Xiao](https://doi.org/10.1016/S0898-1221(03)90004-6) (2003, 6
   schemes up to degree 30),
 * [Taylor-Wingate-Bos](https://arxiv.org/abs/math/0501496) (2005, 5 schemes up
   to degree 14),
 * [Zhang-Cui-Liu](https://www.jstor.org/stable/43693493) (2009, 3 schemes up to
   degree 20),
 * [Xiao-Gimbutas](https://doi.org/10.1016/j.camwa.2009.10.027) (2010, 50
   schemes up to degree 50),
 * [Vioreanu-Rokhlin](https://doi.org/10.1137/110860082) (2014, 20
   schemes up to degree 62),
 * [Williams-Shunn-Jameson](https://doi.org/10.1016/j.cam.2014.01.007) (2014, 8
   schemes up to degree 12),
 * [Witherden-Vincent](https://doi.org/10.1016/j.camwa.2015.03.017) (2015, 19
   schemes up to degree 20),
 * [Papanicolopulos](https://doi.org/10.1016/j.cam.2015.08.001) (2016, 27
   schemes up to degree 25).

Example:
```python
import quadpy

scheme = quadpy.triangle.xiao_gimbutas_05()
scheme.show()
val = scheme.integrate(lambda x: numpy.exp(x[0]), [[0.0, 0.0], [1.0, 0.0], [0.5, 0.7]])
```

### Disk
<img src="https://nschloe.github.io/quadpy/disk.png" width="25%">

 * [Peirce](https://www.jstor.org/stable/2098722) (1957, arbitrary degree)
 * via [Stroud](https://books.google.de/books/about/Approximate_calculation_of_multiple_inte.html?id=L_tQAAAAMAAJ&redir_esc=y):
   - [Radon](https://eudml.org/doc/176796) (1948, degree 5)
   - [Peirce](https://books.google.de/books/about/Numerical_integration_over_planar_region.html?id=WR9SAAAAMAAJ&redir_esc=y)
     (1956, 3 schemes up to degree 11)
   - [Albrecht-Collatz](https://doi.org/10.1002/zamm.19580380102) (1958, degree 3)
   - [Hammer-Stroud](https://doi.org/10.1090/S0025-5718-1958-0102176-6) (1958, 8 schemes up to degree 15)
   - [Albrecht](https://doi.org/10.1002/zamm.19600401014) (1960, 8 schemes up to degree 17)
   - Mysovskih (1964, 3 schemes up to degree 15)
   - [Rabinowitz-Richter](https://doi.org/10.2307/2004962) (1969, 6 schemes up to degree 15)
 * [Lether](https://www.jstor.org/stable/2949473) (1971, arbitrary degree)
 * [Piessens-Haegemans](https://doi.org/10.2307/2005291) (1975, 1 scheme of degree 9)
 * [Haegemans-Piessens](https://www.jstor.org/stable/2156699) (1977, degree 9)
 * [Cools-Haegemans](https://lirias.kuleuven.be/handle/123456789/131870) (1985, 3 schemes up to degree 9)
 * [Wissmann-Becker](https://doi.org/10.1137/0723043) (1986, 3 schemes up to degree 8)
 * [Cools-Kim](https://link.springer.com/article/10.1007/BF03012263) (2000, 3 schemes up to degree 21)

Example:
```python
import numpy
import quadpy

scheme = quadpy.disk.lether(6)
scheme.show()
val = scheme.integrate(lambda x: numpy.exp(x[0]), [0.0, 0.0], 1.0)
```

### Quadrilateral
<img src="https://nschloe.github.io/quadpy/quad.png" width="25%">

 * [Hammer-Stroud](https://doi.org/10.1090/S0025-5718-1958-0102176-6) (1958, 3
   schemes up to degree 7)
 * via [Stroud](https://books.google.de/books/about/Approximate_calculation_of_multiple_inte.html?id=L_tQAAAAMAAJ&redir_esc=y) (1971, 15 schemes up to degree 15):
   - [Maxwell](https://doi.org/10.1017/CBO9780511710377.061) (1890, degree 7)
   - Burnside (1908, degree 5)
   - [Irwin](https://books.google.de/books/about/On_quadrature_and_cubature.html?id=SuruAAAAMAAJ&redir_esc=y) (1923, 3 schemes up to degree 5)
   - [Tyler](https://doi.org/10.4153/CJM-1953-044-1) (1953, 3 schemes up to degree 7)
   - [Albrecht-Collatz](https://doi.org/10.1002/zamm.19580380102) (1958, 4 schemes up to degree 5)
   - [Miller](https://doi.org/10.2307/2003163) (1960, degree 1)
   - [Meister](https://doi.org/10.1093/comjnl/8.4.368) (1966, degree 7)
   - [Phillips](https://doi.org/10.1093/comjnl/10.2.202) (1967, degree 7)
   - [Rabinowitz-Richter](https://doi.org/10.2307/2004962) (1969, 6 schemes up to degree 15)
 * [Franke](https://doi.org/10.1090/S0025-5718-1971-0300440-5) (1971, 10 schemes up to degree 9)
 * [Piessens-Haegemans](https://doi.org/10.2307/2005291) (1975, 2 schemes of degree 9)
 * [Haegemans-Piessens](https://www.jstor.org/stable/2156699) (1977, degree 7)
 * [Schmid](https://eudml.org/doc/132580) (1978, 3 schemes up to degree 6)
 * [Cools-Haegemans](https://lirias.kuleuven.be/handle/123456789/131870) (1985, 3 schemes up to degree 13)
 * [Dunavant](https://doi.org/10.1002/nme.1620211004) (1985, 11 schemes up to degree 19)
 * [Morrow-Patterson](https://doi.org/10.1137/0722071) (1985, 2 schemes up to degree 20, single precision)
 * [Cohen-Gismalla](https://doi.org/10.1080/00207168608803504), (1986, 2 schemes up to degree 3, single precision)
 * [Wissmann-Becker](https://doi.org/10.1137/0723043) (1986, 6 schemes up to degree 8)
 * [Cools-Haegemans](https://doi.org/10.1007/BF02247942) (1988, 2 schemes up to degree 13)
 * [Waldron](http://ftp.cs.wisc.edu/Approx/symmetries.pdf) (1994, infinitely many schemes of degree 3)
 * [Witherden-Vincent](https://doi.org/10.1016/j.camwa.2015.03.017) (2015, 11 schemes up to degree 21)
 * [Sommariva](https://www.math.unipd.it/~alvise/POINTSETS/set_amr_square.m) (2012, 55 schemes up to degree 55)
 * products of line segment schemes
 * all formulas from the n-cube

Example:
```python
import numpy
import quadpy

scheme = quadpy.quadrilateral.stroud_c2_7_2()
val = scheme.integrate(
    lambda x: numpy.exp(x[0]),
    [[[0.0, 0.0], [1.0, 0.0]], [[0.0, 1.0], [1.0, 1.0]]],
    )
```
The points are specified in an array of shape (2, 2, ...) such that `arr[0][0]`
is the lower left corner, `arr[1][1]` the upper right. If your quadrilateral
has its sides aligned with the coordinate axes, you can use the convenience
function
```python
quadpy.quadrilateral.rectangle_points([x0, x1], [y0, y1])
```
to generate the array.


### 2D space with weight function exp(-r)
<img src="https://nschloe.github.io/quadpy/e2r.png" width="25%">

 * via [Stroud](https://books.google.de/books/about/Approximate_calculation_of_multiple_inte.html?id=L_tQAAAAMAAJ&redir_esc=y) (1971):
   - [Stroud-Secrest](https://doi.org/10.1090/S0025-5718-1963-0161473-0) (1963, 2 schemes up to degree 7)
   - [Rabinowitz-Richter](https://doi.org/10.2307/2004962) (1969, 4 schemes up to degree 15)
   - a scheme of degree 4
 * [Haegemans-Piessens](https://www.jstor.org/stable/2156699) (1977, 2 schemes up to degree 9)

Example:
```python
import quadpy

scheme = quadpy.e2r.rabinowitz_richter_5()
scheme.show()
val = scheme.integrate(lambda x: x[0]**2)
```


### 2D space with weight function exp(-r<sup>2</sup>)
<img src="https://nschloe.github.io/quadpy/e2r2.png" width="25%">

 * via [Stroud](https://books.google.de/books/about/Approximate_calculation_of_multiple_inte.html?id=L_tQAAAAMAAJ&redir_esc=y) (1971):
   - [Stroud-Secrest](https://doi.org/10.1090/S0025-5718-1963-0161473-0) (1963, 2 schemes up to degree 7)
   - [Rabinowitz-Richter](https://doi.org/10.2307/2004962) (1969, 5 schemes up to degree 15)
   - 3 schemes up to degree 7
 * [Haegemans-Piessens](https://www.jstor.org/stable/2156699) (1977, 2 schemes of degree 9)

Example:
```python
import quadpy

scheme = quadpy.e2r2.rabinowitz_richter_3()
scheme.show()
val = scheme.integrate(lambda x: x[0]**2)
```


### Sphere
<img src="https://nschloe.github.io/quadpy/sphere.png" width="25%">

 * via [Stroud](https://books.google.de/books/about/Approximate_calculation_of_multiple_inte.html?id=L_tQAAAAMAAJ&redir_esc=y) (1971):
   - [Albrecht-Collatz](https://doi.org/10.1002/zamm.19580380102) (1958, 5
     schemes up to degree 7)
   - [McLaren](https://doi.org/10.1090/S0025-5718-1963-0159418-2) (1963, 10 schemes up to degree 14)
 * [Lebedev](https://en.wikipedia.org/wiki/Lebedev_quadrature) (1976, 34
   schemes up to degree 131)
 * [Bažant-Oh](https://doi.org/10.1002/zamm.19860660108) (1986, 3 schemes up to degree 11)
 * [Heo-Xu](https://doi.org/10.1090/S0025-5718-00-01198-4) (2001, 27 schemes up
   to degree 39, single-precision)
 * [Fliege-Maier](https://www.personal.soton.ac.uk/jf1w07/nodes/nodes.html) (2007, 4 schemes up
   to degree 4, single-precision)

Example:
```python
import numpy
import quadpy

scheme = quadpy.sphere.lebedev_019()
scheme.show()
val = scheme.integrate(lambda x: numpy.exp(x[0]), [0.0, 0.0, 0.0], 1.0)
```
Integration on the sphere can also be done for function defined in spherical
coordinates:
```python
import numpy
import quadpy

scheme = quadpy.sphere.lebedev_019()
val = scheme.integrate_spherical(
    lambda azimuthal, polar: numpy.sin(azimuthal)**2 * numpy.sin(polar),
    )
```

### Ball
<img src="https://nschloe.github.io/quadpy/ball.png" width="25%">

 * [Hammer-Stroud](https://doi.org/10.1090/S0025-5718-1958-0102176-6) (1958, 6 schemes up to degree 7)
 * via: [Stroud](https://books.google.de/books/about/Approximate_calculation_of_multiple_inte.html?id=L_tQAAAAMAAJ&redir_esc=y) (1971):
   - Ditkin (1948, 3 schemes up to degree 7)
   - Mysovskih (1964, degree 7)
   - 2 schemes up to degree 14

Example:
```python
import numpy
import quadpy

scheme = quadpy.ball.hammer_stroud_14_3a()
scheme.show()
val = scheme.integrate(
    lambda x: numpy.exp(x[0]),
    [0.0, 0.0, 0.0], 1.0,
    )
```


### Tetrahedron
<img src="https://nschloe.github.io/quadpy/tet.png" width="25%">

 * [Hammer-Marlowe-Stroud](https://doi.org/10.1090/S0025-5718-1956-0086389-6)
   (1956, 3 schemes up to degree 3, also appearing in [Hammer-Stroud](https://doi.org/10.1090/S0025-5718-1958-0102176-6))
 * open and closed Newton-Cotes (1970, after [Silvester](https://doi.org/10.1090/S0025-5718-1970-0258283-6)) (arbitrary degree)
 * [Stroud](https://cds.cern.ch/record/104291?ln=en) (1971, degree 7)
 * [Grundmann-Möller](https://doi.org/10.1137/0715019) (1978, arbitrary degree),
 * [Yu](https://doi.org/10.1016/0045-7825(84)90072-0) (1984, 5 schemes up to degree 6)
 * [Keast](https://doi.org/10.1016/0045-7825(86)90059-9) (1986, 10 schemes up to degree 8)
 * [Beckers-Haegemans](https://lirias.kuleuven.be/handle/123456789/132648) (1990, degrees 8 and 9)
 * [Gatermann](https://doi.org/10.1007/978-94-011-2646-5_2) (1992, degree 5)
 * [Liu-Vinokur](https://doi.org/10.1006/jcph.1998.5884) (1998, 14 schemes up to
   degree 5)
 * [Walkington](https://www.math.cmu.edu/~nw0z/publications/00-CNA-023/023abs/)
   (2000, 6 schemes up to degree 7)
 * [Zhang-Cui-Liu](https://www.jstor.org/stable/43693493) (2009, 2 schemes up to
   degree 14)
 * [Xiao-Gimbutas](https://doi.org/10.1016/j.camwa.2009.10.027) (2010, 15
   schemes up to degree 15)
 * [Shunn-Ham](https://doi.org/10.1016/j.cam.2012.03.032) (2012, 6 schemes up to
   degree 7)
 * [Vioreanu-Rokhlin](https://doi.org/10.1137/110860082) (2014, 10
   schemes up to degree 13)
 * [Williams-Shunn-Jameson](https://doi.org/10.1016/j.cam.2014.01.007) (2014, 1
   scheme with degree 9)
 * [Witherden-Vincent](https://doi.org/10.1016/j.camwa.2015.03.017) (2015, 9
   schemes up to degree 10)

Example:
```python
import numpy
import quadpy

scheme = quadpy.tetrahedron.keast_10()
scheme.show()
val = scheme.integrate(
    lambda x: numpy.exp(x[0]),
    [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 0.7, 0.0], [0.3, 0.9, 1.0]],
    )
```

### Hexahedron
<img src="https://nschloe.github.io/quadpy/hexa.png" width="25%">

 * Product schemes derived from line segment schemes
 * via: [Stroud](https://books.google.de/books/about/Approximate_calculation_of_multiple_inte.html?id=L_tQAAAAMAAJ&redir_esc=y) (1971):
   - [Sadowsky](https://doi.org/10.2307/2303834) (1940, degree 5)
   - [Tyler](https://doi.org/10.4153/CJM-1953-044-1) (1953, 2 schemes up to degree 5)
   - [Hammer-Wymore](https://doi.org/10.1090/S0025-5718-1957-0087220-6) (1957, degree 7)
   - [Albrecht-Collatz](https://doi.org/10.1002/zamm.19580380102) (1958, degree 3)
   - [Hammer-Stroud](https://doi.org/10.1090/S0025-5718-1958-0102176-6) (1958, 6 schemes up to degree 7)
   - [Mustard-Lyness-Blatt](https://doi.org/10.1093/comjnl/6.1.75) (1963, 6 schemes up to degree 5)
   - [Stroud](https://doi.org/10.1007/BF02162160) (1967, degree 5)
   - [Sarma-Stroud](https://doi.org/10.2307/2004963) (1969, degree 7)
 * all formulas from the n-cube

Example:
```python
import numpy
import quadpy

scheme = quadpy.hexahedron.product(quadpy.line_segment.newton_cotes_closed(3))
scheme.show()
val = scheme.integrate(
    lambda x: numpy.exp(x[0]),
    quadpy.hexahedron.cube_points([0.0, 1.0], [-0.3, 0.4], [1.0, 2.1]),
    )
```

### Pyramid
<img src="https://nschloe.github.io/quadpy/pyra.png" width="25%">

 * [Felippa](https://doi.org/10.1108/02644400410554362) (9 schemes up to degree 5)

Example:
```python
import numpy
import quadpy

scheme = quadpy.pyramid.felippa_5()

val = scheme.integrate(
    lambda x: numpy.exp(x[0]),
    [
      [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 0.7, 0.0], [0.3, 0.9, 0.0],
      [0.0, 0.1, 1.0],
    ]
    )
```

### Wedge
<img src="https://nschloe.github.io/quadpy/wedge.png" width="15%">

 * [Felippa](https://doi.org/10.1108/02644400410554362) (6 schemes up to degree 6)
 * [Kubatko-Yeager-Maggi](https://doi.org/10.1016/j.compfluid.2013.01.002) (21 schemes
   up to degree 9)

Example:
```python
import numpy
import quadpy

scheme = quadpy.wedge.felippa_3
val = quadpy.wedge.integrate(
    lambda x: numpy.exp(x[0]),
    [
      [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 0.7, 0.0]],
      [[0.0, 0.0, 1.0], [1.0, 0.0, 1.0], [0.5, 0.7, 1.0]],
    ]
    )
```


### 3D space with weight function exp(-r)
<img src="https://nschloe.github.io/quadpy/e3r.png" width="25%">

 * via [Stroud](https://books.google.de/books/about/Approximate_calculation_of_multiple_inte.html?id=L_tQAAAAMAAJ&redir_esc=y) (1971):
   - [Stroud-Secrest](https://doi.org/10.1090/S0025-5718-1963-0161473-0) (1963, 5 schemes up to degree 7)

Example:
```python
import quadpy

scheme = quadpy.e3r.stroud_secrest_09()
scheme.show()
val = scheme.integrate(lambda x: x[0]**2)
```


### 3D space with weight function exp(-r<sup>2</sup>)
<img src="https://nschloe.github.io/quadpy/e3r2.png" width="25%">

 * via [Stroud](https://books.google.de/books/about/Approximate_calculation_of_multiple_inte.html?id=L_tQAAAAMAAJ&redir_esc=y) (1971):
   - [Stroud-Secrest](https://doi.org/10.1090/S0025-5718-1963-0161473-0) (1963, 7 schemes up to degree 7)
   - scheme of degree 14

Example:
```python
import quadpy

scheme = quadpy.e3r2.stroud_secrest_10a
scheme.show()
val = scheme.integrate(lambda x: x[0]**2)
```

### n-Simplex
 * via [Stroud](https://books.google.de/books/about/Approximate_calculation_of_multiple_inte.html?id=L_tQAAAAMAAJ&redir_esc=y):
   - [Lauffer](https://doi.org/10.1007/BF01900222) (1955, 5 schemes up to degree 5)
   - [Hammer-Stroud](https://doi.org/10.2307/2002484) (1956, 3 schemes up to degree 3)
   - [Stroud](https://doi.org/10.2307/2002945) (1964, degree 3)
   - [Stroud](https://doi.org/10.1007/BF02165227) (1966, 7 schemes of degree 3)
   - [Stroud](https://doi.org/10.1137/0706009) (1969, degree 5)
 * [Grundmann-Möller](https://doi.org/10.1137/0715019) (1978, arbitrary degree)
 * [Walkington](https://www.math.cmu.edu/~nw0z/publications/00-CNA-023/023abs/) (2000, 5 schemes up to degree 7)

Example:
```python
import numpy
import quadpy

scheme = quadpy.nsimplex.GrundmannMoeller(dim, 3)
dim = 4
val = scheme.integrate(
    lambda x: numpy.exp(x[0]),
    numpy.array([
        [0.0, 0.0, 0.0, 0.0],
        [1.0, 2.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0],
        [0.0, 3.0, 1.0, 0.0],
        [0.0, 0.0, 4.0, 1.0],
        ])
    )
```

### n-Sphere
 * via [Stroud](https://books.google.de/books/about/Approximate_calculation_of_multiple_inte.html?id=L_tQAAAAMAAJ&redir_esc=y) (1971):
   - [Stroud](https://doi.org/10.1137/0704004) (1967, degree 7)
   - [Stroud](https://doi.org/10.1137/0706009) (1969, 3 <= n <= 16, degree 11)
   - 6 schemes up to degree 5
 * [Dobrodeev](https://doi.org/10.1016/0041-5553(70)90084-4) (1978, n >= 2, degree 5)

Example:
```python
import numpy
import quadpy

scheme = quadpy.nsphere.dobrodeev_1978(dim)
dim = 4
val = scheme.integrate(lambda x: numpy.exp(x[0]), numpy.zeros(dim), 1.0)
```


### n-Ball
 * [Dobrodeev](https://doi.org/10.1016/0041-5553(70)90084-4) (1970, n >= 3, degree 7)
 * via [Stroud](https://books.google.de/books/about/Approximate_calculation_of_multiple_inte.html?id=L_tQAAAAMAAJ&redir_esc=y) (1971):
   - [Stroud](https://doi.org/10.2307/2001945) (1957, degree 2)
   - [Hammer-Stroud](https://doi.org/10.1090/S0025-5718-1958-0102176-6) (1958, 2 schemes up to degree 5)
   - [Stroud](https://doi.org/10.1090/S0025-5718-1966-0191094-8) (1966, 4 schemes of degree 5)
   - [Stroud](https://doi.org/10.1007/BF02162160) (1967, 4 <= n <= 7, 2 schemes of degree 5)
   - [Stroud](https://doi.org/10.1137/0704004) (1967, n >= 3, 3 schemes of degree 7)
   - [Stenger](https://www.jstor.org/stable/2004361) (1967, 6 schemes up to degree 11)
 * [Dobrodeev](https://doi.org/10.1016/0041-5553(70)90084-4) (1978, 2 <= n <= 20, degree 5)

Example:
```python
import numpy
import quadpy

scheme = quadpy.nball.dobrodeev_1970(dim)
dim = 4
val = scheme.integrate(lambda x: numpy.exp(x[0]), numpy.zeros(dim), 1.0)
```

### n-Cube
 * [Dobrodeev](https://doi.org/10.1016/0041-5553(70)90084-4) (1970, n >= 5, degree 7)
 * via [Stroud](https://books.google.de/books/about/Approximate_calculation_of_multiple_inte.html?id=L_tQAAAAMAAJ&redir_esc=y) (1971):
    - [Ewing](https://doi.org/doi.org/10.2307/2303604) (1941, degree 3)
    - [Tyler](https://doi.org/10.4153/CJM-1953-044-1) (1953, degree 3)
    - [Stroud](https://doi.org/10.2307/2001945) (1957, 2 schemes up to degree 3)
    - [Hammer-Stroud](https://doi.org/10.1090/S0025-5718-1958-0102176-6) (1958, degree 5)
    - [Mustard-Lyness-Blatt](https://doi.org/10.1093/comjnl/6.1.75) (1963, degree 5)
    - [Thacher](https://doi.org/10.1145/363872.363897) (1964, degree 2)
    - [Stroud](https://doi.org/10.1090/S0025-5718-1966-0191094-8) (1966, 4 schemes of degree 5)
    - [Phillips](https://doi.org/10.1093/comjnl/10.3.297) (1967, degree 7)
    - [Stroud](https://doi.org/10.2307/2004655) (1968, degree 5)
 * [Dobrodeev](https://doi.org/10.1016/0041-5553(70)90084-4) (1978, n >= 2, degree 5)

Example:
```python
import numpy
import quadpy

dim = 4
scheme = quadpy.ncube.stroud_cn_3_3(dim)
quadpy.ncube.integrate(
    lambda x: numpy.exp(x[0]),
    quadpy.ncube.ncube_points(
        [0.0, 1.0], [0.1, 0.9], [-1.0, 1.0], [-1.0, -0.5]
        )
    )
```

### nD space with weight function exp(-r)

 * via [Stroud](https://books.google.de/books/about/Approximate_calculation_of_multiple_inte.html?id=L_tQAAAAMAAJ&redir_esc=y) (1971):
   - [Stroud-Secrest](https://doi.org/10.1090/S0025-5718-1963-0161473-0) (1963, 4 schemes up to degree 5)
   - 2 schemes up to degree 5

Example:
```python
import quadpy

dim = 4
scheme = quadpy.enr.stroud_5_4(dim)
val = scheme.integrate(lambda x: x[0]**2)
```

### nD space with weight function exp(-r<sup>2</sup>)

 * via [Stroud](https://books.google.de/books/about/Approximate_calculation_of_multiple_inte.html?id=L_tQAAAAMAAJ&redir_esc=y) (1971):
   - [Stroud-Secrest](https://doi.org/10.1090/S0025-5718-1963-0161473-0) (1963, 4 schemes up to degree 5)
   - [Stroud](https://doi.org/10.1007/BF02162160) (1967, 2 schemes of degree 5)
   - [Stroud](https://doi.org/10.1137/0704004) (1967, 3 schemes of degree 7)
   - [Stenger](https://www.jstor.org/stable/2004361) (1971, 6 schemes up to degree 11, varying dimensionality restrictions)
   - 5 schemes up to degree 5

Example:
```python
import quadpy

dim = 4
scheme = quadpy.enr2.stroud_5_2(dim)
val = scheme.integrate(lambda x: x[0]**2)
```


### Installation

quadpy is [available from the Python Package Index](https://pypi.org/project/quadpy/), so with
```
pip3 install quadpy --user
```
you can install.

### Testing

To run the tests, just check out this repository and type
```
MPLBACKEND=Agg pytest
```

### License
quadpy is published under the [MIT license](https://en.wikipedia.org/wiki/MIT_License).

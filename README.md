# quadpy

Your one-stop shop for numerical integration in Python.

[![CircleCI](https://img.shields.io/circleci/project/github/nschloe/quadpy/master.svg)](https://circleci.com/gh/nschloe/quadpy/tree/master)
[![codecov](https://codecov.io/gh/nschloe/quadpy/branch/master/graph/badge.svg)](https://codecov.io/gh/nschloe/quadpy)
[![PyPi Version](https://img.shields.io/pypi/v/quadpy.svg)](https://pypi.python.org/pypi/quadpy)
[![awesome](https://img.shields.io/badge/awesome-yes-brightgreen.svg)](https://img.shields.io/badge/awesome-yes-brightgreen.svg)
[![GitHub stars](https://img.shields.io/github/stars/nschloe/quadpy.svg?style=social&label=Stars&maxAge=2592000)](https://github.com/nschloe/quadpy)

Hundreds of numerical integration schemes for
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
1D/2D/3D/nD spaces with weight functions exp(-r) and exp(-r<sup>2</sup>).

To numerically integrate any function over any given triangle, do
```python
import numpy
import quadpy

def f(x):
    return numpy.sin(x[0]) * numpy.sin(x[1])

triangle = numpy.array([[0.0, 0.0], [1.0, 0.0], [0.7, 0.5]])

val = quadpy.triangle.integrate(f, triangle, quadpy.triangle.Strang(9))
```
This uses Strang's rule of degree 6.

quadpy is fully vectorized, so if you like to compute the integral of a
function on many domains at once, you can provide them all in one `integrate()`
call, e.g.,
```python
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

### Adaptive quadrature

quadpy can do adaptive quadrature for certain domains.
Again, everything is fully vectorized, so you can provide multiple intervals
and vector-valued functions.

#### Line segments
```python
val, error_estimate = quadpy.line_segment.integrate_adaptive(
        lambda x: x * sin(5 * x),
        [0.0, pi],
        1.0e-10
        )
```

#### tanh-sinh quadrature

The more modern tanh-sinh quadrature is different from all other methods in
quadpy in that it doesn't exactly integrate any function exactly, not even
polynomials of low degree. Its tremendous usefulness rather comes from the fact
that a wide variety of function, even seemingly difficult ones with
(integrable) singularities at the end points, can be integrated with
_arbitrary_ precision.
```python
from mpmath import mp

mp.dps = 50

val, error_estimate = quadpy.line_segment.tanh_sinh(
        lambda x: mp.exp(t) * sympy.cos(t),
        [0, mp.pi/2],
        1.0e-50  # !
        )
```
Note the usage of `mpmath` here for arbirtrary precision arithmetics.

If the function has a singularity at a boundary, it needs to be shifted such
that the singularity is at 0. If there are singularities at both ends, the
function can be shifted both ways and be handed off to `tanh_sinh_lr`:
```
tanh_sinh_lr(f_left, f_right, interval_length, tol)
```


#### Triangles
```python
val, error_estimate = quadpy.triangle.integrate_adaptive(
        lambda x: x[0] * sin(5 * x[1]),
        [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]],
        1.0e-10
        )
```
_ProTip:_ You can provide many triangles that together form a domain to get an
approximation of the integral over the domain.

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
 * Gauss-Legendre (via
   [NumPy](https://docs.scipy.org/doc/numpy/reference/generated/numpy.polynomial.legendre.leggauss.html), arbitrary degree)
 * Gauss-Lobatto (arbitrary degree)
 * Gauss-Kronrod (after [Laurie](https://doi.org/10.1090/S0025-5718-97-00861-2), arbitrary degree)
 * [Gauss-Patterson](https://doi.org/10.1090/S0025-5718-68-99866-9) (7 schemes up to degree 191)
 * Gauss-Radau (arbitrary degree)
 * closed Newton-Cotes (arbitrary degree)
 * open Newton-Cotes (arbitrary degree)
 * [tanh-sinh quadrature](https://en.wikipedia.org/wiki/Tanh-sinh_quadrature)
   (see above)

You can use [orthopy](https://github.com/nschloe/orthopy) to generate Gauss
formulas for your own weight functions.

Example:
```python
val = quadpy.line_segment.integrate(
    lambda x: numpy.exp(x),
    [0.0, 1.0],
    quadpy.line_segment.GaussPatterson(5)
    )
```

### 1D half-space with weight function exp(-r)
<img src="https://nschloe.github.io/quadpy/e1r.png" width="50%">

 * Gauss-Laguerre (via
   [NumPy](https://docs.scipy.org/doc/numpy/reference/generated/numpy.polynomial.laguerre.laggauss.html), arbitrary degree)

Example:
```python
val = quadpy.e1r.integrate(
    lambda x: x**2,
    quadpy.e1r.GaussLaguerre(5)
    )
```


### 1D space with weight function exp(-r<sup>2</sup>)
<img src="https://nschloe.github.io/quadpy/e1r2.png" width="50%">

 * Gauss-Hermite (via
   [NumPy](https://docs.scipy.org/doc/numpy/reference/generated/numpy.polynomial.hermite.hermgauss.html), arbitrary degree)

Example:
```python
val = quadpy.e1r2.integrate(
    lambda x: x**2,
    quadpy.e1r2.GaussHermite(5)
    )
```

### Circle
<img src="https://nschloe.github.io/quadpy/circle.png" width="25%">

 * [Krylov](https://books.google.de/books/about/Approximate_Calculation_of_Integrals.html?id=ELeRwR27IRIC&redir_esc=y) (1959, arbitrary degree)

Example:
```python
val = quadpy.circle.integrate(
    lambda x: numpy.exp(x[0]),
    [0.0, 0.0], 1.0,
    quadpy.circle.Krylov(7)
    )
```

### Triangle
<img src="https://nschloe.github.io/quadpy/triangle.png" width="25%">

Apart from the classical centroid, vertex, and seven-point schemes we have

 * [Hammer-Marlowe-Stroud](https://doi.org/10.1090/S0025-5718-1956-0086389-6)
   (1956, 5 schemes up to degree 5),
 * [Hammer-Stroud](https://doi.org/10.1090/S0025-5718-1958-0102176-6) (1958, 2 schemes up to degree 3)
 * open and closed Newton-Cotes schemes (1970, after [Silvester](https://doi.org/10.1090/S0025-5718-1970-0258283-6), arbitrary degree),
 * via [Stroud](https://books.google.de/books/about/Approximate_calculation_of_multiple_inte.html?id=L_tQAAAAMAAJ&redir_esc=y) (1971):
   - [Albrecht-Collatz](https://doi.org/10.1002/zamm.19580380102) (1958, degree 3)
   - conical product scheme (degree 7)
 * [Strang](http://bookstore.siam.org/wc08/)/[Cowper](https://doi.org/10.1002/nme.1620070316) (1973, 10 schemes up to
   degree 7),
 * [Lyness-Jespersen](https://doi.org/10.1093/imamat/15.1.19) (1975, 21
   schemes up to degree 11),
 * [Lether](https://doi.org/10.1016/0771-050X(76)90008-5) (1976, degree 2n-2, arbitrary n, not symmetric)
 * [Hillion](https://doi.org/10.1002/nme.1620110504) (1977, 10 schemes up to
   degree 3),
 * [Grundmann-Möller](https://doi.org/10.1137/0715019) (1978, arbitrary degree),
 * [Laursen-Gellert](https://doi.org/10.1002/nme.1620120107) (1978, 17
   schemes up to degree 10),
 * [CUBTRI](http://dl.acm.org/citation.cfm?id=356001) (Laurie, 1982, degree 8),
 * [TRIEX](http://dl.acm.org/citation.cfm?id=356070) (de Doncker-Robinson, 1984, degrees 9 and 11),
 * [Dunavant](https://doi.org/10.1002/nme.1620210612) (1985, 20 schemes up
   to degree 20),
 * [Cools-Haegemans](https://lirias.kuleuven.be/handle/123456789/131869) (1987,
   degrees 8 and 11),
 * [Gatermann](https://doi.org/10.1007/BF02251251) (1988, degree 7)
 * Berntsen-Espelid (1990, 4 schemes of degree 13, the first one being
   [DCUTRI](http://dl.acm.org/citation.cfm?id=131772)),
 * [Liu-Vinokur](https://doi.org/10.1006/jcph.1998.5884) (1998, 13 schemes
   up to degree 5),
 * [Walkington](http://www.math.cmu.edu/~nw0z/publications/00-CNA-023/023abs/)
   (2000, 5 schemes up to degree 5),
 * [Wandzura-Xiao](https://doi.org/10.1016/S0898-1221(03)90004-6) (2003, 6
   schemes up to degree 30),
 * [Taylor-Wingate-Bos](https://arxiv.org/abs/math/0501496) (2005, 5 schemes up
   to degree 14),
 * [Zhang-Cui-Liu](http://www.jstor.org/stable/43693493) (2009, 3 schemes up to
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
val = quadpy.triangle.integrate(
    lambda x: numpy.exp(x[0]),
    [[0.0, 0.0], [1.0, 0.0], [0.5, 0.7]],
    quadpy.triangle.XiaoGimbutas(5)
    )
```

### Disk
<img src="https://nschloe.github.io/quadpy/disk.png" width="25%">

 * [Peirce](http://www.jstor.org/stable/2098722) (1957, arbitrary degree)
 * via [Stroud](https://books.google.de/books/about/Approximate_calculation_of_multiple_inte.html?id=L_tQAAAAMAAJ&redir_esc=y):
   - [Radon](https://eudml.org/doc/176796) (1948, degree 5)
   - [Peirce](https://books.google.de/books/about/Numerical_integration_over_planar_region.html?id=WR9SAAAAMAAJ&redir_esc=y)
     (1956, 3 schemes up to degree 11)
   - [Albrecht-Collatz](https://doi.org/10.1002/zamm.19580380102) (1958, degree 3)
   - [Hammer-Stroud](https://doi.org/10.1090/S0025-5718-1958-0102176-6) (1958, 8 schemes up to degree 15)
   - [Albrecht](https://doi.org/10.1002/zamm.19600401014) (1960, 8 schemes up to degree 17)
   - Mysovskih (1964, 3 schemes up to degree 15)
   - [Rabinowitz-Richter](https://doi.org/10.2307/2004962) (1969, 6 schemes up to degree 15)
 * [Lether](http://www.jstor.org/stable/2949473) (1971, arbitrary degree)
 * [Cools-Haegemans](https://lirias.kuleuven.be/handle/123456789/131870) (1985, 3 schemes up to degree 9)
 * [Wissmann-Becker](https://doi.org/10.1137/0723043) (1986, 3 schemes up to degree 8)
 * [Cools-Kim](https://link.springer.com/article/10.1007/BF03012263) (2000, 3 schemes up to degree 21)

Example:
```python
val = quadpy.disk.integrate(
    lambda x: numpy.exp(x[0]),
    [0.0, 0.0], 1.0,
    quadpy.disk.Lether(6)
    )
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
 * [Cools-Haegemans](https://lirias.kuleuven.be/handle/123456789/131870) (1985, 3 schemes up to degree 13)
 * [Dunavant](https://doi.org/10.1002/nme.1620211004) (1985, 11 schemes up to degree 19)
 * [Morrow-Patterson](https://doi.org/10.1137/0722071) (1985, 2 schemes up to degree 20, single precision)
 * [Wissmann-Becker](https://doi.org/10.1137/0723043) (1986, 6 schemes up to degree 8)
 * [Cools-Haegemans](https://doi.org/10.1007/BF02247942) (1988, 2 schemes up to degree 13)
 * products of line segment schemes
 * all formulas from the n-cube

Example:
```python
val = quadpy.quadrilateral.integrate(
    lambda x: numpy.exp(x[0]),
    [[[0.0, 0.0], [1.0, 0.0]], [[0.0, 1.0], [1.0, 1.0]]],
    quadpy.quadrilateral.Stroud('C2 7-2')
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

Example:
```python
val = quadpy.e2r.integrate(
    lambda x: x[0]**2,
    quadpy.e2r.RabinowitzRichter(5)
    )
```


### 2D space with weight function exp(-r<sup>2</sup>)
<img src="https://nschloe.github.io/quadpy/e2r2.png" width="25%">

 * via [Stroud](https://books.google.de/books/about/Approximate_calculation_of_multiple_inte.html?id=L_tQAAAAMAAJ&redir_esc=y) (1971):
   - [Stroud-Secrest](https://doi.org/10.1090/S0025-5718-1963-0161473-0) (1963, 2 schemes up to degree 7)
   - [Rabinowitz-Richter](https://doi.org/10.2307/2004962) (1969, 5 schemes up to degree 15)
   - 3 schemes up to degree 7

Example:
```python
val = quadpy.e2r2.integrate(
    lambda x: x[0]**2,
    quadpy.e2r2.RabinowitzRichter(3)
    )
```


### Sphere
<img src="https://nschloe.github.io/quadpy/sphere.png" width="25%">

 * via [Stroud](https://books.google.de/books/about/Approximate_calculation_of_multiple_inte.html?id=L_tQAAAAMAAJ&redir_esc=y) (1971):
   - [Albrecht-Collatz](https://doi.org/10.1002/zamm.19580380102) (1958, 5
     schemes up to degree 7)
   - [McLaren](https://doi.org/10.1090/S0025-5718-1963-0159418-2) (1963, 10 schemes up to degree 14)
 * [Lebedev](https://en.wikipedia.org/wiki/Lebedev_quadrature) (1976, 32
   schemes up to degree 131)
 * [Heo-Xu](https://doi.org/10.1090/S0025-5718-00-01198-4) (2001, 27 schemes up
   to degree 39, single-precision)

Example:
```python
val = quadpy.sphere.integrate(
    lambda x: numpy.exp(x[0]),
    [0.0, 0.0, 0.0], 1.0,
    quadpy.sphere.Lebedev(19)
    )
```
Integration on the sphere can also be done for function defined in spherical
coordinates:
```python
val = quadpy.sphere.integrate_spherical(
    lambda azimuthal, polar: numpy.sin(azimuthal)**2 * numpy.sin(polar),
    rule=quadpy.sphere.Lebedev(19)
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
val = quadpy.ball.integrate(
    lambda x: numpy.exp(x[0]),
    [0.0, 0.0, 0.0], 1.0,
    quadpy.ball.HammerStroud('14-3a')
    )
```


### Tetrahedron
<img src="https://nschloe.github.io/quadpy/tet.png" width="25%">

 * [Hammer-Marlowe-Stroud](https://doi.org/10.1090/S0025-5718-1956-0086389-6)
   (1956, 3 schemes up to degree 3)
 * [Hammer-Stroud](https://doi.org/10.1090/S0025-5718-1958-0102176-6) (1958, 2 schemes up to degree 3)
 * open and closed Newton-Cotes (1970, after [Silvester](https://doi.org/10.1090/S0025-5718-1970-0258283-6)) (arbitrary degree)
 * [Stroud](https://cds.cern.ch/record/104291?ln=en) (1971, degree 7)
 * [Grundmann-Möller](https://doi.org/10.1137/0715019) (1978, arbitrary degree),
 * [Yu](https://doi.org/10.1016/0045-7825(84)90072-0) (1984, 5 schemes up to degree 6)
 * [Keast](https://doi.org/10.1016/0045-7825(86)90059-9) (1986, 11 schemes up to
   degree 8)
 * [Beckers-Haegemans](https://lirias.kuleuven.be/handle/123456789/132648) (1990, degrees 8 and 9)
 * [Gatermann](https://doi.org/10.1007/978-94-011-2646-5_2) (1992, degree 5)
 * [Liu-Vinokur](https://doi.org/10.1006/jcph.1998.5884) (1998, 14 schemes up to
   degree 5)
 * [Walkington](http://www.math.cmu.edu/~nw0z/publications/00-CNA-023/023abs/)
   (2000, 6 schemes up to degree 7)
 * [Zienkiewicz](http://www.sciencedirect.com/science/book/9780750664318)
   (2005, 2 schemes up to degree 3)
 * [Zhang-Cui-Liu](http://www.jstor.org/stable/43693493) (2009, 2 schemes up to
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
val = quadpy.tetrahedron.integrate(
    lambda x: numpy.exp(x[0]),
    [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 0.7, 0.0], [0.3, 0.9, 1.0]],
    quadpy.tetrahedron.Keast(10)
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
val = quadpy.hexahedron.integrate(
    lambda x: numpy.exp(x[0]),
    quadpy.hexahedron.cube_points([0.0, 1.0], [-0.3, 0.4], [1.0, 2.1]),
    quadpy.hexahedron.Product(quadpy.line_segment.NewtonCotesClosed(3))
    )
```

### Pyramid
<img src="https://nschloe.github.io/quadpy/pyra.png" width="25%">

 * [Felippa](https://doi.org/10.1108/02644400410554362) (9 schemes
   up to degree 5)

Example:
```python
val = quadpy.pyramid.integrate(
    lambda x: numpy.exp(x[0]),
    [
      [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 0.7, 0.0], [0.3, 0.9, 0.0],
      [0.0, 0.1, 1.0],
    ],
    quadpy.pyramid.Felippa(5)
    )
```

### Wedge
<img src="https://nschloe.github.io/quadpy/wedge.png" width="15%">

 * [Felippa](https://doi.org/10.1108/02644400410554362) (6 schemes
   up to degree 6)

Example:
```python
val = quadpy.wedge.integrate(
    lambda x: numpy.exp(x[0]),
    [
      [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 0.7, 0.0]],
      [[0.0, 0.0, 1.0], [1.0, 0.0, 1.0], [0.5, 0.7, 1.0]],
    ],
    quadpy.wedge.Felippa(3)
    )
```


### 3D space with weight function exp(-r)
<img src="https://nschloe.github.io/quadpy/e3r.png" width="25%">

 * via [Stroud](https://books.google.de/books/about/Approximate_calculation_of_multiple_inte.html?id=L_tQAAAAMAAJ&redir_esc=y) (1971):
   - [Stroud-Secrest](https://doi.org/10.1090/S0025-5718-1963-0161473-0) (1963, 5 schemes up to degree 7)

Example:
```python
val = quadpy.e2r.integrate(
    lambda x: x[0]**2,
    quadpy.e2r.StroudSecrest('IX')
    )
```


### 3D space with weight function exp(-r<sup>2</sup>)
<img src="https://nschloe.github.io/quadpy/e3r2.png" width="25%">

 * via [Stroud](https://books.google.de/books/about/Approximate_calculation_of_multiple_inte.html?id=L_tQAAAAMAAJ&redir_esc=y) (1971):
   - [Stroud-Secrest](https://doi.org/10.1090/S0025-5718-1963-0161473-0) (1963, 7 schemes up to degree 7)
   - scheme of degree 14

Example:
```python
val = quadpy.e2r2.integrate(
    lambda x: x[0]**2,
    quadpy.e2r2.RabinowitzRichter(3)
    )
```

### n-Simplex
 * via [Stroud](https://books.google.de/books/about/Approximate_calculation_of_multiple_inte.html?id=L_tQAAAAMAAJ&redir_esc=y):
   - [Lauffer](https://doi.org/10.1007/BF01900222) (1955, 5 schemes up to degree 5)
   - [Hammer-Stroud](https://doi.org/10.2307/2002484) (1956, 3 schemes up to degree 3)
   - [Stroud](https://doi.org/10.2307/2002945) (1964, degree 3)
   - [Stroud](https://doi.org/10.1007/BF02165227) (1966, 7 schemes of degree 3)
   - [Stroud](https://doi.org/10.1137/0706009) (1969, degree 5)
 * [Grundmann-Möller](https://doi.org/10.1137/0715019) (1978, arbitrary degree)
 * [Walkington](http://www.math.cmu.edu/~nw0z/publications/00-CNA-023/023abs/) (2000, 5 schemes up to degree 7)

Example:
```python
dim = 4
val = quadpy.simplex.integrate(
    lambda x: numpy.exp(x[0]),
    numpy.array([
        [0.0, 0.0, 0.0, 0.0],
        [1.0, 2.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0],
        [0.0, 3.0, 1.0, 0.0],
        [0.0, 0.0, 4.0, 1.0],
        ]),
    quadpy.simplex.GrundmannMoeller(dim, 3)
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
dim = 4
quadpy.nsphere.integrate(
    lambda x: numpy.exp(x[0]),
    numpy.zeros(dim), 1.0,
    quadpy.nsphere.Dobrodeev1978(dim)
    )
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
dim = 4
quadpy.nball.integrate(
    lambda x: numpy.exp(x[0]),
    numpy.zeros(dim), 1.0,
    quadpy.nball.Dobrodeev1970(dim)
    )
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
dim = 4
quadpy.ncube.integrate(
    lambda x: numpy.exp(x[0]),
    quadpy.ncube.ncube_points(
        [0.0, 1.0], [0.1, 0.9], [-1.0, 1.0], [-1.0, -0.5]
        ),
    quadpy.ncube.Stroud(dim, 'Cn 3-3')
    )
```

### nD space with weight function exp(-r)

 * via [Stroud](https://books.google.de/books/about/Approximate_calculation_of_multiple_inte.html?id=L_tQAAAAMAAJ&redir_esc=y) (1971):
   - [Stroud-Secrest](https://doi.org/10.1090/S0025-5718-1963-0161473-0) (1963, 4 schemes up to degree 5)
   - 2 schemes up to degree 5

Example:
```python
dim = 4
val = quadpy.enr.integrate(
    lambda x: x[0]**2,
    quadpy.enr.Stroud(dim, '5-4')
    )
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
dim = 4
val = quadpy.enr2.integrate(
    lambda x: x[0]**2,
    quadpy.enr2.Stroud(dim, '5-2')
    )
```

### Installation

quadpy is [available from the Python Package Index](https://pypi.python.org/pypi/quadpy/), so with
```
pip install -U quadpy
```
you can install/upgrade.

### Testing

To run the tests, just check out this repository and type
```
MPLBACKEND=Agg pytest
```

### Distribution

To create a new release

1. bump the `__version__` number,

2. publish to PyPi and GitHub:
    ```
    $ make publish
    ```

### License
quadpy is published under the [MIT license](https://en.wikipedia.org/wiki/MIT_License).

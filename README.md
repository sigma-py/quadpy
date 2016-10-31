# quadrature

Your one-stop shop for numerical integration in Python.

[![Build Status](https://travis-ci.org/nschloe/quadrature.svg?branch=master)](https://travis-ci.org/nschloe/quadrature)
[![Code Health](https://landscape.io/github/nschloe/quadrature/master/landscape.png)](https://landscape.io/github/nschloe/quadrature/master)
[![codecov](https://codecov.io/gh/nschloe/quadrature/branch/master/graph/badge.svg)](https://codecov.io/gh/nschloe/quadrature)
[![PyPi Version](https://img.shields.io/pypi/v/quadrature.svg)](https://pypi.python.org/pypi/quadrature)
[![GitHub stars](https://img.shields.io/github/stars/nschloe/quadrature.svg?style=social&label=Star&maxAge=2592000)](https://github.com/nschloe/quadrature)

Hundreds of numerical integration schemes for line segments, circles, disks,
triangles, quadrilaterals, spheres, tetrahedra, hexahedra, wedges, pyramids.

To numerically integrate any function over any given triangle, do
```python
import numpy
import quadrature

def f(x):
    return numpy.sin(x[0]) * numpy.sin(x[1])

triangle = numpy.array([[0.0, 0.0], [1.0, 0.0], [0.7, 0.5]])

val = quadrature.triangle.integrate(f, triangle, quadrature.triangle.Strang(9))
```
This uses Strang's rule of degree 6.

## Schemes

### Line segment
![](https://nschloe.github.io/quadrature/line.svg)

 * Chebyshev-Gauß (both variants, arbitrary order)
 * Clenshaw-Curtis (after
   [Waldvogel](https://dx.doi.org/10.1007/s10543-006-0045-4), arbitrary order)
 * Fejér-type-1 (after
   [Waldvogel](https://dx.doi.org/10.1007/s10543-006-0045-4), arbitrary order)
 * Fejér-type-2 (after
   [Waldvogel](https://dx.doi.org/10.1007/s10543-006-0045-4), arbitrary order)
 * Gauß-Hermite (via
   [NumPy](https://docs.scipy.org/doc/numpy/reference/generated/numpy.polynomial.hermite.hermgauss.html), arbitrary order)
 * Gauß-Laguerre (via
   [NumPy](https://docs.scipy.org/doc/numpy/reference/generated/numpy.polynomial.laguerre.laggauss.html), arbitrary order)
 * Gauß-Legendre (via
   [NumPy](https://docs.scipy.org/doc/numpy/reference/generated/numpy.polynomial.legendre.leggauss.html), arbitrary order)
 * Gauß-Lobatto (arbitrary order)
 * Gauß-Patterson (7 schemes up to degree 191)
 * Gauß-Radau (arbitrary order)
 * closed Newton-Cotes (arbitray order)
 * open Newton-Cotes (arbitray order)

### Circle
![](https://nschloe.github.io/quadrature/circle.png)

 * equidistant points

### Triangle
![](https://nschloe.github.io/quadrature/tri.png)

Apart from the classical centroid, vertex, and seven-point schemes we have

 * [Hammer-Marlowe-Stroud](https://doi.org/10.1090/S0025-5718-1956-0086389-6)
   (1956, 5 schemes up to degree 5),
 * open and closed Newton-Cotes schemes (1970, after [Silvester](https://doi.org/10.1090/S0025-5718-1970-0258283-6), arbitrary degree),
 * [Strang](http://bookstore.siam.org/wc08/)/[Cowper](https://dx.doi.org/10.1002/nme.1620070316) (1973, 10 schemes up to
   degree 7),
 * [Lyness-Jespersen](https://dx.doi.org/10.1093/imamat/15.1.19) (1975, 21
   schemes up to degree 11),
 * [Hillion](https://dx.doi.org/10.1002/nme.1620110504) (1977),
 * [Laursen-Gellert](https://dx.doi.org/10.1002/nme.1620120107) (1978, 17
   schemes up to degree 10),
 * [CUBTRI](http://dl.acm.org/citation.cfm?id=356001) (1982, degree 8),
 * [TRIEX](http://dl.acm.org/citation.cfm?id=356070) (1984, degrees 9 and 11),
 * [Dunavant](https://dx.doi.org/10.1002/nme.1620210612) (1985, 20 schemes up
   to degree 20),
 * [Cools-Haegemans](https://lirias.kuleuven.be/handle/123456789/131869) (1987,
   degrees 8 and 11),
 * Berntsen-Espelid (1990, 4 schemes of degree 13, the first one being
   [DCUTRI](http://dl.acm.org/citation.cfm?id=131772)),
 * [Liu-Vinokur](https://dx.doi.org/10.1006/jcph.1998.5884) (1998, 13 schemes
   up to degree 5),
 * [Wandzura-Xiao](https://dx.doi.org/10.1016/S0898-1221(03)90004-6) (2003, 6
   schemes up to degree 30),
 * [Taylor-Wingate-Bos](https://arxiv.org/abs/math/0501496) (2005, 5 schemes up
   to degree 14),
 * [Zhang-Cui-Liu](http://www.jstor.org/stable/43693493) (2009, 3 schemes up to
   degree 20),
 * [Xiao-Gimbutas](http://dx.doi.org/10.1016/j.camwa.2009.10.027) (2010, 50
   schemes up to degree 50).

### Disk
![](https://nschloe.github.io/quadrature/disk.png)

 * [Peirce](http://www.jstor.org/stable/2098722) (1957, arbitrary degree)
 * [Lether](http://www.jstor.org/stable/2949473) (1971, arbitrary degree)

### Quadrilateral
![](https://nschloe.github.io/quadrature/quad.png)

 * Product schemes derived from line segment schemes
 * [Stroud's schemes](https://books.google.de/books/about/Approximate_calculation_of_multiple_inte.html?id=L_tQAAAAMAAJ&redir_esc=y) (6 schemes up to degree 15)

### Tetrahedron
![](https://nschloe.github.io/quadrature/tet.png)

 * [Hammer-Marlowe-Stroud](https://doi.org/10.1090/S0025-5718-1956-0086389-6)
   (1956, 3 schemes up to degree 3)
 * open and closed Newton-Cotes (1970, after [Silvester](https://doi.org/10.1090/S0025-5718-1970-0258283-6)) (arbitrary degree)
 * [Yu](http://dx.doi.org/10.1016/0045-7825(84)90072-0) (1984, 5 schemes up to degree 6)
 * [Keast](http://dx.doi.org/10.1016/0045-7825(86)90059-9) (1986, 11 schemes up to
   degree 8)
 * [Liu-Vinokur](http://dx.doi.org/10.1006/jcph.1998.5884) (1998, 14 schemes up to
   degree 5)
 * [Zienkiewicz](http://www.sciencedirect.com/science/book/9780750664318)
   (2005, 2 schemes up to degree 3)
 * [Zhang-Cui-Liu](http://www.jstor.org/stable/43693493) (2009, 2 schemes up to
   degree 14)
 * [Xiao-Gimbutas](http://dx.doi.org/10.1016/j.camwa.2009.10.027) (2010, 15
   schemes up to degree 15)
 * [Shunn-Ham](http://dx.doi.org/10.1016/j.cam.2012.03.032) (2012, 6 schemes up to
   degree 7)

### Hexahedron
![](https://nschloe.github.io/quadrature/hexa.png)

 * Product schemes derived from line segment schemes

### Pyramid
![](https://nschloe.github.io/quadrature/pyra.png)

 * [Felippa's schemes](http://dx.doi.org/10.1108/02644400410554362) (9 schemes
   up to degree 5)

### Wedge
![](https://nschloe.github.io/quadrature/wedge.png)

 * [Felippa's schemes](http://dx.doi.org/10.1108/02644400410554362) (6 schemes
   up to degree 6)

### Sphere
![](https://nschloe.github.io/quadrature/sphere.png)

 * [Lebedev's schemes](https://en.wikipedia.org/wiki/Lebedev_quadrature) (32
   schemes up to degree 131)

### Installation

#### Python Package Index

quadrature is [available from the Python Package Index](https://pypi.python.org/pypi/quadrature/), so with
```
pip install -U quadrature
```
you can install/upgrade.

#### Manual installation

Download quadrature from
[the Python Package Index](https://pypi.python.org/pypi/quadrature/).
Place the quadrature script in a directory where Python can find it (e.g.,
`$PYTHONPATH`). You can install it system-wide with
```
python setup.py install
```

### Testing

To run the tests, just check out this repository and type
```
pytest
```

### Distribution

To create a new release

1. bump the `__version__` number,

2. publish to PyPi and GitHub:
    ```
    $ make publish
    ```

### License
quadrature is published under the [MIT license](https://en.wikipedia.org/wiki/MIT_License).

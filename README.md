# quadrature

Your one-stop shop for numerical integration in Python.

[![Build Status](https://travis-ci.org/nschloe/quadrature.svg?branch=master)](https://travis-ci.org/nschloe/quadrature)
[![Code Health](https://landscape.io/github/nschloe/quadrature/master/landscape.png)](https://landscape.io/github/nschloe/quadrature/master)
[![codecov](https://codecov.io/gh/nschloe/quadrature/branch/master/graph/badge.svg)](https://codecov.io/gh/nschloe/quadrature)
[![PyPi Version](https://img.shields.io/pypi/v/quadrature.svg)](https://pypi.python.org/pypi/quadrature)
[![GitHub stars](https://img.shields.io/github/stars/nschloe/quadrature.svg?style=social&label=Star&maxAge=2592000)](https://github.com/nschloe/quadrature)

![](https://nschloe.github.io/quadrature/s9.png)

Numerical integration schemes for lines, triangles, tetrahedra,
quadrilaterals, hexahedra, wedges, pyramids, circles, and spheres.

To numerically integrate any function over any given triangle, do
```python
import numpy
import quadrature

def f(x): return numpy.sin(x[0]) * numpy.sin(x[1])

triangle = numpy.array([[0.0, 0.0], [1.0, 0.0], [0.7, 0.5]])

val = quadrature.triangle.integrate(f, triangle, quadrature.triangle.Strang9())
```
This uses Strang's rule of degree 6 (see picture above); [many more are
implemented](https://people.sc.fsu.edu/~jburkardt/datasets/quadrature_rules_tri/quadrature_rules_tri.html).

### Rules

##### Line segment
![](https://nschloe.github.io/quadrature/line.png)

##### Triangle
![](https://nschloe.github.io/quadrature/tri.png)

 * [Strang's schemes](http://bookstore.siam.org/wc08/),
 * [CUBTRI](http://dl.acm.org/citation.cfm?id=356001),
 * [TRIEX](http://dl.acm.org/citation.cfm?id=356070&CFID=836775288&CFTOKEN=89206835),
 * [DCUTRI](http://dl.acm.org/citation.cfm?id=131772),
 * [Dunavant's schemes](https://dx.doi.org/10.1002/nme.1620210612),
 * [Zhang-Cui-Liu](http://www.jstor.org/stable/43693493),
 * [Wandzura-Xiao](https://dx.doi.org/10.1016/S0898-1221(03)90004-6),
 * [Lyness-Jespersen](https://dx.doi.org/10.1093/imamat/15.1.19),
 * open and closed Newton-Cotes schemes (arbitrary degree),
 * [Taylor-Wingate-Bos](https://arxiv.org/abs/math/0501496),
 * Berntsen-Espelid (three degree-13 schemes),
 * [Hammer-Marlowe-Stroud](https://doi.org/10.1090/S0025-5718-1956-0086389-6),
 * [Cowper](https://dx.doi.org/10.1002/nme.1620070316),
 * [Liu-Vinokur](https://dx.doi.org/10.1006/jcph.1998.5884),
 * [Hillion](https://dx.doi.org/10.1002/nme.1620110504),
 * [Cools-Haegemans](https://lirias.kuleuven.be/handle/123456789/131869),
 * [Laursen-Gellert](https://dx.doi.org/10.1002/nme.1620120107).

##### Circle
![](https://nschloe.github.io/quadrature/circle.png)

##### Quadrilateral
![](https://nschloe.github.io/quadrature/quad.png)

##### Tetrahedron
![](https://nschloe.github.io/quadrature/tet.png)

##### Hexahedron
![](https://nschloe.github.io/quadrature/hexa.png)

##### Pyramid
![](https://nschloe.github.io/quadrature/pyra.png)

##### Sphere
![](https://nschloe.github.io/quadrature/sphere.png)

##### Wedge
![](https://nschloe.github.io/quadrature/wedge.png)

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

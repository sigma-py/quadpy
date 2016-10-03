# quadrature

[![Build Status](https://travis-ci.org/nschloe/quadrature.svg?branch=master)](https://travis-ci.org/nschloe/quadrature)
[![Code Health](https://landscape.io/github/nschloe/quadrature/master/landscape.png)](https://landscape.io/github/nschloe/quadrature/master)
[![codecov](https://codecov.io/gh/nschloe/quadrature/branch/master/graph/badge.svg)](https://codecov.io/gh/nschloe/quadrature)
[![PyPi Version](https://img.shields.io/pypi/v/quadrature.svg)](https://pypi.python.org/pypi/quadrature)
[![GitHub stars](https://img.shields.io/github/stars/nschloe/quadrature.svg?style=social&label=Star&maxAge=2592000)](https://github.com/nschloe/quadrature)

![](https://nschloe.github.io/quadrature/s9.png)

Numerical integration schemes for triangles.

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
nosetests
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

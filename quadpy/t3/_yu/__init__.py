import pathlib
from sympy import Rational as frac
from sympy import sqrt

from ...helpers import article
from .._helpers import _read, T3Scheme, concat, s4, s31

source = article(
    authors=["Yu Jinyun"],
    title="Symmetyric Gaussian quadrature formulae for tetrahedronal regions",
    journal="Computer Methods in Applied Mechanics and Engineering",
    volume="43",
    year="1984",
    pages="349-353",
    url="https://doi.org/10.1016/0045-7825(84)90072-0",
)

this_dir = pathlib.Path(__file__).resolve().parent


def yu_2():
    degree = 2
    weights, points = s31([frac(1, 4), frac(1, 4) - sqrt(5) * frac(1, 20)])
    return T3Scheme("Yu 1", weights, points, degree, source)


def yu_3():
    degree = 3
    weights, points = concat(s4(-frac(4, 5)), s31([frac(9, 20), frac(1, 6)]))
    return T3Scheme("Yu 2", weights, points, degree, source)


def yu_4():
    return _read(this_dir / "yu_4.json", source)


def yu_5():
    return _read(this_dir / "yu_5.json", source)


def yu_6():
    return _read(this_dir / "yu_6.json", source)

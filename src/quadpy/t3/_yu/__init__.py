import pathlib

from sympy import Rational as frac
from sympy import sqrt

from ...helpers import article
from .._helpers import T3Scheme, _read, register

source = article(
    authors=["Yu Jinyun"],
    title="Symmetyric Gaussian quadrature formulae for tetrahedronal regions",
    journal="Computer Methods in Applied Mechanics and Engineering",
    volume="43",
    year="1984",
    pages="349-353",
    url="https://doi.org/10.1016/0045-7825%2884%2990072-0",
)

this_dir = pathlib.Path(__file__).resolve().parent


def yu_2():
    degree = 2
    d = {"s31": [[frac(1, 4)], [frac(1, 4) - sqrt(5) * frac(1, 20)]]}
    return T3Scheme("Yu 1", d, degree, source)


def yu_3():
    degree = 3
    d = {"s4": [[-frac(4, 5)]], "s31": [[frac(9, 20)], [frac(1, 6)]]}
    return T3Scheme("Yu 2", d, degree, source)


def yu_4():
    return _read(this_dir / "yu_4.json", source)


def yu_5():
    return _read(this_dir / "yu_5.json", source)


def yu_6():
    return _read(this_dir / "yu_6.json", source)


register([yu_2, yu_3, yu_4, yu_5, yu_6])

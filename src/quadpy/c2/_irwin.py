from sympy import Rational as frac

from ..helpers import book
from ._helpers import C2Scheme, register

source = book(
    authors=["Joseph Oscar Irwin"],
    title="On quadrature and cubature",
    publisher="Cambridge University Press",
    year="1923",
    url="https://books.google.de/books/about/On_quadrature_and_cubature.html?id=SuruAAAAMAAJ",
)


def irwin_1():
    d = {"d4_aa": [[frac(14, 48)], [1]], "d4_ab": [[-frac(1, 48)], [3], [1]]}
    return C2Scheme("Irwin 1", d, 3, source)


def irwin_2():
    d = {
        "d4_aa": [[frac(889, 2880), frac(5, 2880)], [1, 3]],
        "d4_ab": [[-frac(98, 2880), frac(11, 2880)], [3, 5], [1, 1]],
    }
    return C2Scheme("Irwin 2", d, 5, source, 5.685e-14)


register([irwin_1, irwin_2])

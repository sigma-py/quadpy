from sympy import Rational as frac

from ..helpers import book
from ._helpers import C2Scheme, concat, symm_s, symm_s_t

source = book(
    authors=["Joseph Oscar Irwin"],
    title="On quadrature and cubature",
    publisher="Cambridge University Press",
    year="1923",
    url="https://books.google.de/books/about/On_quadrature_and_cubature.html?id=SuruAAAAMAAJ",
)


def irwin_1():
    weights, points = concat(symm_s([frac(14, 48), 1]), symm_s_t([-frac(1, 48), 3, 1]))
    return C2Scheme("Irwin 1", weights, points, 3, source)


def irwin_2():
    weights, points = concat(
        symm_s([frac(889, 2880), 1], [frac(5, 2880), 3]),
        symm_s_t([-frac(98, 2880), 3, 1], [frac(11, 2880), 5, 1]),
    )
    return C2Scheme("Irwin 2", weights, points, 5, source, 5.685e-14)

import numpy
import sympy

from ..helpers import article
from ._helpers import S2Scheme, expand_symmetries

_source = article(
    authors=["J. Radon"],
    title="Zur mechanischen Kubatur",
    journal="Monatshefte für Mathematik",
    year="1948",
    volume="52",
    pages="286-300",
    issn="0026-9255",
    issne="1436-5081/e",
    url="https://eudml.org/doc/176796",
)

frac = sympy.Rational
pi = sympy.pi
sqrt = numpy.vectorize(sympy.sqrt)


def radon(alpha):
    r = sqrt(frac(alpha + 4, alpha + 6))
    s = sqrt(frac(alpha + 4, 4 * (alpha + 6)))
    t = sqrt(frac(3 * (alpha + 4), 4 * (alpha + 6)))

    A = frac(4, (alpha + 4) ** 2)
    B = frac((alpha + 2) * (alpha + 6), 6 * (alpha + 4) ** 2)

    d = {
        "zero": [[A]],
        "c2_a0": [[B], [r]],
        # ERR Stroud is missing +- in front of t.
        "pm2": [[B], [s], [t]],
    }
    points, weights = expand_symmetries(d)
    return S2Scheme(f"Radon({alpha})", weights, points, 5, _source)

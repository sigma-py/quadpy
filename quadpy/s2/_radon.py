import numpy as np
import sympy

from ..helpers import article
from ._helpers import S2Scheme, register

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
sqrt = np.vectorize(sympy.sqrt)


def radon(alpha=0):
    r = sqrt(frac(alpha + 4, alpha + 6))
    s = sqrt(frac(alpha + 4, 4 * (alpha + 6)))
    t = sqrt(frac(3 * (alpha + 4), 4 * (alpha + 6)))

    A = frac(4, (alpha + 4) ** 2)
    B = frac((alpha + 2) * (alpha + 6), 6 * (alpha + 4) ** 2)

    d = {
        "zero2": [[A]],
        "c2_a0": [[B], [r]],
        # ERR Stroud is missing +- in front of t.
        "sxy": [[B], [s], [t]],
    }
    return S2Scheme(f"Radon({alpha})", d, 5, _source)


register([radon])

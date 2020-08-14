import pathlib

from sympy import Rational as frac
from sympy import sqrt

from ...helpers import article
from .._helpers import C2Scheme, _read, register

source = article(
    authors=["H.J. Schmid"],
    title="On cubature formulae with a minimal number of knots",
    journal="Numerische Mathematik",
    month="sep",
    year="1978",
    volume="31",
    number="3",
    pages="281–297",
    url="https://eudml.org/doc/132580",
)

this_dir = pathlib.Path(__file__).resolve().parent


def schmid_2():
    d = {
        "plain": [
            [frac(1, 4), frac(1, 4), frac(1, 2)],
            [-sqrt(frac(1, 3)), -sqrt(frac(1, 3)), +sqrt(frac(1, 3))],
            [+sqrt(frac(2, 3)), -sqrt(frac(2, 3)), 0],
        ]
    }
    return C2Scheme("Schmid 2", d, 2, source, 4.441e-16)


def schmid_4():
    d = {
        "plain": [
            [
                frac(2, 9) - 2 * sqrt(5) / 45,
                frac(2, 9) + 2 * sqrt(5) / 45,
                frac(5, 36) + 5 * sqrt(29) / 18 / 29,
                frac(5, 36) + 5 * sqrt(29) / 18 / 29,
                frac(5, 36) - 5 * sqrt(29) / 18 / 29,
                frac(5, 36) - 5 * sqrt(29) / 18 / 29,
            ],
            [0, 0, +sqrt(15) / 5, -sqrt(15) / 5, +sqrt(15) / 5, -sqrt(15) / 5],
            [
                (sqrt(3) + sqrt(15)) / 6,
                (sqrt(3) - sqrt(15)) / 6,
                (+sqrt(87) - 2 * sqrt(3)) / 15,
                (+sqrt(87) - 2 * sqrt(3)) / 15,
                (-sqrt(87) - 2 * sqrt(3)) / 15,
                (-sqrt(87) - 2 * sqrt(3)) / 15,
            ],
        ]
    }
    return C2Scheme("Schmid 4", d, 4, source, 4.441e-16)


def schmid_6():
    return _read(this_dir / "schmid_6.json", source)


register([schmid_2, schmid_4, schmid_6])

import numpy as np
import sympy

from ..helpers import article
from ._albrecht_collatz import albrecht_collatz as hammer_stroud_11_2
from ._helpers import S2Scheme, register
from ._peirce_1956 import peirce_1956_1, peirce_1956_3
from ._radon import radon

_source = article(
    authors=["Preston C. Hammer", "Arthur H. Stroud"],
    title="Numerical Evaluation of Multiple Integrals II",
    journal="Math. Comp.",
    volume="12",
    number="64",
    year="1958",
    month="oct",
    pages="272-280",
    url="https://doi.org/10.1090/S0025-5718-1958-0102176-6",
)

frac = sympy.Rational
pi = sympy.pi
sqrt = np.vectorize(sympy.sqrt)
pm_ = np.array([+1, -1])
cos = np.vectorize(sympy.cos)
sin = np.vectorize(sympy.sin)


def hammer_stroud_12_2():
    d = {
        "zero2": [[frac(1, 6)]],
        "d4_a0": [[frac(1, 6)], [sqrt(frac(1, 2))]],
        "sxy": [[frac(1, 24)], [sqrt(frac(1, 2))], [sqrt(frac(1, 2))]],
    }
    return S2Scheme("Hammer-Stroud 12-2", d, 5, _source)


def hammer_stroud_13_2():
    return peirce_1956_1()


def hammer_stroud_17():
    # DUP ENH This is Radon's formula.
    return radon(0)


def hammer_stroud_18():
    # ENH The article only gives floats, but really this is the spherical-product gauss
    # formula as described in Strouds book, S2 7-2.
    #
    # data = [
    #     (frac(1, 16), fs([0.4247082002778669, 0.1759198966061612])),
    #     (frac(1, 16), fs([0.8204732385702833, 0.3398511429799874])),
    # ]
    r1, r2 = sqrt((3 - pm_ * sqrt(3)) / 6)

    d = {"d8.0": [[frac(1, 16), frac(1, 16)], [r1, r2]]}
    return S2Scheme("Hammer-Stroud 18", d, 7, _source)


def hammer_stroud_19():
    sqrt6 = sqrt(6)
    alpha1 = (16 + sqrt6) / 288
    alpha2 = (137 - 32 * sqrt6) / 1818
    alpha3 = (520 + 155 * sqrt6) / 3636 / 8

    a = sqrt((6 + sqrt6) / 10)

    d = {
        "zero2": [[frac(1, 9)]],
        "d4_ab": [
            [alpha1, alpha3],
            [0.5505043204538557, 0.7932084745126058],
            [0.2280263556769715, 0.4645097310495256],
        ],
        "d4_a0": [[alpha2], [a]],
    }
    return S2Scheme("Hammer-Stroud 19", d, 9, _source)


def hammer_stroud_20():
    # ENH Also Peirce's formula, even given symbolically.
    return peirce_1956_3()


def hammer_stroud_21():
    alpha0 = 0.0341505695624825 / np.pi
    alpha1 = 0.0640242008621985 / np.pi
    alpha2 = 0.0341505695624825 / np.pi

    d = {
        "d4_ab": [
            [alpha0, alpha1, alpha1, alpha1, alpha1, alpha2, alpha2, alpha2],
            [
                0.2584361661674054,
                0.5634263397544869,
                0.4776497869993547,
                0.8028016728473508,
                0.6805823955716280,
                0.2190916025980981,
                0.9461239423417719,
                0.8020851487551318,
            ],
            [
                0.0514061496288813,
                0.1120724670846205,
                0.3191553840796721,
                0.1596871812824163,
                0.4547506180649039,
                0.1463923286035535,
                0.1881957532057769,
                0.5359361621905023,
            ],
        ],
    }
    return S2Scheme("Hammer-Stroud 21", d, 15, _source)


register(
    [
        hammer_stroud_11_2,
        hammer_stroud_12_2,
        hammer_stroud_13_2,
        hammer_stroud_17,
        hammer_stroud_18,
        hammer_stroud_19,
        hammer_stroud_20,
        hammer_stroud_21,
    ]
)

import pathlib

from sympy import Rational as frac
from sympy import sqrt

from ...helpers import article
from .._helpers import C2Scheme, concat, symm_r0, symm_s, symm_s_t, zero, _read

source = article(
    authors=["D.A. Dunavant"],
    title="Economical symmetrical quadrature rules for complete polynomials over a square domain",
    journal="Numerical Methods in Engineering",
    volume="21",
    number="10",
    month="oct",
    year="1985",
    pages="1777–1784",
    url="https://doi.org/10.1002/nme.1620211004",
)

this_dir = pathlib.Path(__file__).resolve().parent


def dunavant_00():
    weights, points = zero(1)
    return C2Scheme("Dunavant 0", weights, points, 1, source)


def dunavant_01():
    weights, points = symm_s([1, sqrt(frac(1, 3))])
    weights /= 4
    return C2Scheme("Dunavant 1", weights, points, 3, source)


def dunavant_02():
    weights, points = concat(
        symm_r0([frac(40, 49), sqrt(frac(7, 15))]),
        symm_s([frac(9, 49), sqrt(frac(7, 9))]),
    )
    weights /= 4
    return C2Scheme("Dunavant 2", weights, points, 5, source)


def dunavant_03():
    weights, points = concat(
        symm_r0([frac(98, 405), sqrt(frac(6, 7))]),
        symm_s(
            [0.237431774690630, 0.805979782918599],
            [0.520592916667394, 0.380554433208316],
        ),
    )
    weights /= 4
    return C2Scheme("Dunavant 3", weights, points, 7, source)


def dunavant_04():
    return _read(this_dir / "dunavant_04.json", source)


def dunavant_05():
    return _read(this_dir / "dunavant_05.json", source)


def dunavant_06():
    return _read(this_dir / "dunavant_06.json", source)


def dunavant_07():
    return _read(this_dir / "dunavant_07.json", source)


def dunavant_08():
    return _read(this_dir / "dunavant_08.json", source)


def dunavant_09():
    return _read(this_dir / "dunavant_09.json", source)


def dunavant_10():
    return _read(this_dir / "dunavant_10.json", source)
    weights, points = concat(
        symm_r0(
            [0.019503841092684, 0.980883148832881],
            [0.089012127744268, 0.678152700336576],
            [0.114568584702749, 0.240599282275864],
        ),
        symm_s(
            [0.007463627359106, 0.965176994929162],
            [0.050585943594705, 0.749698539312765],
            [0.074613865184212, 0.568983925500818],
        ),
        symm_s_t(
            [0.023501091310143, 0.971086142843168, 0.355832132274584],
            [0.011588562644144, 0.983453947854968, 0.645588139196562],
            [0.023073245798171, 0.933927707027213, 0.821920249234369],
            [0.001570221774472, 1.014086498915039, 0.862185099566557],
            [0.049102258016277, 0.877914842155496, 0.168914072450263],
            [0.042512352239126, 0.882246882640128, 0.568113580166780],
            [0.067270936863160, 0.741324453314596, 0.371360260002223],
            [0.103507336515645, 0.469570217710647, 0.237333359193547],
        ),
    )
    weights /= 4
    return C2Scheme("Dunavant 10", weights, points, 21, source, 2.246e-14)

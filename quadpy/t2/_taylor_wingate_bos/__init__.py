import pathlib

from sympy import Rational as frac

from ...helpers import article
from .._helpers import T2Scheme, expand_symmetries, _read

source = article(
    authors=["Mark A. Taylor", "Beth A. Wingate", "Len P. Bos"],
    title="Several new quadrature formulas for polynomial integration in the triangle",
    journal="arXiv Mathematics e-prints",
    year="2005",
    month="jan",
    url="https://arxiv.org/abs/math/0501496",
)

this_dir = pathlib.Path(__file__).resolve().parent


# TODO missing Taylor-Wingate-Bos schemes


def taylor_wingate_bos_1():
    d = {
        "s2": [[frac(2, 3)], [frac(1, 6)]]
    }
    points, weights = expand_symmetries(d)
    weights /= 2
    return T2Scheme("Taylor-Wingate-Bos 1", weights, points, 2, source)


def taylor_wingate_bos_2():
    return _read(this_dir / "taylor_wingate_bos_2.json", source)


def taylor_wingate_bos_4():
    return _read(this_dir / "taylor_wingate_bos_4.json", source)


def taylor_wingate_bos_5():
    return _read(this_dir / "taylor_wingate_bos_5.json", source)


def taylor_wingate_bos_8():
    return _read(this_dir / "taylor_wingate_bos_8.json", source)

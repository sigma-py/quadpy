import pathlib

from sympy import Rational as frac

from ...helpers import article
from .._helpers import T2Scheme, _read, register

source = article(
    authors=["D.M. Williams", "L. Shunn", "A. Jameson"],
    title="Symmetric quadrature rules for simplexes based on sphere close packed lattice arrangements",
    journal="Journal of Computational and Applied Mathematics",
    volume="266",
    year="2014",
    pages="18–38",
    url="https://doi.org/10.1016/j.cam.2014.01.007",
)

this_dir = pathlib.Path(__file__).resolve().parent


def williams_shunn_jameson_1():
    d = {"centroid": [[1]]}
    return T2Scheme("Williams-Shunn-Jameson 1", d, 1, source)


def williams_shunn_jameson_2():
    d = {"d3_aa": [[frac(1, 3)], [frac(1, 6)]]}
    return T2Scheme("Williams-Shunn-Jameson 2", d, 2, source)


def williams_shunn_jameson_3():
    return _read(this_dir / "williams_shunn_jameson_3.json", source)


def williams_shunn_jameson_4():
    return _read(this_dir / "williams_shunn_jameson_4.json", source)


def williams_shunn_jameson_5():
    return _read(this_dir / "williams_shunn_jameson_5.json", source)


def williams_shunn_jameson_6():
    return _read(this_dir / "williams_shunn_jameson_6.json", source)


def williams_shunn_jameson_7():
    return _read(this_dir / "williams_shunn_jameson_7.json", source)


def williams_shunn_jameson_8():
    return _read(this_dir / "williams_shunn_jameson_8.json", source)


register(
    [
        williams_shunn_jameson_1,
        williams_shunn_jameson_2,
        williams_shunn_jameson_3,
        williams_shunn_jameson_4,
        williams_shunn_jameson_5,
        williams_shunn_jameson_6,
        williams_shunn_jameson_7,
        williams_shunn_jameson_8,
    ]
)

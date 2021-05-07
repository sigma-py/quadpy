import pathlib

from sympy import Rational as frac
from sympy import sqrt

from ...helpers import article
from .._helpers import T3Scheme, _read, register

source = article(
    authors=["Lee Shunn", "Frank Ham"],
    title="Symmetric quadrature rules for tetrahedra based on a cubic close-packed lattice arrangement",
    journal="Journal of Computational and Applied Mathematics",
    year="2012",
    url="https://doi.org/10.1016/j.cam.2012.03.032",
)

this_dir = pathlib.Path(__file__).resolve().parent


def shunn_ham_1():
    d = {"s4": [[1]]}
    return T3Scheme("Shunn-Ham 2", d, 1, source)


def shunn_ham_2():
    d = {"s31": [[frac(1, 4)], [frac(1, 4) - sqrt(5) / 20]]}
    return T3Scheme("Shunn-Ham 2", d, 2, source)


def shunn_ham_3():
    return _read(this_dir / "shunn_ham_3.json", source)


def shunn_ham_4():
    return _read(this_dir / "shunn_ham_4.json", source)


def shunn_ham_5():
    return _read(this_dir / "shunn_ham_5.json", source)


def shunn_ham_6():
    return _read(this_dir / "shunn_ham_6.json", source)


register([shunn_ham_1, shunn_ham_2, shunn_ham_3, shunn_ham_4, shunn_ham_5, shunn_ham_6])

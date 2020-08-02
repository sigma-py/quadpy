import math
import pathlib

from ...helpers import article
from .._helpers import T3Scheme, _read, s4, s31

source = article(
    authors=["Lee Shunn", "Frank Ham"],
    title="Symmetric quadrature rules for tetrahedra based on a cubic close-packed lattice arrangement",
    journal="Journal of Computational and Applied Mathematics",
    year="2012",
    url="https://doi.org/10.1016/j.cam.2012.03.032",
)

this_dir = pathlib.Path(__file__).resolve().parent


def shunn_ham_1():
    weights, points = s4(1)
    return T3Scheme("Shunn-Ham 2", weights, points, 1, source)


def shunn_ham_2():
    weights, points = s31([1 / 4, 1 / (5 + math.sqrt(5))])
    return T3Scheme("Shunn-Ham 2", weights, points, 2, source)


def shunn_ham_3():
    return _read(this_dir / "shunn_ham_3.json", source)


def shunn_ham_4():
    return _read(this_dir / "shunn_ham_4.json", source)


def shunn_ham_5():
    return _read(this_dir / "shunn_ham_5.json", source)


def shunn_ham_6():
    return _read(this_dir / "shunn_ham_6.json", source)

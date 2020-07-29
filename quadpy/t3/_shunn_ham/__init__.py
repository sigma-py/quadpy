import math
import pathlib

from ...helpers import article
from .._helpers import T3Scheme, _read, concat, s4, s31, s211

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
    degree = 8
    weights, points = concat(
        s31(
            [0.0010373112336140, 0.0149520651530592],
            [0.0366291366405108, 0.1344783347929940],
        ),
        s211(
            [0.0096016645399480, 0.0340960211962615, 0.1518319491659370],
            [0.0164493976798232, 0.0462051504150017, 0.5526556431060170],
            [0.0153747766513310, 0.2281904610687610, 0.0055147549744775],
            [0.0293520118375230, 0.3523052600879940, 0.0992057202494530],
        ),
    )
    return T3Scheme("Shunn-Ham 6", weights, points, degree, source)

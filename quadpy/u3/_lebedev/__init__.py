import json
import os
import re

from ...helpers import article
from .._helpers import U3Scheme, cartesian_to_spherical, untangle2

# Sphere integration schemes from a series of publications, in chronological order
# <https://en.wikipedia.org/wiki/Lebedev_quadrature>
# <https://people.sc.fsu.edu/~jburkardt/datasets/sphere_lebedev_rule/sphere_lebedev_rule.html>
sources = [
    article(
        authors="V.I. Lebedev",
        title="Values of the nodes and weights of ninth to seventeenth order Gauss-Markov quadrature formulae invariant under the octahedron group with inversion",
        journal="Computational Mathematics and Mathematical Physics",
        volume="15",
        year="1975",
        pages="44-51",
        url="https://doi.org/10.1016/0041-5553(75)90133-0",
    ),
    article(
        authors=["V.I. Lebedev"],
        year="1976",
        title="Quadratures on a sphere",
        journal="Zh. Vȳchisl. Mat. Mat. Fiz.",
        volume="16",
        number="2",
        pages="293–306",
        url="https://doi.org/10.1016/0041-5553(76)90100-2",
    ),
    article(
        authors=["V.I. Lebedev"],
        title="Spherical quadrature formulas exact to orders 25-29",
        journal="Siberian Mathematical Journal",
        volume="18",
        year="1977",
        pages="99-107",
        url="https://doi.org/10.1007/BF00966954",
    ),
    article(
        authors=["V.I. Lebedev", "A.L. Skorokhodov"],
        title="Quadrature formulas of orders 41, 47, and 53 for the sphere",
        journal="Russian Acad. Sci. Dokl. Math.",
        volume="45",
        year="1992",
        pages="587-592",
    ),
    article(
        authors=["V.I. Lebedev"],
        title="A quadrature formula for the sphere of 59th algebraic order of accuracy",
        journal="Russian Acad. Sci. Dokl. Math.",
        volume="50",
        year="1995",
        pages="283-286",
    ),
    article(
        authors=["V.I. Lebedev", "D.N. Laikov"],
        title="A quadrature formula for the sphere of the 131st algebraic order of accuracy",
        journal="Doklady Mathematics",
        volume="59",
        number="3",
        year="1999",
        pages="477-481",
    ),
]


def _read(index, source):
    name = f"Lebedev({index})"

    this_dir = os.path.dirname(os.path.realpath(__file__))

    m = re.match("([0-9]+)([a-z]*)", index)
    filename = "lebedev_{:03d}{}.json".format(int(m.group(1)), m.group(2))
    with open(os.path.join(this_dir, filename), "r") as f:
        data = json.load(f)

    degree = data.pop("degree")

    points, weights = untangle2(data)
    azimuthal_polar = cartesian_to_spherical(points)
    return U3Scheme(name, weights, points, azimuthal_polar, degree, source)


def lebedev_003a():
    return _read("003a", sources[0])


def lebedev_003b():
    return _read("003b", sources[0])


def lebedev_003c():
    return _read("003c", sources[0])


def lebedev_005():
    return _read("005", sources[0])


def lebedev_007():
    return _read("007", sources[0])


def lebedev_009():
    return _read("009", sources[0])


def lebedev_011():
    return _read("011", sources[0])


def lebedev_013():
    return _read("013", sources[0])


def lebedev_015():
    return _read("015", sources[0])


def lebedev_017():
    return _read("017", sources[0])


def lebedev_019():
    return _read("019", sources[1])


def lebedev_021():
    return _read("021", sources[1])


def lebedev_023():
    return _read("023", sources[1])


def lebedev_025():
    return _read("025", sources[2])


def lebedev_027():
    return _read("027", sources[2])


def lebedev_029():
    return _read("029", sources[2])


def lebedev_031():
    return _read("031", sources[3])


def lebedev_035():
    return _read("035", sources[3])


def lebedev_041():
    return _read("041", sources[3])


def lebedev_047():
    return _read("047", sources[3])


def lebedev_053():
    return _read("053", sources[3])


def lebedev_059():
    return _read("059", sources[4])


def lebedev_065():
    return _read("065", sources[5])


def lebedev_071():
    return _read("071", sources[5])


def lebedev_077():
    return _read("077", sources[5])


def lebedev_083():
    return _read("083", sources[5])


def lebedev_089():
    return _read("089", sources[5])


def lebedev_095():
    return _read("095", sources[5])


def lebedev_101():
    return _read("101", sources[5])


def lebedev_107():
    return _read("107", sources[5])


def lebedev_113():
    return _read("113", sources[5])


def lebedev_119():
    return _read("119", sources[5])


def lebedev_125():
    return _read("125", sources[5])


def lebedev_131():
    return _read("131", sources[5])

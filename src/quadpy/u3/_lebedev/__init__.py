import pathlib

from ...helpers import article
from .._helpers import _read, register

# Sphere integration schemes from a series of publications, in chronological order
# <https://en.wikipedia.org/wiki/Lebedev_quadrature>
# <https://people.sc.fsu.edu/~jburkardt/datasets/sphere_lebedev_rule/sphere_lebedev_rule.html>
sources = [
    article(
        authors=["V.I. Lebedev"],
        title="Values of the nodes and weights of ninth to seventeenth order Gauss-Markov quadrature formulae invariant under the octahedron group with inversion",
        journal="Computational Mathematics and Mathematical Physics",
        volume="15",
        year="1975",
        pages="44-51",
        url="https://doi.org/10.1016/0041-5553%2875%2990133-0",
    ),
    article(
        authors=["V.I. Lebedev"],
        year="1976",
        title="Quadratures on a sphere",
        journal="Zh. Vȳchisl. Mat. Mat. Fiz.",
        volume="16",
        number="2",
        pages="293–306",
        url="https://doi.org/10.1016/0041-5553%2876%2990100-2",
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

this_dir = pathlib.Path(__file__).resolve().parent


def lebedev_003a():
    return _read(this_dir / "lebedev_003a.json", sources[0])


def lebedev_003b():
    return _read(this_dir / "lebedev_003b.json", sources[0])


def lebedev_003c():
    return _read(this_dir / "lebedev_003c.json", sources[0])


def lebedev_005():
    return _read(this_dir / "lebedev_005.json", sources[0])


def lebedev_007():
    return _read(this_dir / "lebedev_007.json", sources[0])


def lebedev_009():
    return _read(this_dir / "lebedev_009.json", sources[0])


def lebedev_011():
    return _read(this_dir / "lebedev_011.json", sources[0])


def lebedev_013():
    return _read(this_dir / "lebedev_013.json", sources[0])


def lebedev_015():
    return _read(this_dir / "lebedev_015.json", sources[0])


def lebedev_017():
    return _read(this_dir / "lebedev_017.json", sources[0])


def lebedev_019():
    return _read(this_dir / "lebedev_019.json", sources[1])


def lebedev_021():
    return _read(this_dir / "lebedev_021.json", sources[1])


def lebedev_023():
    return _read(this_dir / "lebedev_023.json", sources[1])


def lebedev_025():
    return _read(this_dir / "lebedev_025.json", sources[2])


def lebedev_027():
    return _read(this_dir / "lebedev_027.json", sources[2])


def lebedev_029():
    return _read(this_dir / "lebedev_029.json", sources[2])


def lebedev_031():
    return _read(this_dir / "lebedev_031.json", sources[3])


def lebedev_035():
    return _read(this_dir / "lebedev_035.json", sources[3])


def lebedev_041():
    return _read(this_dir / "lebedev_041.json", sources[3])


def lebedev_047():
    return _read(this_dir / "lebedev_047.json", sources[3])


def lebedev_053():
    return _read(this_dir / "lebedev_053.json", sources[3])


def lebedev_059():
    return _read(this_dir / "lebedev_059.json", sources[4])


def lebedev_065():
    return _read(this_dir / "lebedev_065.json", sources[5])


def lebedev_071():
    return _read(this_dir / "lebedev_071.json", sources[5])


def lebedev_077():
    return _read(this_dir / "lebedev_077.json", sources[5])


def lebedev_083():
    return _read(this_dir / "lebedev_083.json", sources[5])


def lebedev_089():
    return _read(this_dir / "lebedev_089.json", sources[5])


def lebedev_095():
    return _read(this_dir / "lebedev_095.json", sources[5])


def lebedev_101():
    return _read(this_dir / "lebedev_101.json", sources[5])


def lebedev_107():
    return _read(this_dir / "lebedev_107.json", sources[5])


def lebedev_113():
    return _read(this_dir / "lebedev_113.json", sources[5])


def lebedev_119():
    return _read(this_dir / "lebedev_119.json", sources[5])


def lebedev_125():
    return _read(this_dir / "lebedev_125.json", sources[5])


def lebedev_131():
    return _read(this_dir / "lebedev_131.json", sources[5])


register(
    [
        lebedev_003a,
        lebedev_003b,
        lebedev_003c,
        lebedev_005,
        lebedev_007,
        lebedev_009,
        lebedev_011,
        lebedev_013,
        lebedev_015,
        lebedev_017,
        lebedev_019,
        lebedev_021,
        lebedev_023,
        lebedev_025,
        lebedev_027,
        lebedev_029,
        lebedev_031,
        lebedev_035,
        lebedev_041,
        lebedev_047,
        lebedev_053,
        lebedev_059,
        lebedev_065,
        lebedev_071,
        lebedev_077,
        lebedev_083,
        lebedev_089,
        lebedev_095,
        lebedev_101,
        lebedev_107,
        lebedev_113,
        lebedev_119,
        lebedev_125,
        lebedev_131,
    ]
)

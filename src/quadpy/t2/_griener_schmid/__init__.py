import pathlib

from ...helpers import article
from .._helpers import _read, register

source = article(
    authors=["Bernhard Griener", "Hans Joachim Schmid"],
    title="An interactive tool to visualize common zeros of two-dimensional polynomials",
    journal="Journal of Computational and Applied Mathematics",
    volume="112",
    year="1999",
    pages="83-94",
    url="https://doi.org/10.1016/S0377-0427%2899%2900215-0",
)

# source2 = article(
#     authors="G.G. Rasputin",
#     title="Construction of cubature formulas containing prespecied knots",
#     journal="Metody Vychisl.",
#     volume="13",
#     year="1983",
#     pages="122–129",
#     note="in Russian."
# )

this_dir = pathlib.Path(__file__).resolve().parent


def griener_schmid_1():
    # According to the article, this scheme appeared earlier in `source2`.
    return _read(this_dir / "griener_schmid_1.json", source)


def griener_schmid_2():
    return _read(this_dir / "griener_schmid_2.json", source)


register([griener_schmid_1, griener_schmid_2])

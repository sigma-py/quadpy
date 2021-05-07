import pathlib

from ...helpers import article
from .._helpers import _read, register

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


def williams_shunn_jameson():
    return _read(this_dir / "williams_shunn_jameson.json", source)


register([williams_shunn_jameson])

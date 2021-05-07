import pathlib

from ...helpers import techreport
from .._helpers import _read, register

source = techreport(
    authors=["M. Beckers", "A. Haegemans"],
    title="The construction of cubature formulae for the tetrahedron",
    institution="Dept. of Computer Science, K.U. Leuven",
    year="1990",
    note="Report TW 128",
    url="https://lirias.kuleuven.be/handle/123456789/132648",
)

this_dir = pathlib.Path(__file__).resolve().parent


def beckers_haegemans_8():
    return _read(this_dir / "beckers_haegemans_8.json", source)


def beckers_haegemans_9():
    return _read(this_dir / "beckers_haegemans_9.json", source)


register([beckers_haegemans_8, beckers_haegemans_9])

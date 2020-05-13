# TODO There are more schemes in the technical report
from ..helpers import techreport, untangle
from ._helpers import S2Scheme, _s4, _s8, _s40

_citation = techreport(
    authors=["R. Cools", "A. Haegemans"],
    title="Construction of fully symmetric cubature formulae of degree 4k-3 for fully symmetric planar regions",
    year="1985",
    institution="Dept. of Computer Science, KU Leuven",
    number="Report TW 71",
    url="https://lirias.kuleuven.be/bitstream/123456789/131870/1/TW71.pdf",
)


def cools_haegemans_1():
    data = [
        (0.233253175473, _s4(0.459700843381)),
        (0.167468245269e-01, _s40(0.125592606040e01)),
    ]
    points, weights = untangle(data)
    return S2Scheme("Cools-Haegemans 1", weights, points, 5, _citation)


def cools_haegemans_2():
    data = [
        (0.567209601536e-01, _s8(0.243244191752, 0.809458260086)),
        (0.109948866164e00, _s4(0.302217386264)),
        (0.261900192462e-01, _s4(0.664341348594)),
        (0.419194282996e-03, _s40(0.134279080737e01)),
    ]
    points, weights = untangle(data)
    return S2Scheme("Cools-Haegemans 2", weights, points, 9, _citation)


def cools_haegemans_3():
    data = [
        (0.123447696401e-01, _s8(0.343855345294, 0.944778017142)),
        (0.932719633554e-01, _s4(0.277496500297)),
        (0.589496783783e-01, _s4(0.592355387396)),
        (0.730888189861e-01, _s40(0.778610819923)),
    ]
    points, weights = untangle(data)
    return S2Scheme("Cools-Haegemans 3", weights, points, 9, _citation)

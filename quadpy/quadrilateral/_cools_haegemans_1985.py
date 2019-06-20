# -*- coding: utf-8 -*-
#
# TODO There are three more schemes in the technical report
from ._helpers import QuadrilateralScheme, concat, symm_s_t, symm_r0, s4a
from ..helpers import techreport

citation = techreport(
    authors=["R. Cools", "A. Haegemans"],
    title="Construction of fully symmetric cubature formulae of degree 4k-3 for fully symmetric planar regions",
    year="1985",
    institution="Dept. of Computer Science, KU Leuven",
    number="Report TW 71",
    url="https://lirias.kuleuven.be/bitstream/123456789/131870/1/TW71.pdf",
)


def cools_haegemans_1985_1():
    weights, points = concat(
        symm_s_t([0.361130558151e-01, 0.344872025364, 0.918620441057]),
        s4a([0.535500902317e-01, 0.690880550486], [0.106828079664e-01, 0.939655258097]),
        symm_r0([0.113540990172e00, 0.488926856974]),
    )
    weights *= 4
    return QuadrilateralScheme("Cools-Haegemans 1985-1", weights, points, 9, citation)


def cools_haegemans_1985_2():
    weights, points = concat(
        symm_s_t(
            [0.348818790231e-01, 0.266676738695e-01, 0.377724312590],
            [0.344998496602e-01, 0.235988332487, 0.793396171109],
            [0.987441946914e-02, 0.265486560241, 0.978761747825],
            [0.203490805188e-01, 0.702141598362, 0.913909457030],
        ),
        s4a([0.475325029082e-01, 0.551473280570], [0.325703974952e-02, 0.968340720218]),
    )
    weights *= 4
    return QuadrilateralScheme("Cools-Haegemans 1985-2", weights, points, 13, citation)


def cools_haegemans_1985_3():
    weights, points = concat(
        symm_s_t(
            [0.197386321888e-01, 0.168234947696, 0.914794463441],
            [0.484363166325e-01, 0.252666976106, 0.591294378163],
            [0.421281899422e-02, 0.584047706043, 0.102139695463e01],
            [0.287255968895e-01, 0.586713014973, 0.826081709475],
        ),
        s4a([0.361061434781e-01, 0.178898689064], [0.116671271121e-01, 0.914197956909]),
    )
    weights *= 4
    return QuadrilateralScheme("Cools-Haegemans 1985-3", weights, points, 13, citation)

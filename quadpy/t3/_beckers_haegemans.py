from ..helpers import techreport
from ._helpers import T3Scheme, concat, s4, s22, s31, s211

citation = techreport(
    authors=["M. Beckers", "A. Haegemans"],
    title="The construction of cubature formulae for the tetrahedron",
    institution="Dept. of Computer Science, K.U. Leuven",
    year="1990",
    note="Report TW 128",
    url="https://lirias.kuleuven.be/handle/123456789/132648",
)


def beckers_haegemans_8():
    degree = 8
    weights, points = concat(
        s4(-0.020500188658639915),
        s31(
            [+0.014250305822866901, 0.20682993161067320],
            [+1.9670333131339009e-3, 0.082103588310546723],
            [+1.6983410909288737e-4, 5.7819505051979972e-3],
        ),
        s22([+4.5796838244672818e-3, 0.44946725998110577]),
        s211(
            [+5.7044858086819185e-3, 0.22906653611681113, 0.506227344977843697],
            [+2.1405191411620925e-3, 0.036607749553197423, 0.19048604193463345],
        ),
    )
    weights *= 6
    return T3Scheme("Beckers-Haegemans", weights, points, degree, citation)


def beckers_haegemans_9():
    degree = 9
    weights, points = concat(
        s4(-0.13779903832610864),
        s31(
            [+1.8653365690852895e-3, 0.048351038549736740],
            [+4.3094239694934006e-3, 0.32457928011788236],
            [-0.090184766481201525, 0.11461654022399521],
            [+0.044672576202511444, 0.22548995191151391],
        ),
        s211(
            [+0.034700405884550761, 0.13162780924686980, 0.083664701617184967],
            [+3.3525839026606469e-3, 0.43395146141140677, 0.10776985954942861],
            [+4.3162887555699692e-4, -1.3762773181382007e-3, 0.27655347263680734],
        ),
    )
    weights *= 6
    return T3Scheme("Beckers-Haegemans", weights, points, degree, citation)

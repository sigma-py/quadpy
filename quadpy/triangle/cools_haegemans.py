# -*- coding: utf-8 -*-
#
"""
R. Cools, A. Haegemans,
Construction of minimal cubature formulae for the square and the triangle using
invariant theory,
Department of Computer Science, K.U.Leuven,
TW Reports vol:TW96, Sept. 1987,
<https://lirias.kuleuven.be/handle/123456789/131869>.
"""
from mpmath import mp

from .helpers import rot_ab, TriangleScheme


def cools_haegemans_1(mpmath=False):
    flt = mp.mpf if mpmath else float

    mp.dps = 20
    weights, bary = rot_ab(
        [
            flt("0.16058343856681218798E-09"),
            flt("0.34579201116826902882E+00"),
            flt("0.36231682215692616667E+01"),
        ],
        [
            flt("0.26530624434780379347E-01"),
            flt("0.65101993458939166328E-01"),
            flt("0.87016510156356306078E+00"),
        ],
        [
            flt("0.29285717640155892159E-01"),
            flt("0.65177530364879570754E+00"),
            flt("0.31347788752373300717E+00"),
        ],
        [
            flt("0.43909556791220782402E-01"),
            flt("0.31325121067172530696E+00"),
            flt("0.63062143431895614010E+00"),
        ],
        [
            flt("0.66940767639916174192E-01"),
            flt("0.51334692063945414949E+00"),
            flt("0.28104124731511039057E+00"),
        ],
    )
    points = bary[:, 1:]
    weights *= 2
    return TriangleScheme("Cools-Haegemans 1", 8, weights, points, bary)


# TODO find error
# def cools_haegemans_2(mpmath=False):
#     flt = mp.mpf if mpmath else float
#
#     mp.dps = 20
#     data = [
#         (
#             0.15319130036758557631e-06,
#             _r3(+0.58469201683584513031e-01, -0.54887778772527519316e00),
#         ),
#         (
#             0.13260526227928785221e-01,
#             _r3(0.50849285064031410705e-01, 0.90799059794957813439e00),
#         ),
#         (
#             0.15646439344539042136e-01,
#             _r3(0.51586732419949574487e00, 0.46312452842927062902e00),
#         ),
#         (
#             0.21704258224807323311e-01,
#             _r3(0.24311033191739048230e00, 0.72180595182371959467e-00),
#         ),
#         (
#             0.21797613600129922367e-01,
#             _r3(0.75397765920922660134e-00, 0.20647569839132397633e00),
#         ),
#         (
#             0.38587913508193459468e-01,
#             _r3(0.42209207910846960294e-00, 0.12689533413411127327e00),
#         ),
#         (
#             0.39699584282594413022e-01,
#             _r3(0.19823878346663354068e00, 0.62124412566393319745e00),
#         ),
#         (0.47910534861520060665e-01, numpy.array([[1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0]])),
#     ]
#
#     bary, weights = untangle2(data)
#     points = bary[:, 1:]
#     weights *= 2
#     return TriangleScheme("Cools-Haegemans 1", 10, weights, points, bary)


CoolsHaegemans = {
    1: cools_haegemans_1,
    # 2: cools_haegemans_2,
}

import pathlib

from ...helpers import techreport
from .._helpers import _read, register

source = techreport(
    authors=["R. Cools", "A. Haegemans"],
    title="Construction of minimal cubature formulae for the square and the triangle using invariant theory",
    institution="Department of Computer Science, K.U.Leuven",
    # TW Reports vol:TW96,
    month="sep",
    year="1987",
    url="https://lirias.kuleuven.be/handle/123456789/131869",
)

this_dir = pathlib.Path(__file__).resolve().parent


def cools_haegemans_1(mpmath=False):
    return _read(this_dir / "cools_haegemans_1.json", source)


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
#         (0.47910534861520060665e-01, np.array([[1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0]])),
#     ]
#
#     points, weights = untangle2(data)
#     weights *= 2
#     return T2Scheme("Cools-Haegemans 1", 10, weights, points)


register([cools_haegemans_1])

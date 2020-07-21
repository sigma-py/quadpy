# TODO give all weights in terms in scientific notation to avoid round-off error
import math
import pathlib

from ...helpers import article, pm, untangle
from .._helpers import S2Scheme, _read

_source = article(
    authors=["Zhongxuan Luo", "Zhaoliang Meng"],
    title="Cubature formulas over the n-sphere",
    journal="Journal of Computational and Applied Mathematics",
    year="2007",
    volume="202",
    pages="511-522",
    url="https://doi.org/10.1016/j.cam.2006.03.004",
)

this_dir = pathlib.Path(__file__).resolve().parent


def luo_meng_1():
    return _read(this_dir / "luo_meng_1.json", _source, weight_factor=1 / math.pi)


def luo_meng_2():
    return _read(this_dir / "luo_meng_2.json", _source, weight_factor=1 / math.pi)


def luo_meng_3():
    return _read(this_dir / "luo_meng_3.json", _source, weight_factor=1 / math.pi)


def luo_meng_4():
    return _read(this_dir / "luo_meng_4.json", _source, weight_factor=1 / math.pi)


def luo_meng_5():
    return _read(this_dir / "luo_meng_5.json", _source, weight_factor=1 / math.pi)


# TODO find mistake
# def luo_meng_6():
#     data = [
#         (0.12566370614359, [[0, 0]]),
#         (0.05840846471928, pm([0.09675812994165, 0.36110625670843])),
#         (0.05840846471928, pm([0.26434812694957, 0.26434812676097])),
#         (0.05840846471928, pm([0.36110625674316, 0.09675812981205])),
#         (0.06992018223251, pm([0.14347639334166, 0.62914552544285])),
#         (0.06996062112205, pm([0.40212812179518, 0.50468063296759])),
#         (0.07001072267677, pm([0.58130075071609, 0.28017673859189])),
#         (0.07003290482486, pm([0.64529804558133, 0])),
#         (0.05451221994169, pm([0.16359877475730, 0.83450130452851])),
#         (0.05497050073244, pm([0.46782321680688, 0.71013972158889])),
#         (0.05555796700074, pm([0.70405774447146, 0.47692733074510])),
#         (0.05593580988433, pm([0.83361584135279, 0.16805241863048])),
#         (0.02373604417773, pm([0.15912091676987, 0.95790205017621])),
#         (0.02450554997213, pm([0.46530163489622, 0.85228527644703])),
#         (0.02542637441947, pm([0.72779904913934, 0.64280972920247])),
#         (0.02606318434828, pm([0.90714987061331, 0.34637395417630])),
#         (0.02628242756629, pm([0.97102821992230, 0])),
#     ]
#     points, weights = untangle(data)
#     weights /= math.pi
#     return S2Scheme("Luo-Meng 6", weights, points, 17, _source)


# TODO find error
# def luo_meng_7():
#     data = [
#         (0.62831853071796, [[0, 0]]),
#         (0.15739501981117, pm([0.22313923447538, 0.83276696021106])),
#         (0.15739501981117, pm([0.60962772574453, 0.60962772573512])),
#         (0.15739501981117, pm([0.83276696021278, 0.22313923446896])),
#         (0.04163923078518, pm([0.35648911914249, 1.56350584123606])),
#         (0.04166872260747, pm([0.99921385294824, 1.25427536200095])),
#         (0.04170520866795, pm([1.44454064039602, 0.69637442937729])),
#         (0.04172134401374, pm([1.60363181798263, 0])),
#         (0.00248954249688, pm([0.45703259734803, 2.34995743719975])),
#         (0.00252313408619, pm([1.31014474018574, 2.00367150787068])),
#         (0.00256386275246, pm([1.97776433468844, 1.34893550183970])),
#         (0.00258875098551, pm([2.34618078213053, 0.47604042816812])),
#         (1.709615741554363e-5, pm([0.45304937170440, 3.27851194591105])),
#         (2.053886363962915e-5, pm([1.44579291938731, 2.97717603559694])),
#         (2.222665327035189e-5, pm([2.39757845514132, 2.28155908626652])),
#         (2.292272152266965e-5, pm([3.06875713971016, 1.23960636097153])),
#         (2.312500019706606e-5, pm([3.30966679783376, 0])),
#     ]
#     points, weights = untangle(data)
#     weights /= math.pi
#     return S2Scheme("Luo-Meng 6", weights, points, 17, _source)

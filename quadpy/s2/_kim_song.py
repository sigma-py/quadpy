import math
import warnings

import numpy

from ..helpers import article, fsd, untangle, z
from ._helpers import S2Scheme

_source = article(
    authors=["KyoungJoong Kim", "ManSuk Song"],
    title="Symmetric quadrature formulas over a unit disk",
    journal="Korean J. Comp. & Appl. Math.",
    year="1997",
    volume="4",
    pages="179-192",
    url="https://doi.org/10.1007/BF03011388",
)


def kim_song_1():
    data = [
        (1.0, z(2)),
    ]
    points, weights = untangle(data)
    points = points.astype(numpy.float)
    return S2Scheme("Kim-Song 1", weights, points, 1, _source)


def kim_song_2():
    data = [
        (0.25, fsd(2, (math.sqrt(0.5), 1))),
    ]
    points, weights = untangle(data)
    return S2Scheme("Kim-Song 2", weights, points, 3, _source)


def kim_song_3():
    data = [
        (0.25, fsd(2, (0.5, 2))),
    ]
    points, weights = untangle(data)
    return S2Scheme("Kim-Song 3", weights, points, 3, _source)


def kim_song_4():
    data = [
        (0.732786462492640, fsd(2, (0.650115167343736, 1))),
        (0.052611700904808, fsd(2, (0.888073833977115, 2))),
    ]
    points, weights = untangle(data)
    weights /= math.pi
    return S2Scheme("Kim-Song 4", weights, points, 5, _source, 1.047e-14)


def kim_song_5():
    data = [
        (1 / 4, z(2)),
        (3 / 32, fsd(2, (0.754344479484572, 1), (0.312459714103782, 1))),
    ]
    points, weights = untangle(data)
    return S2Scheme("Kim-Song 5", weights, points, 5, _source)


def kim_song_6():
    data = [
        (2 / 27 * math.pi, fsd(2, (0.866025403784439, 1))),
        (0.387077796006226, fsd(2, (0.322914992067400, 2))),
        (0.165609800458645, fsd(2, (0.644171310389465, 2))),
    ]
    points, weights = untangle(data)
    weights /= math.pi
    return S2Scheme("Kim-Song 6", weights, points, 7, _source)


def kim_song_7():
    data = [
        (0.071488826617391, fsd(2, (0.952458896434417, 1))),
        (0.327176874928167, fsd(2, (0.415187657878755, 1))),
        (0.005591341512851, fsd(2, (0.834794942216211, 2))),
        (0.190570560169519, fsd(2, (0.740334457173511, 1), (0.379016937530835, 1))),
    ]
    points, weights = untangle(data)
    weights /= math.pi
    return S2Scheme("Kim-Song 7", weights, points, 9, _source, 7.364e-14)


# TODO find issue
def kim_song_8():
    data = [
        (-0.103220620849181, z(2)),
        (0.000561120485457, fsd(2, (1.081713912649978, 2))),
        (0.234185472691706, fsd(2, (0.329855087221954, 2))),
        (0.123731369524117, fsd(2, (0.163309595229403, 2))),
        (0.194337435615691, fsd(2, (0.582963679648472, 2))),
        (0.129193960146342, fsd(2, (0.862551409331466, 1), (0.191706749096887, 1))),
    ]
    points, weights = untangle(data)
    weights /= math.pi
    # ERR article claims degree 11
    warnings.warn("Kim-Song claim degree 11, but the scheme is only degree 1.")
    return S2Scheme("Kim-Song 8", weights, points, 1, _source, 2.001e-13)


def kim_song_9():
    data = [
        (0.051310052712355, fsd(2, (0.680167267076408, 2))),
        (0.208368275231940, fsd(2, (0.232463234651158, 2))),
        (0.113628206510048, fsd(2, (0.547722557505169, 2))),
        (0.126977836503225, fsd(2, (0.652159581445885, 1), (0.174745633184644, 1))),
        (0.079067977968328, fsd(2, (0.904823085572323, 1), (0.242446615072141, 1))),
    ]
    points, weights = untangle(data)
    weights /= math.pi
    return S2Scheme("Kim-Song 9", weights, points, 11, _source, 1.742e-14)


def kim_song_10():
    data = [
        (0.156239869837333, fsd(2, (0.283402832348825, 1))),
        (0.069447409266977, fsd(2, (0.920575474036251, 1))),
        (0.138145910534871, fsd(2, (0.649007230578002, 1))),
        (0.043376583887155, fsd(2, (0.677355106028069, 2))),
        (0.150228334140010, fsd(2, (0.412672216763289, 2))),
        (0.018199321926118, fsd(2, (0.938694583835123, 1), (0.335767600829539, 1))),
        (0.095780705939433, fsd(2, (0.752042776803954, 1), (0.379717011170077, 1))),
    ]
    points, weights = untangle(data)
    weights /= math.pi
    return S2Scheme("Kim-Song 10", weights, points, 13, _source, 3.073e-13)


def kim_song_11():
    data = [
        (0.125290208564338, fsd(2, (0.252863797091293, 1))),
        (0.016712625496982, fsd(2, (0.989746802511614, 1))),
        (0.109500391126365, fsd(2, (0.577728928444958, 1))),
        (0.066237455796397, fsd(2, (0.873836956645035, 1))),
        (0.026102860184358, fsd(2, (0.689299380791136, 2))),
        (0.066000934661100, fsd(2, (0.597614304667208, 2))),
        (0.127428372681720, fsd(2, (0.375416824626170, 2))),
        (0.042523065826681, fsd(2, (0.883097111318591, 1), (0.365790800400663, 1))),
        (0.081539591616413, fsd(2, (0.707438744960070, 1), (0.293030722710664, 1))),
    ]
    points, weights = untangle(data)
    weights /= math.pi
    return S2Scheme("Kim-Song 11", weights, points, 15, _source, 1.910e-12)


# TODO find issue
def kim_song_12():
    data = [
        (0.034218732123473, fsd(2, (0.000061843487605, 1))),
        (0.016529016194279, fsd(2, (0.980790725733799, 1))),
        (0.128441621643862, fsd(2, (0.536380801743208, 1))),
        (0.057523498415081, fsd(2, (0.859692027602019, 1))),
        (0.000019445945395, fsd(2, (0.890507932479849, 2))),
        (0.029882788585530, fsd(2, (0.597728327225623, 2))),
        (0.086625257876566, fsd(2, (0.472638508474430, 2))),
        (0.135289763196368, fsd(2, (0.256583368615187, 2))),
        (0.030442743336826, fsd(2, (0.915666628257231, 1), (0.290750921163427, 1))),
        (0.036355668706617, fsd(2, (0.775130682141990, 1), (0.550400039905684, 1))),
        (0.081635607665003, fsd(2, (0.722293685844502, 1), (0.274942296086914, 1))),
    ]
    points, weights = untangle(data)
    weights /= math.pi
    # ERR article claims degree 17
    warnings.warn("Kim-Song claim degree 17, but the scheme is only degree 3.")
    return S2Scheme("Kim-Song 12", weights, points, 3, _source)


def kim_song_13():
    data = [
        (0.107540143418461, fsd(2, (0.233779702361770, 1))),
        (0.021291625083247, fsd(2, (0.974967664200431, 1))),
        (0.082028519637405, fsd(2, (0.670603246539245, 1))),
        (0.000480959036869, fsd(2, (0.767149887105600, 2))),
        (0.055918006547614, fsd(2, (0.529910286773986, 2))),
        (0.020154425991986, fsd(2, (0.920197160917920, 1), (0.324880277222329, 1))),
        (0.056453982349758, fsd(2, (0.696315960030871, 1), (0.353373943382497, 1))),
        (0.038487694687371, fsd(2, (0.763756057536127, 1), (0.554819904211117, 1))),
        (0.096415834450115, fsd(2, (0.455714401760160, 1), (0.222878364006045, 1))),
        (0.047557517357696, fsd(2, (0.855907714936672, 1), (0.165103693017075, 1))),
    ]
    points, weights = untangle(data)
    weights /= math.pi
    return S2Scheme("Kim-Song 13", weights, points, 17, _source, 5.952e-12)


# TODO find issue
def kim_song_14():
    data = [
        (-0.336959088794964, fsd(2, (0.860125674956782, 1))),
        (-0.252102704609804, fsd(2, (0.860125119919931, 1))),
        (+0.131164048626801, fsd(2, (0.547071550894734, 1))),
        (-0.010919895825604, fsd(2, (0.927565701076356, 2))),
        (+0.010744805554708, fsd(2, (0.928162274936277, 2))),
        (+0.004798012482265, fsd(2, (0.755935467576794, 2))),
        (+0.072711702889623, fsd(2, (0.618205357502056, 2))),
        (+0.115883691032049, fsd(2, (0.398705530364909, 2))),
        (+0.123519504105529, fsd(2, (0.175541035301314, 2))),
        (+0.022426281825809, fsd(2, (0.959568314399471, 1), (0.155454100986327, 1))),
        (+0.033906878046312, fsd(2, (0.851375493057155, 1), (0.435023471372536, 1))),
        (+0.332072154179621, fsd(2, (0.858667011279674, 1), (0.030948726948493, 1))),
        (+0.074873729916680, fsd(2, (0.711569636216547, 1), (0.300178547400372, 1))),
    ]
    points, weights = untangle(data)
    weights /= math.pi
    # ERR article claims degree 19
    warnings.warn("Kim-Song claim degree 19, but the scheme is only degree 3.")
    return S2Scheme("Kim-Song 14", weights, points, 3, _source, 3.829e-14)


# TODO find issue
def kim_song_15():
    data = [
        (0.082558858859169, fsd(2, (0.204668989256100, 1))),
        (0.009721593541193, fsd(2, (0.992309839464756, 1))),
        (0.061920685878045, fsd(2, (0.740931035494388, 1))),
        (0.079123279187043, fsd(2, (0.477987648986077, 1))),
        (0.087526733002317, fsd(2, (0.306138805262459, 2))),
        (0.057076811471306, fsd(2, (0.524780156099700, 2))),
        (0.020981864256888, fsd(2, (0.921806074110042, 1), (0.310920075968188, 1))),
        (0.015226392255721, fsd(2, (0.790235832571934, 1), (0.579897645710646, 1))),
        (0.033136884897617, fsd(2, (0.725790566968788, 1), (0.525045580895713, 1))),
        (0.044853730819348, fsd(2, (0.788230650371813, 1), (0.290244481132460, 1))),
        (0.065321481701811, fsd(2, (0.584894890453686, 1), (0.264317463415838, 1))),
        (0.024214746797802, fsd(2, (0.909637445684200, 1), (0.092571132370888, 1))),
    ]
    points, weights = untangle(data)
    weights /= math.pi
    # ERR article claims degree 19
    warnings.warn("Kim-Song claim degree 19, but the scheme is only degree 9.")
    return S2Scheme("Kim-Song 15", weights, points, 9, _source, 6.563e-12)

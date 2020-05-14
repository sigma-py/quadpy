from ..helpers import article
from ._helpers import T3Scheme, concat, s4, s22, s31, s211

source = article(
    authors=["Lee Shunn", "Frank Ham"],
    title="Symmetric quadrature rules for tetrahedra based on a cubic close-packed lattice arrangement",
    journal="Journal of Computational and Applied Mathematics",
    year="2012",
    url="https://doi.org/10.1016/j.cam.2012.03.032",
)


def shunn_ham_1():
    degree = 1
    weights, points = concat(s4(1))
    return T3Scheme("Shunn-Ham 1", weights, points, degree, source)


def shunn_ham_2():
    degree = 2
    weights, points = s31([0.25, 0.1381966011250110])
    return T3Scheme("Shunn-Ham 2", weights, points, degree, source)


def shunn_ham_3():
    degree = 3
    weights, points = concat(
        s31([0.0476331348432089, 0.0738349017262234]),
        s22([0.1349112434378610, 0.0937556561159491]),
    )
    return T3Scheme("Shunn-Ham 3", weights, points, degree, source)


def shunn_ham_4():
    degree = 5
    weights, points = concat(
        s31(
            [0.0070670747944695, 0.0323525947272439],
            [0.1019369182898680, 0.3097693042728620],
        ),
        s211([0.0469986689718877, 0.0603604415251421, 0.2626825838877790]),
    )
    return T3Scheme("Shunn-Ham 4", weights, points, degree, source)


def shunn_ham_5():
    degree = 6
    weights, points = concat(
        s4(0.0931745731195340),
        s31([0.0021900463965388, 0.0267367755543735]),
        s22([0.0250305395686746, 0.0452454000155172]),
        s211(
            [0.0143395670177665, 0.0391022406356488, 0.7477598884818090],
            [0.0479839333057554, 0.2232010379623150, 0.0504792790607720],
        ),
    )
    return T3Scheme("Shunn-Ham 5", weights, points, degree, source)


def shunn_ham_6():
    degree = 8
    weights, points = concat(
        s31(
            [0.0010373112336140, 0.0149520651530592],
            [0.0366291366405108, 0.1344783347929940],
        ),
        s211(
            [0.0096016645399480, 0.0340960211962615, 0.1518319491659370],
            [0.0164493976798232, 0.0462051504150017, 0.5526556431060170],
            [0.0153747766513310, 0.2281904610687610, 0.0055147549744775],
            [0.0293520118375230, 0.3523052600879940, 0.0992057202494530],
        ),
    )
    return T3Scheme("Shunn-Ham 6", weights, points, degree, source)
